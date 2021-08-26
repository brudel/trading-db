extern "C" {
#include "postgres.h"
#include "fmgr.h" //Module magic
#include "funcapi.h" //Return row
#include "utils/array.h" //Array stuff
#include "utils/lsyscache.h" //get_typlenbyvalalign
#include "executor/spi.h" //SPI
#include <lapacke.h> // Eigenvectors and eigenvalues
// Commented because do not improve execution time
//#include <omp.h> // Parallelization

PG_MODULE_MAGIC;
}

#include <unordered_map>
#include "utils.h"

//# Retirar indicadores de erro antes da entrega final
//# Consulta ao bd significativamente mal otimizada para continent, pois
// poderia agrupar e não o faz.
//# Documentate
//# Retirar members vazios?

//#define PROFILE
#define COUNTRY_SIZE 3
#define PRODUCT_SIZE 6
#define REAL_VALUE "COALESCE((exp_val + imp_val) / 2, exp_val, imp_val, 0)"

#ifdef PROFILE
extern "C" {
#include <time.h>
}
clock_t init_t;
#define initick() init_t = clock();
#define tick(label) elog(INFO, "Tempo " label ": %lf ms", (double) (clock() - init_t) * 1000 / CLOCKS_PER_SEC);
#else
#define initick()
#define tick(label)
#endif


#define BG_FUNCTION_INFO_V1(func) extern "C" { \
PG_FUNCTION_INFO_V1(func); \
}

#define SBI_getText(tuple, tupdesc, col, is_null) DatumGetTextPP( \
		SPI_getbinval(tuple, tupdesc, col, is_null) \
	)

#define SBI_getString(tuple, tupdesc, col, is_null) ((char*) \
	VARDATA_ANY(SBI_getText(tuple, tupdesc, col, is_null)))

#define SBI_getDouble(tuple, tupdesc, col, is_null) DatumGetFloat8( \
		SPI_getbinval(tuple, tupdesc, col, is_null) \
	)

#define SBI_getInt(tuple, tupdesc, col, is_null) DatumGetInt32( \
		SPI_getbinval(tuple, tupdesc, col, is_null) \
	)

#define PMYSTRING_INIT(mstr, n) mstr = (mystring*) palloc(n + sizeof(mystring)); \
mstr->len = 0;

#define ARR_DIM(v) ArrayGetNItems(ARR_NDIM(v), ARR_DIMS(v))

#define HeapTupleHeaderHasNulls(hth) \
	(((hth)->t_infomask & HEAP_HASNULL) != 0)


typedef struct perm_mem {
	int n_indexes; // pci
	text** prods; // pci e eci_pci
	text** cntrs; // eci e eci_pci
	double* indexes; // eci e pci
	Datum vecs_tuple; // eci_pci
} perm_mem;

typedef std::unordered_map<text*, int, t_aux, t_aux> t_map;

enum index_t {
	ECI = 1,
	PCI,
	ECI_PCI
};


void create_common_countrs_map(t_map *countrs_map, ArrayType* groups, perm_mem* pm, index_t index)
{
	int count = 0;
	Datum daux;
	HeapTupleHeader hth;
	ArrayType* members;
	bool is_null;
	text* taux;

	if (index & ECI)
		pm->cntrs = (text**) palloc(sizeof(*pm->cntrs) * ARR_DIM(groups));

	//Retirar grupos com 0 elementos
	ArrayIterator itv = array_create_iterator(groups, 0, NULL);
	while (array_iterate(itv, &daux, &is_null))
	{
		hth = DatumGetHeapTupleHeader(daux);

		if (HeapTupleHeaderHasNulls(hth))
			elog(ERROR, "groups has null elements");

		members = DatumGetArrayTypeP(GetAttributeByNum(hth, 2, &is_null));

		if (ARR_HASNULL(members))
			elog(ERROR, "cgroup members has null elements");

		ArrayIterator itm = array_create_iterator(members, 0, NULL);

		// if (!ARR_NDIM(members))
		// 	continue;

		if (index & ECI)
			pm->cntrs[count] = DatumGetTextPP(GetAttributeByNum(hth, 1, &is_null));

		while (array_iterate(itm, &daux, &is_null))
		{
			taux = DatumGetTextPP(daux);
			if (VARSIZE_ANY_EXHDR(taux) != COUNTRY_SIZE)
				elog(ERROR, "country code text must be of size %d", COUNTRY_SIZE);

			if (!countrs_map->insert({taux, count}).second)
				elog(ERROR, "country belong to more than one group");
		}

		++count;
	}
}

void create_groups_countrs_map(t_map *countrs_map, ArrayType* groups, perm_mem* pm, index_t index)
{
	int count = 0;
	Datum daux;
	bool is_null;
	text* taux;

	if (index & ECI)
		pm->cntrs = (text**) palloc(sizeof(*pm->cntrs) * ARR_DIM(groups));

	ArrayIterator itv = array_create_iterator(groups, 0, NULL);
	while (array_iterate(itv, &daux, &is_null))
	{
		taux = DatumGetTextPP(daux);

		if (VARSIZE_ANY_EXHDR(taux) != COUNTRY_SIZE)
			elog(ERROR, "country code text must be of size %d", COUNTRY_SIZE);

		if (index & ECI)
			pm->cntrs[count] = taux;

		countrs_map->insert({taux, count++});
	}
}

int create_prods_map(t_map *prods_map, int hs_digits, perm_mem* pm, index_t index)
{
	int n;
	TupleDesc td;
	HeapTuple tuple;
	text* prod;
	bool is_null;

#define Q "SELECT DISTINCT left(hs_code, 0) FROM product"
	char* query = (char*) palloc(sizeof(Q));
	memcpy(query, Q, sizeof(Q));
	query[30] |= hs_digits;
#undef Q

	int status = SPI_execute(query, true, 0);
	pfree(query);

	if (status <= 0 || SPI_tuptable == NULL)
		elog(ERROR, "can't successfully access needed data on database");

	td = SPI_tuptable->tupdesc;
	prods_map->reserve(SPI_tuptable->numvals * 1.3);

	if (index & PCI)
		pm->prods = (text**) SPI_palloc(sizeof(*pm->prods) * SPI_tuptable->numvals);

	n = SPI_tuptable->numvals;
	for (int i = 0; i < n; i++)
	{
		tuple = SPI_tuptable->vals[i];

		if (index & PCI)
		{
			prod = (text*) SPI_palloc(hs_digits + VARHDRSZ);
			pm->prods[i] = prod;
		}
		else
			prod = (text*) palloc(hs_digits + VARHDRSZ);

		SET_VARSIZE(prod, hs_digits + VARHDRSZ);
		memcpy(VARDATA_ANY(prod), SBI_getString(tuple, td, 1, &is_null), hs_digits);

		prods_map->insert({prod, i});

	}
	SPI_freetuptable(SPI_tuptable);

	return n;
}

void calc_X_common_query(int s_year, int f_year, int hs_digits, bool c_groups)
{
	int status;
	mystring* query;

//GRUP BY make it a lot faster
#define QEXP "exporter"
#define QC_G "c_group"

#define QSELECT "SELECT "
#define QSCOLUMNS ", left(product, 0), sum(" REAL_VALUE ") FROM transaction "
#define QJOIN "JOIN country_group_belonging ON " \
	"(exporter = country AND year >= entry_year AND (exit_year IS NULL OR year < exit_year)) "
//# Possivelmente adicionar restrição c_group in na cláusula do join
#define QWHERE "WHERE year "
#define QWEQ "= "
#define QWGE ">= "
#define QWFYEAR " AND year <= "
#define QGROUP " GROUP BY "
#define QGPROD ", left(product, 0)\0"

	PMYSTRING_INIT(query, sizeof(QSELECT) + 2 * MAX(sizeof(QEXP) - 1,
		sizeof(QC_G) - 1) + sizeof(QSCOLUMNS) + sizeof(QJOIN) + sizeof(QWHERE)
		+ MAX(sizeof(QWEQ) + MAXINTSIZE - 1, sizeof(QWGE) + sizeof(QWFYEAR) +
		2 * MAXINTSIZE - 2) + sizeof(QGROUP) + sizeof(QGPROD) - 6);
	// -1 para cada string literal sem \0

	query->litcat(QSELECT);

	if (c_groups)
		query->litcat(QC_G);
	else
		query->litcat(QEXP);

	query->litcat(QSCOLUMNS);
	query->data[query->len - 81] |= hs_digits;

	if (c_groups)
		query->litcat(QJOIN);

	query->litcat(QWHERE);

	if (f_year == 0)
	{
		query->litcat(QWEQ);
		query->concat(s_year);
	}
	else
	{
		query->litcat(QWGE);
		query->concat(s_year);
		query->litcat(QWFYEAR);
		query->concat(f_year);
	}

	query->litcat(QGROUP);

	if (c_groups)
		query->litcat(QC_G);
	else
		query->litcat(QEXP);

	query->litcat(QGPROD);
	query->data[query->len - 3] |= hs_digits;

#undef QEXP
#undef QC_G
#undef QSELECT
#undef QSCOLUMNS
#undef QJOIN
#undef QWHERE
#undef QWEQ
#undef QWGE
#undef QWFYEAR
#undef QGROUP
#undef QGPROD

	status = SPI_execute(query->data, true, 0);
	pfree(query);

	if (status <= 0 || SPI_tuptable == NULL)
		elog(ERROR, "can't successfully access needed data on database");
}

int calc_X(ArrayType* groups, int s_year, int f_year, int hs_digits,
	double*** _X, double** _Xp, double** _Xc, perm_mem* pm, index_t index)
{
	double *Xc, *Xp;
	int c, p, n_prods, n_groups = ARR_DIM(groups);
	bool is_null, c_groups = groups->elemtype == TEXTOID;
	t_map* countrs_map;
	t_map* prods_map = new t_map(0, t_aux(hs_digits), t_aux(hs_digits));

	if (c_groups)
	{
		countrs_map = new t_map(241 * 1.3);
		create_groups_countrs_map(countrs_map, groups, pm, index);
	}
	else
	{
		countrs_map = new t_map(241 * 1.3, t_aux(COUNTRY_SIZE), t_aux(COUNTRY_SIZE));
		create_common_countrs_map(countrs_map, groups, pm, index);
	}
	tick("countrs_map");

	SPI_connect();
	n_prods = create_prods_map(prods_map, hs_digits, pm, index);
	tick("prods_map");

	double (*X)[n_prods] = (decltype(X)) SPI_palloc(sizeof(*X) * n_groups);
	*_X = (double**) X;
	*_Xc = Xc = (double*) SPI_palloc(sizeof(*Xc) * n_groups);
	*_Xp = Xp = (double*) SPI_palloc(sizeof(*Xp) * n_prods);


	//Primeira iteração do for em seguida, mas iniciando Xp
	for (int i = 0; i < n_prods; ++i)
	{
		X[0][i] = 0;
		Xp[i] = 0;
	}
	Xc[0] = 0;

	//Inicia X e Xc
	for (int i = 1; i < n_groups; ++i)
	{
		for (int j = 0; j < n_prods; ++j)
			X[i][j] = 0;

		Xc[i] = 0;
	}
	tick("calc_X set");

	calc_X_common_query(s_year, f_year, hs_digits, c_groups);
	tick("SPI_execute");

	TupleDesc tupdesc = SPI_tuptable->tupdesc;

	if (!SPI_tuptable->numvals)
		elog(ERROR, "no transactions were found with these parameters");

	for (int i = 0; i < SPI_tuptable->numvals; i++)
	{
		HeapTuple tuple = SPI_tuptable->vals[i];

		auto it = countrs_map->find(SBI_getText(tuple, tupdesc, 1, &is_null));
		if (it == countrs_map->end())
			continue;

		c = it->second;
		p = (*prods_map)[SBI_getText(tuple, tupdesc, 2, &is_null)];

		double value = SBI_getDouble(tuple, tupdesc, 3, &is_null);

		X[c][p] += value;
		//# Mais operações do que calcular em x[c][p] pronto para grupos, mas não precisa percorrer dnv
		Xc[c] += value;
		Xp[p] += value;
	}
	SPI_freetuptable(SPI_tuptable);

	SPI_finish();
	delete countrs_map, prods_map;

	return n_prods;
}

//# Se usar valor mínimo > 0 vai precisar ajustar Xp pra tirar os produtos excluídos
int filter_groups(double* Xc, double** _X, text** cntrs, int n_groups, int n_prods, index_t index)
{
	double (*X)[n_prods] = (decltype(X)) _X;

	while (--n_groups >= 0 && Xc[n_groups] == 0);

	if (n_groups < 0)
		elog(ERROR, "no transactions were found with non-zero value for these parameters");
	n_groups++;

	//Identifica
	for (int i = 0; i < n_groups; ++i)
		if (Xc[i] == 0)
		{
			while (Xc[--n_groups] == 0);
			if (n_groups < i) // It's never equal
			{
				n_groups = i;
				break;
			}

			Xc[i] = Xc[n_groups];

			for (int j = 0; j < n_prods; ++j)
				X[i][j] = X[n_groups][j];

			if (index & ECI)
				cntrs[i] = cntrs[n_groups];
		}

	if (n_groups < 2)
		elog(ERROR, "less than two groups have transactions with non-zero value");

	return n_groups;
}

//# Se usar valor mínimo > 0 vai precisar ajustar Xc pra tirar os produtos excluídos
int filter_products(double* Xp, double** _X, text** prods, int n_groups, int n_prods, index_t index)
{
	double (*X)[n_prods] = (decltype(X)) _X;

	for (int i = 0; i < n_prods; ++i)
		if (Xp[i] == 0)
		{
			while (Xp[--n_prods] == 0);
			if (n_prods < i) // It's never equal
			{
				n_prods = i;
				break;
			}

			Xp[i] = Xp[n_prods];

			for (int j = 0; j < n_groups; ++j)
				X[j][i] = X[j][n_prods];

			if (index & PCI)
			{
				SPI_pfree(prods[i]);
				prods[i] = prods[n_prods];
			}
		}

	if (n_prods < 2)
		elog(ERROR, "less than two products have transactions with non-zero value");

	return n_prods;
}

void calc_M(double** _X, double* Xp, double* Xc, int n_groups, int n_prods, int n_total_prods,
	double X_total, bool**_M, double* Mc, double* Mp)
{
	double (*X)[n_total_prods] = (decltype(X)) _X;
	bool (*M)[n_prods] = (decltype(M)) _M;

	Mc[0] = 0;
	for (int j = 0; j < n_prods; ++j)
	{
		M[0][j] = (X[0][j] * X_total / (Xc[0] * Xp[j]) >= 1 ? 1 : 0);
		Mc[0] += M[0][j];
		Mp[j] = M[0][j];
	}

	for (int i = 1; i < n_groups; ++i)
	{
		Mc[i] = 0;
		for (int j = 0; j < n_prods; ++j)
		{
			M[i][j] = (X[i][j] * X_total / (Xc[i] * Xp[j]) >= 1 ? 1 : 0);
			Mc[i] += M[i][j];
			Mp[j] += M[i][j];
		}
	}
}

void calc_W(bool**_M, double* Mc, double* Mp, int n_groups, int n_prods, double** _W)
{
	bool (*M)[n_prods] = (decltype(M)) _M;
	double (*W)[n_groups] = (decltype(W)) _W;

	// #pragma omp parallel for
	//Dava pra pular produtos que um país não tem especialidade
	for (int i = 0; i < n_groups; ++i)
		for (int j = 0; j < n_groups; ++j)
		{
			W[i][j] = 0;
			for (int p = 0; p < n_prods; ++p)
				//# Da pra fazer condicional (bool** M;) ou conversão
				W[i][j] += (M[i][p] & M[j][p]) / Mp[p];
			W[i][j] /= Mc[i];
		}
}

void calc_Kc(double** W, int n_groups, double* K)
{
	double *avlr, *avli, mean = 0, sum_quad = 0, stdev;
	int max[2];

	double (*avtr)[n_groups] = (decltype(avtr)) palloc(sizeof(*avtr) * n_groups);
	avlr = (double*) palloc(sizeof(*avlr) * n_groups);
	avli = (double*) palloc(sizeof(*avli) * n_groups);

	int info = LAPACKE_dgeev(LAPACK_ROW_MAJOR, 'N', 'V', n_groups, (double*)W,
		n_groups, avlr, avli, NULL, n_groups, (double*)avtr, n_groups);

	if (avlr[0] > avlr[1])
	{
		max[0] = 0;
		max[1] = 1;
	}
	else
	{
		max[0] = 1;
		max[1] = 0;
	}

	for (int i = 2; i < n_groups; ++i)
		if (avlr[i] > avlr[max[1]])
			if (avlr[i] > avlr[max[0]])
			{
				max[1] = max[0];
				max[0] = i;
			}
			else
				max[1] = i;

	for (int i = 0; i < n_groups; ++i)
		K[i] = avtr[i][max[1]];

	pfree(avlr);
	pfree(avli);
	pfree(avtr);
}

void calc_Kp(double* Kc, int n_groups, int n_prods, bool** _M, double* Mp, double* Kp)
{
	bool (*M)[n_prods] = (decltype(M)) _M;

	for (int j = 0; j < n_prods; ++j)
		Kp[j] = M[0][j] ? Kc[0] : 0;

	for (int i = 1; i < n_groups; ++i)
		for (int j = 0; j < n_prods; ++j)
			Kp[j] += M[i][j] ? Kc[i] : 0;

	for (int j = 0; j < n_prods; ++j)
		Kp[j] /= Mp[j];
}

void pack_indexes(double* Kc, int n_groups, double* Kp, int n_prods, perm_mem* pm,
	TupleDesc call_td)
{
	short elmlen;
	bool elmbyval;
	char elmalign;
	Datum data[2], *elems;
	bool isnull[] = {false, false};
	ArrayType* arr;
	TupleDesc td;

	elems = (Datum*) palloc(sizeof(*elems) * MAX(n_groups, n_prods));
	td = RelationNameGetTupleDesc("eciout");

	for (int i = 0; i < n_groups; ++i)
	{
		data[0] = PointerGetDatum(pm->cntrs[i]);
		data[1] = Float8GetDatumFast(Kc[i]);
		elems[i] = PointerGetDatum(heap_form_tuple(td, data, isnull)->t_data);
	}

	get_typlenbyvalalign(td->tdtypeid, &elmlen, &elmbyval, &elmalign);
	arr = construct_array(elems, n_groups, td->tdtypeid, elmlen, elmbyval, elmalign);

	td = RelationNameGetTupleDesc("pciout");

	for (int i = 0; i < n_prods; ++i)
	{
		data[0] = PointerGetDatum(pm->prods[i]);
		data[1] = Float8GetDatumFast(Kp[i]);
		elems[i] = PointerGetDatum(heap_form_tuple(td, data, isnull)->t_data);
	}

	get_typlenbyvalalign(td->tdtypeid, &elmlen, &elmbyval, &elmalign);
	data[1] = PointerGetDatum(construct_array(elems, n_prods, td->tdtypeid,
		elmlen, elmbyval, elmalign));

	data[0] = PointerGetDatum(arr);

	pm->vecs_tuple = HeapTupleGetDatum(heap_form_tuple(call_td, data, isnull));
}

void calc_indexes(FunctionCallInfo fcinfo, perm_mem* pm, index_t index, TupleDesc call_td)
{
	//# Mc e Mp podem virar int*, se fossem convertidos para double no uso
	double **X, *Xc, *Xp, X_total = 0, *Mc, *Mp, **W, *Kc, *Kp;
	bool** M;
	int n_total_prods, n_prods, n_groups = ARR_DIM(PG_GETARG_ARRAYTYPE_P(0));

	n_total_prods = calc_X(PG_GETARG_ARRAYTYPE_P(0), PG_GETARG_INT32(1), PG_GETARG_INT32(2),
		PG_GETARG_INT32(3) << 1, &X, &Xp, &Xc, pm, index);
	tick("calc_X");
	n_groups = filter_groups(Xc, X, pm->cntrs, n_groups, n_total_prods, index);
	tick("filter_groups");
	n_prods = filter_products(Xp, X, pm->prods, n_groups, n_total_prods, index);
	tick("filter_products");

	for (int i = 0; i < n_groups; ++i)
		X_total += Xc[i];

	M = (bool**) palloc(sizeof(*M) * n_groups * n_prods);
	Mc = (double*) palloc(sizeof(*Mc) * n_groups);
	Mp = (double*) palloc(sizeof(*Mp) * n_prods);

	calc_M(X, Xp, Xc, n_groups, n_prods, n_total_prods, X_total, M, Mc, Mp);
	pfree(X);
	pfree(Xc);
	pfree(Xp);
	tick("calc_M");

	W = (double**) palloc(sizeof(*W) * n_groups * n_prods);

	calc_W(M, Mc, Mp, n_groups, n_prods, W);
	pfree(Mc);
	if (index == ECI)
	{
		pfree(M);
		pfree(Mp);
	}
	tick("calc_W");

	Kc = (double*) palloc(sizeof(*Kc) * n_groups);

	calc_Kc(W, n_groups, Kc);
	pfree(W);
	tick("calc_Kc");

	if (index == ECI)
	{
		z_transform(Kc, n_groups);
		pm->indexes = Kc;

		pm->n_indexes = n_groups;
		return;
	}

	Kp = (double*) palloc(sizeof(*Kp) * n_prods);

	calc_Kp(Kc, n_groups, n_prods, M, Mp, Kp);
	tick("calc_Kp");
	pfree(M);
	pfree(Mp);

	if (index == PCI)
	{
		pfree(Kc);

		z_transform(Kp, n_prods);
		pm->indexes = Kp;

		pm->n_indexes = n_prods;
		return;
	}

	//ECI_PCI
	z_transform(Kc, n_groups);
	z_transform(Kp, n_prods);
	pack_indexes(Kc, n_groups, Kp, n_prods, pm, call_td);
}

void common_index_init(FunctionCallInfo fcinfo, index_t index)
{
	FuncCallContext *funcctx;
	perm_mem* pm;
	TupleDesc td;
	ArrayType* groups = PG_GETARG_ARRAYTYPE_P(0);
	int hs_digits = PG_GETARG_INT32(3);
	MemoryContext original_context;

	//Validate args
	if (ARR_HASNULL(groups))
		elog(ERROR, "groups has null elements");

	if (ARR_DIM(groups) < 2)
		elog(ERROR, "groups must have at least 2 elements");

	if (hs_digits > 3 || hs_digits < 1)
		elog(ERROR, "hs_digit_pairs must be 1, 2 or 3");

	// Cria contexto de chamada
	funcctx = SRF_FIRSTCALL_INIT();
	original_context = MemoryContextSwitchTo(funcctx->multi_call_memory_ctx);

	// Identifica tipo de retorno
	if (get_call_result_type(fcinfo, NULL, &td) != TYPEFUNC_COMPOSITE)
		ereport(ERROR,
		(
			errcode(ERRCODE_FEATURE_NOT_SUPPORTED),
			errmsg("function returning record called in context that cannot accept type record")
		));
	funcctx->tuple_desc = BlessTupleDesc(td);

	// Cria memória multichamada
	pm = (perm_mem*) palloc(sizeof(perm_mem));
	funcctx->user_fctx = pm;

	/*		Calcula índices      */

	initick();
	calc_indexes(fcinfo, pm, index, funcctx->tuple_desc);
	tick("calc_indexes");

	MemoryContextSwitchTo(original_context);
}

Datum return_table(FunctionCallInfo fcinfo, index_t index)
{
	FuncCallContext *funcctx;
	HeapTuple ht;
	HeapTupleHeader hth;
	Datum dt[2], ret, daux;
	bool isnull[2] = {false, false}, isNullAux;
	perm_mem* pm;

	//Recupera contexto da chamada
	funcctx = SRF_PERCALL_SETUP();
	pm = (perm_mem*) funcctx->user_fctx;

	//Encerra a última chamada
	if (funcctx->call_cntr == pm->n_indexes)
		SRF_RETURN_DONE(funcctx);

	//Cria resultados da chamada
	if (index == ECI)
		dt[0] = PointerGetDatum(pm->cntrs[funcctx->call_cntr]);
	else
		dt[0] = PointerGetDatum(pm->prods[funcctx->call_cntr]);

	dt[1] = Float8GetDatumFast(pm->indexes[funcctx->call_cntr]);

	//Cria e retorna tupla
	ht = heap_form_tuple(funcctx->tuple_desc, dt, isnull);
	ret = HeapTupleGetDatum(ht);
	SRF_RETURN_NEXT(funcctx, ret);
}

BG_FUNCTION_INFO_V1(common_eci);

// Args: groups (cgroup[] ou text[] para groups_), start_year integer, end_year integer, hs_digit_pairs integer
Datum common_eci(PG_FUNCTION_ARGS)
{
	if (SRF_IS_FIRSTCALL())
		common_index_init(fcinfo, ECI);

	return return_table(fcinfo, ECI);
}

BG_FUNCTION_INFO_V1(common_pci);

Datum common_pci(PG_FUNCTION_ARGS)
{
	if (SRF_IS_FIRSTCALL())
		common_index_init(fcinfo, PCI);

	return return_table(fcinfo, PCI);
}

BG_FUNCTION_INFO_V1(common_eci_pci);

Datum common_eci_pci(PG_FUNCTION_ARGS)
{
	common_index_init(fcinfo, ECI_PCI);

	perm_mem* pm = (perm_mem*) SRF_PERCALL_SETUP()->user_fctx;

	PG_RETURN_DATUM(pm->vecs_tuple);
}

text* add_country(mystring* query, ArrayIterator itv)
{
	Datum daux;
	bool is_null;
	text* country;

	array_iterate(itv, &daux, &is_null);
	country = DatumGetTextPP(daux);

	if (VARSIZE_ANY_EXHDR(country) != COUNTRY_SIZE)
		elog(ERROR, "country code text must be of size %d", COUNTRY_SIZE);

	query->concat(VARDATA_ANY(country), COUNTRY_SIZE);

	return country;
}

t_map* add_contries(mystring* query, ArrayType* g1, ArrayType* g2)
{
	int n;
	Datum daux;
	bool is_null;
	text* country;
	ArrayIterator itv;
	t_map* countrs_map;

	if (ARR_HASNULL(g1) || ARR_HASNULL(g2) || !(ARR_NDIM(g1) && ARR_NDIM(g2)))
		elog(ERROR, "country array is empty or has null elements");
	
	countrs_map = new t_map(ARR_DIM(g1) + ARR_DIM(g2) * 1.3, t_aux(COUNTRY_SIZE),
		t_aux(COUNTRY_SIZE));

#define QSEP "', '"

	itv = array_create_iterator(g1, 0, NULL);
	n = ARR_DIM(g1);
	for (int i = 0; i < n; ++i)
	{
		(*countrs_map)[add_country(query, itv)] = 1;
		query->litcat(QSEP);
	}

	itv = array_create_iterator(g2, 0, NULL);
	n = ARR_DIM(g2) - 1;
	for (int i = 0; i < n; ++i)
	{
		(*countrs_map)[add_country(query, itv)] |= 2;
		query->litcat(QSEP);
	}

#undef QSEP

	(*countrs_map)[add_country(query, itv)] |= 2;

	return countrs_map;
}

t_map* query_series(ArrayType* g1, ArrayType* g2, int start_yi, int end_yi, VarChar* prod)
{
	int len = VARSIZE_ANY_EXHDR(prod);
	mystring* query;
	t_map* countrs_map;
	int status;

	if (len < 0 || len > 6 || len & 1)
		elog(ERROR, "hs_code must have 0, 2, 4 or 6 digits");

#define QSELECT "SELECT exporter, product, year, sum(" REAL_VALUE ")"\
	" FROM transaction WHERE exporter IN ('"
#define QFCLOSE "')"
#define QWSYEAR " AND year >= "
#define QWFYEAR " AND year <= "
#define QWPROD " AND left(product, 0) = '"
#define QGROUP " GROUP BY exporter, product, year ORDER BY year, product\0"

	len += sizeof(QSELECT) + sizeof(QFCLOSE) + sizeof(QWSYEAR) + sizeof(QWFYEAR) + sizeof(QWPROD) +
		+ sizeof(QGROUP) + 2 * MAXINTSIZE - 5;
	//5 = 1(\') - 6 (literals \0s)

	PMYSTRING_INIT(query, VARHDRSZ + len);
	query->litcat(QSELECT);

	countrs_map = add_contries(query, g1, g2);

	query->litcat(QFCLOSE);

	if (start_yi != 0)
	{
		query->litcat(QWSYEAR);
		query->concat(start_yi);
	}

	if (end_yi != 0)
	{
		query->litcat(QWFYEAR);
		query->concat(end_yi);
	}

	if (VARSIZE_ANY_EXHDR(prod) > 0)
	{
		query->litcat(QWPROD);
		query->data[query->len - 6] += VARSIZE_ANY_EXHDR(prod);
		query->concat(VARDATA_ANY(prod), VARSIZE_ANY_EXHDR(prod));
		query->concat('\'');
	}

	query->litcat(QGROUP);

#undef QSELECT
#undef QFCLOSE
#undef QWSYEAR
#undef QWFYEAR
#undef QWPROD
#undef QGROUP

	status = SPI_execute(query->data, true, 0);
	pfree(query);

	if (status <= 0 || SPI_tuptable == NULL)
		elog(ERROR, "can't successfully access needed data on database");

	return countrs_map;
}

BG_FUNCTION_INFO_V1(euclidean_distance);

//country_1 text, country_2 text, start_year integer, end_year integer, hs_code varchar(6)
Datum euclidean_distance(PG_FUNCTION_ARGS)
{
	bool is_null;
	t_map* countrs_map;
	TupleDesc td;
	t_aux p_compare(PRODUCT_SIZE);
	HeapTuple first, current;
	double dist = 0, aux, g1, g2;
	int mask;

	SPI_connect();

	countrs_map = query_series(PG_GETARG_ARRAYTYPE_P(0), PG_GETARG_ARRAYTYPE_P(1),
		PG_GETARG_INT32(2), PG_GETARG_INT32(3), PG_GETARG_VARCHAR_PP(4));

	td = SPI_tuptable->tupdesc;

	if (!SPI_tuptable->numvals)
	{
		SPI_freetuptable(SPI_tuptable);
		SPI_finish();
		PG_RETURN_FLOAT8(0.0);
	}

	for (int i = 0; i < SPI_tuptable->numvals;)
	{
		current = first = SPI_tuptable->vals[i];

		g1 = g2 = 0;

		do
		{
			aux = SBI_getDouble(current, td, 4, &is_null);
			mask = (*countrs_map)[SBI_getText(current, td, 1, &is_null)];

			if (mask & 1)
				g1 += aux;

			if (mask & 2)
				g2 += aux;

			current = SPI_tuptable->vals[i];
		}
		while (++i < SPI_tuptable->numvals
			&& p_compare(SBI_getText(current = SPI_tuptable->vals[i], td, 2, &is_null),
				SBI_getText(first, td, 2, &is_null))
			&& SBI_getInt(current, td, 3, &is_null) == SBI_getInt(first, td, 3, &is_null));

		aux = g1 - g2;

		dist += aux * aux;
	}

	SPI_freetuptable(SPI_tuptable);
	SPI_finish();

	PG_RETURN_FLOAT8(sqrt(dist));
}