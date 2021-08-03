extern "C" {
#include "postgres.h"
#include "fmgr.h" //Module magic
#include "funcapi.h" //Return row
#include "utils/array.h" //Coisas de array
#include "utils/lsyscache.h" //get_typlenbyvalalign
#include "executor/spi.h" //SPI
#include <lapacke.h> //AlgeLin
#include <omp.h> //Paralelização

PG_MODULE_MAGIC;
}

#include <unordered_map>
#include "utils.h"

//# Retirar indicadores de erro antes da entrega final

//#define PROFILE

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

typedef struct perm_mem {
	ArrayIterator ar_it; // eci e eci_pci
	int n_prods; // pci e eci_pci
	text** prods; // pci e eci_pci
	double* indexes; // eci e pci
	Datum vecs_tuple; // eci_pci
} perm_mem;

#define SBI_getString(tuple, tupdesc, col, is_null) (char*) \
	VARDATA_ANY( DatumGetTextP( \
		SPI_getbinval(tuple, tupdesc, col, is_null) \
	))

#define SBI_getDouble(tuple, tupdesc, col, is_null) DatumGetFloat8( \
		SPI_getbinval(tuple, tupdesc, col, is_null) \
	)

#define SBI_getInt(tuple, tupdesc, col, is_null) DatumGetInt32( \
		SPI_getbinval(tuple, tupdesc, col, is_null) \
	)

#define PMYSTRING_INIT(mstr, n) mstr = (mystring*) palloc(n); \
mstr->len = 0;

typedef std::unordered_map<char*, int, l_str, l_str> l_map;

enum index_t {
	ECI = 1,
	PCI,
	ECI_PCI
};

void create_countrs_map(l_map *countrs_map, ArrayType* groups)
{
	int count = 0;
	Datum daux;
	HeapTupleHeader hth;
	ArrayType* members;
	bool is_null;
	text* taux;

	ArrayIterator itv = array_create_iterator(groups, 0, NULL);
	while (array_iterate(itv, &daux, &is_null))
	{
		hth = DatumGetHeapTupleHeader(daux);
		members = DatumGetArrayTypeP(GetAttributeByNum(hth, 2, &is_null));

		ArrayIterator itm = array_create_iterator(members, 0, NULL);
		while (array_iterate(itm, &daux, &is_null))
		{
			taux = DatumGetTextP(daux);
			countrs_map->insert({{VARDATA_ANY(taux), count}});
		}

		++count;
	}
}

int create_prods_map(l_map *prods_map, int hs_digits, perm_mem* pm, index_t index)
{
	int n;
	bool is_null;

#define Q "SELECT DISTINCT left(hs_code, 0) FROM product"
	char* query = (char*) palloc(sizeof(Q));
	memcpy(query, Q, sizeof(Q));
	query[30] |= hs_digits;
#undef Q

	//elog(INFO, "query: %s", query);
	int status = SPI_execute(query, true, 0);
	pfree(query);

	if (status > 0 && SPI_tuptable != NULL)
	{
		TupleDesc tupdesc = SPI_tuptable->tupdesc;
		prods_map->reserve(SPI_tuptable->numvals * 1.3);

		if (index & PCI)
			pm->prods = (text**) SPI_palloc(sizeof(*pm->prods) * SPI_tuptable->numvals);

		//elog(INFO, "create prods: %d", SPI_tuptable->numvals);
		n = SPI_tuptable->numvals;
		for (int i = 0; i < n; i++)
		{
			HeapTuple tuple = SPI_tuptable->vals[i];
			char* prod = SBI_getString(tuple, tupdesc, 1, &is_null);

			//elog(INFO, "i = %d", i);
			prods_map->insert({{prod, i}});
			//elog(INFO, "%d", i);

			if (index & PCI)
			{//# Usar para o mapa também
				text* prod_t = (text*) SPI_palloc(hs_digits + VARHDRSZ);
				SET_VARSIZE(prod_t, hs_digits + VARHDRSZ);
				memcpy(VARDATA_ANY(prod_t), prod, hs_digits);

				pm->prods[i] = prod_t;
				//prods_map->insert({VARDATA(prod_t), i});
			}
		}
	}
	SPI_freetuptable(SPI_tuptable);


	//for (auto it = prods_map->begin(); it != prods_map->end(); ++it)
	//	elog(INFO, "prods[%.6s] = %d", it->first, it->second);

	//for (int i = 0; i < prods_map->bucket_count(); ++i)
	//	elog(INFO, "prods[%d] = %d", i, prods_map->bucket_size(i));

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
#define QSCOLUMNS ", left(product, 0),"\
	" sum(COALESCE((exp_val + imp_val) / 2, exp_val, imp_val)) FROM transaction "
#define QJOIN "JOIN country_group_belonging ON " \
	"(exporter = country AND year >= entry_year AND (exit_year IS NULL OR year < exit_year)) "
// Possivelmente adicionar restrição c_name in na cláusula do join
#define QWHERE "WHERE year "
#define QWEQ "= "
#define QWGE ">= "
#define QWFYEAR " AND year <= "
#define QGROUP " GROUP BY "
#define QGPROD ", left(product, 0)"

	PMYSTRING_INIT(query, sizeof(QSELECT) + 2 * MAX(sizeof(QEXP) - 1,
		sizeof(QC_G) - 1) + sizeof(QSCOLUMNS) + sizeof(QJOIN) + sizeof(QWHERE)
		+ MAX(sizeof(QWEQ) + MAXINTSIZE - 1, sizeof(QWGE) + sizeof(QWFYEAR) +
		2 * MAXINTSIZE - 2) + sizeof(QGROUP) + sizeof(QGPROD) - 5);
	// -1 para cada string literal sem \0

	query->litcat(QSELECT);

	if (c_groups)
		query->litcat(QC_G);
	else
		query->litcat(QEXP);

	query->litcat(QSCOLUMNS);
	query->data[query->len - 78] |= hs_digits;

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

	query->concat(QGPROD, sizeof(QGPROD));
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

	//elog(INFO, "query: %s", query->data);
	status = SPI_execute(query->data, true, 0);
	pfree(query);

	if (status <= 0 && SPI_tuptable == NULL)
		elog(ERROR, "falha ao consultar o bd");
}

int calc_X(ArrayType* groups, int s_year, int f_year, int hs_digits,
	double*** _X, double** _Xp, double** _Xc, perm_mem* pm, index_t index)
{
	double *Xc, *Xp;
	int c, p, n_prods, n_groups = ARR_DIMS(groups)[0];
	bool is_null, c_groups = groups->elemtype == TEXTOID;

	l_map* countrs_map = new l_map(241 * 1.3, l_str(3), l_str(3));
	l_map* prods_map = new l_map(0, l_str(hs_digits), l_str(hs_digits));
	create_countrs_map(countrs_map, groups);
	tick("countrs_map");
	n_prods = create_prods_map(prods_map, hs_digits, pm, index);
	tick("prods_map");

	double (*X)[n_prods] = (decltype(X)) palloc(sizeof(*X) * n_groups);
	*_X = (double**) X;
	*_Xc = Xc = (double*) palloc(sizeof(*Xc) * n_groups);
	*_Xp = Xp = (double*) palloc(sizeof(*Xp) * n_prods);


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
		//elog(INFO, "Bucket[%d] > %d", i, countrs_map->bucket_size(i));
		for (int j = 0; j < n_prods; ++j)
			X[i][j] = 0;

		Xc[i] = 0;
	}
	tick("calc_X set");

	//elog(INFO, "Total %d", countrs_map->bucket_count());

	calc_X_common_query(s_year, f_year, hs_digits, c_groups);
	tick("SPI_execute");

	TupleDesc tupdesc = SPI_tuptable->tupdesc;
/*
	auto coiso = countrs_map->hash_function();
	elog(INFO, "%ld, %ld, %d", coiso("gbr"), coiso("gbr"), coiso("gbr", "gbr"));

	auto itt = countrs_map->find("gbr");
	if (itt == countrs_map->end())
		elog(INFO, "erroww");
	else
		elog(INFO, "sec: %d", itt->second);*/


	//elog(INFO, "n: %d", SPI_tuptable->numvals);
	for (int i = 0; i < SPI_tuptable->numvals; i++)
	{
		//elog(INFO, "a");
		HeapTuple tuple = SPI_tuptable->vals[i];

		auto it = countrs_map->find(SBI_getString(tuple, tupdesc, 1, &is_null));
		//elog(INFO, "truple[%d]='%3s'", i, SBI_getString(tuple, tupdesc, 1, &is_null));
		if (it == countrs_map->end())
			continue;

		c = it->second;
		p = (*prods_map)[SBI_getString(tuple, tupdesc, 2, &is_null)];
		//elog(INFO, "%d %d", i, p);
		//elog(INFO, "b");

		double value = SBI_getDouble(tuple, tupdesc, 3, &is_null);

		X[c][p] += value;
		//# Mais operações do que calcular em x[c][p] pronto, mas não precisa percorrer dnv
		Xc[c] += value;
		Xp[p] += value;
		//elog(INFO, "c");
	}
	SPI_freetuptable(SPI_tuptable);

	return n_prods;
}

//# Não vai precisar mudar o mapa
//# Se usar valor mínimo > 0 vai precisar ajustar Xc pra tirar os produtos excluídos
int filter_products(double* Xp, double** _X, text** prods, int n_groups, int n_prods, index_t index)
{
	//elog(INFO, "Entro");
	int count = 0, idx, m_count, aux;
	double (*X)[n_prods] = (decltype(X)) _X;
	int* eliminated = (int*) palloc(sizeof(*eliminated)* n_prods), *moved;

	//Identifica
	for (int i = 0; i < n_prods; ++i)
		if (Xp[i] == 0)
			eliminated[count++] = i;

	if (count == 0)
		return n_prods;
	//elog(INFO, "Identifico");

		// Prepara variáveis para o loops
	m_count = 0;
	aux = n_prods - count; //Menor a ser movido

	idx = count;
	while (eliminated[--idx] >= aux); //Index do menor removido do intervalo dos movidos
	
	// Pula removidos entre os últimos
	// idx itera sobre os removidos entre os últios
	// aux itera sobre todos os últimos
	// m_count itera sobre os removidos
	for (; ++idx < count; ++aux)
		while (aux != eliminated[idx])
		{
			Xp[eliminated[m_count]] = Xp[aux];

			if (index & PCI)
			{
				SPI_pfree(prods[eliminated[m_count]]);
				prods[eliminated[m_count]] = prods[aux];
			}

			for (int i = 0; i < n_groups; ++i)
				X[i][eliminated[m_count]] = X[i][aux];

			m_count++;
			aux++;
		}

	// Nesse loop não existem mais removiso entre os últimos
	while (aux < n_prods)
	{
		Xp[eliminated[m_count]] = Xp[aux];

		if (index & PCI)
		{
			SPI_pfree(prods[eliminated[m_count]]);
			prods[eliminated[m_count]] = prods[aux];
		};

		for (int i = 0; i < n_groups; ++i)
			X[i][eliminated[m_count]] = X[i][aux];

		m_count++;
		aux++;
	}

	//Debug
	for (int i = 0; i < count; ++i)
		if (Xp[i] == 0)
			elog(INFO, "Erro na remoção de produtos não comercializados.");

	//elog(INFO, "NULLs: %d", count);
	pfree(eliminated);

	return n_prods - count;
}

void calc_M(double** _X, double* Xp, double* Xc, int n_groups, int n_prods, int n_total_prods,
	double X_total, char**_M, double* Mc, double* Mp)
{
	double (*X)[n_total_prods] = (decltype(X)) _X;
	char (*M)[n_prods] = (decltype(M)) _M;

	Mc[0] = 0;
	for (int j = 0; j < n_prods; ++j)
	{
		//elog(INFO, "M[%d][%d] = %lf * %lf / (%lf * %lf)", 0, j, X[0][j], X_total, Xc[0],  Xp[j]);
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
			//elog(INFO, "M[%d][%d] %d = %lf * %lf / (%lf * %lf)", i, j, M[i][j], X[i][j], X_total, Xc[i], Xp[j]);
			Mc[i] += M[i][j];
			Mp[j] += M[i][j];
		}
		//elog(INFO, "i: %d", i);
	}
}

void calc_W(char**_M, double* Mc, double* Mp, int n_groups, int n_prods, double** _W)
{
	char (*M)[n_prods] = (decltype(M)) _M;
	double (*W)[n_groups] = (decltype(W)) _W;

	#pragma omp parallel for
	//Dava pra pular produtos que um país não tem especialidade
	for (int i = 0; i < n_groups; ++i)
	{
		for (int j = 0; j < n_groups; ++j)
		{
			W[i][j] = 0;
			for (int p = 0; p < n_prods; ++p)
				//# Da pra fazer condicional (bool** M;) ou conversão
				W[i][j] += (M[i][p] & M[j][p]) / Mp[p];
			W[i][j] /= Mc[i];
		}
		//elog(INFO, "i: %d", i);
	}
}

void calc_Kc(double** W, int n_groups, double* K)
{
	double mean = 0, sum_quad = 0, stdev;
	double* avlr = (double*) palloc(sizeof(*avlr) * n_groups);
	double* avli = (double*) palloc(sizeof(*avli) * n_groups);
	double (*avtr)[n_groups] = (decltype(avtr)) palloc(sizeof(*avtr) * n_groups);

	int info = LAPACKE_dgeev(LAPACK_ROW_MAJOR, 'N', 'V', n_groups, (double*)W,
		n_groups, avlr, avli, NULL, n_groups, (double*)avtr, n_groups);
	//elog(INFO, "info %d", info);

	int max = 1;

	//# Se o 1 sempre for o primeiro
	for (int i = 2; i < n_groups; ++i)
		if (avlr[i] > avlr[max])
			max = i;
/*
	for (int i = 0; i < n_groups; ++i)
		elog(INFO, "Eigenvalue[%d]: %lf + %lfi", i, avlr[i], avli[i]);

	elog(INFO, "Maior autovalor: %d\n", max);

	for (int i = 0; i < n_groups; ++i)
		elog(INFO, "Eigenvvector[%d][%d]: %lf", i, 0, avtr[i][0]);
*/	

	for (int i = 0; i < n_groups; ++i)
		K[i] = avtr[i][max];

	pfree(avlr);
	pfree(avli);
	pfree(avtr);
}

void calc_Kp(double* Kc, int n_groups, int n_prods, char** _M, double* Mp, double* Kp)
{
	char (*M)[n_prods] = (decltype(M)) _M;

	for (int j = 0; j < n_prods; ++j)
		Kp[j] = M[0][j] ? Kc[0] : 0;

	for (int i = 1; i < n_groups; ++i)
		for (int j = 0; j < n_prods; ++j)
			Kp[j] += M[i][j] ? Kc[i] : 0;

	for (int j = 0; j < n_prods; ++j)
		Kp[j] /= Mp[j];
}

/*void teste()
{
	ArrayType* arr = palloc0(ARR_OVERHEAD_NONULLS(1) + n_groups * );
	arr->ndim = 1;
	arr->dataoffset = 0;
	arr->elemtype = td->Oid;
	ARR_DIMS(arr)[0] = n_groups;
	ARR_DIMS(arr)[1] = 1; // lower bound[0] hard coded
	char* data = ARR_DIMS(arr) + 2;

	len = VARSIZE_ANY(ht);
	memcpy(data, ht, len);
	data += len;
	data = att_align_nominal(data, MAXIMUM_ALIGNOF);
}*/

void pack_indexes(double* Kc, int n_groups, double* Kp, int n_prods,
		perm_mem* pm, TupleDesc call_td)
{
	elog(INFO, "pack");
	short elmlen;
	bool elmbyval;
	char elmalign;
	Datum data[2], daux, *elems;
	int len;
	HeapTupleHeader hth;
	bool isnull[] = {false, false}, isNullAux;
	ArrayType* arr;
	TupleDesc td;

	elems = (Datum*) palloc(sizeof(*elems) * MAX(n_groups, n_prods));
	td = RelationNameGetTupleDesc("eciout");
	elog(INFO, "%d", td->tdtypeid);

	for (int i = 0; i < n_groups; ++i)
	{
		array_iterate(pm->ar_it, &daux, &isNullAux);
		hth = DatumGetHeapTupleHeader(daux);

		data[0] = GetAttributeByNum(hth, 1, &isNullAux);
		data[1] = Float8GetDatumFast(Kc[i]);
		elems[i] = PointerGetDatum(heap_form_tuple(td, data, isnull));
	}

	get_typlenbyvalalign(td->tdtypeid, &elmlen, &elmbyval, &elmalign);
	arr = construct_array(elems, n_groups, td->tdtypeid, elmlen, elmbyval, elmalign);
	//arr = construct_empty_array(td->tdtypeid);
	elog(INFO, "arr 0");

	td = RelationNameGetTupleDesc("pciout");
	elog(INFO, "%d", td->tdtypeid);

	for (int i = 0; i < n_prods; ++i)
	{
		data[0] = PointerGetDatum(pm->prods[i]);
		data[1] = Float8GetDatumFast(Kp[i]);
		elems[i] = PointerGetDatum(heap_form_tuple(td, data, isnull));
	}

	get_typlenbyvalalign(td->tdtypeid, &elmlen, &elmbyval, &elmalign);
	data[1] = PointerGetDatum(construct_array(elems, n_prods, td->tdtypeid,
		elmlen, elmbyval, elmalign));
	elog(INFO, "arr 1");

	data[0] = PointerGetDatum(arr);
	//data[1] = PointerGetDatum(construct_empty_array(td->tdtypeid));
	elog(INFO, "%d", call_td->tdtypeid);
	pm->vecs_tuple = HeapTupleGetDatum(heap_form_tuple(call_td, data, isnull));
	elog(INFO, "packed");
}

void calc_indexes(FunctionCallInfo fcinfo, perm_mem* pm, index_t index,
	FuncCallContext* funcctx, MemoryContext original_context)
{
	//# Mc e Mp podem virar int*
	double **X, *Xc, *Xp, X_total = 0, *Mc, *Mp, *Kc, *Kp;
	int n_total_prods, n_prods, n_groups = ARR_DIMS(PG_GETARG_ARRAYTYPE_P(0))[0];

	n_total_prods = calc_X(PG_GETARG_ARRAYTYPE_P(0), PG_GETARG_INT32(1), PG_GETARG_INT32(2),
		PG_GETARG_INT32(3) << 1, &X, &Xp, &Xc, pm, index);
	tick("calc_X");
	n_prods = filter_products(Xp, X, pm->prods, n_groups, n_total_prods, index);
	tick("filter_products");

	for (int i = 0; i < n_groups; ++i)
		X_total += Xc[i];

	char (*M)[n_prods] = (decltype(M)) palloc(sizeof(*M) * n_groups);
	Mc = (double*) palloc(sizeof(*Mc) * n_groups);
	Mp = (double*) palloc(sizeof(*Mp) * n_prods);
	//# Conferir valores de Xc e Xp

	calc_M(X, Xp, Xc, n_groups, n_prods, n_total_prods, X_total, (char**) M, Mc, Mp);
	pfree(X);
	pfree(Xc);
	pfree(Xp);
	tick("calc_M");

	/*for (int i = 0; i < n_groups; ++i)
		elog(INFO, "Mc[%d]: %lf", i, Mc[i]);*/

	double (*W)[n_groups] = (decltype(W)) palloc(sizeof(*W) * n_groups);

	calc_W((char**) M, Mc, Mp, n_groups, n_prods, (double**) W);
	pfree(Mc);
	if (index == ECI)
	{
		pfree(M);
		pfree(Mp);
	}
	tick("calc_W");

	Kc = (double*) SPI_palloc(sizeof(*Kc) * n_groups);

	calc_Kc((double**) W, n_groups, Kc);
	pfree(W);
	tick("calc_Kc");

	if (index == ECI)
	{
		z_transform(Kc, n_groups);
		pm->indexes = Kc;
		return;
	}

	Kp = (double*) SPI_palloc(sizeof(*Kp) * n_prods);

	calc_Kp(Kc, n_groups, n_prods, (char**) M, Mp, Kp);
	tick("calc_Kp");
	pfree(M);
	pfree(Mp);

	if (index == PCI)
	{
		SPI_pfree(Kc);

		z_transform(Kp, n_prods);
		pm->indexes = Kp;

		pm->n_prods = n_prods;
		return;
	}

	//ECI_PCI
	z_transform(Kc, n_groups);
	z_transform(Kp, n_prods);
	//MemoryContext old_contex = MemoryContextSwitchTo(original_context);
	pack_indexes(Kc, n_groups, Kp, n_prods, pm, funcctx->tuple_desc);
	//MemoryContextSwitchTo(old_contex);
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
	if (hs_digits > 3 || hs_digits < 1)
		ereport(ERROR,
		(
			errcode(ERRCODE_FEATURE_NOT_SUPPORTED),
			errmsg("hs_digit_pairs must be 1, 2 or 3")
		));

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

	if (index & ECI)
		pm->ar_it = array_create_iterator(groups, 0, NULL);

	/*		Calcula índices      */
	SPI_connect();

	initick();
	calc_indexes(fcinfo, pm, index, funcctx, original_context);
	tick("calc_indexes");

	SPI_finish();
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

	if (index == ECI)
	{
		//Encerra a última chamada
		if (funcctx->call_cntr == *ARR_DIMS(PG_GETARG_ARRAYTYPE_P(0)))
		{
			array_free_iterator(pm->ar_it);
			SRF_RETURN_DONE(funcctx);
		}

		//Cria resultados da chamada
		array_iterate(pm->ar_it, &daux, &isNullAux);
		hth = DatumGetHeapTupleHeader(daux);
		dt[0] = GetAttributeByNum(hth, 1, &isNullAux);
	}
	else
	{
		//Encerra a última chamada
		if (funcctx->call_cntr == pm->n_prods)
			SRF_RETURN_DONE(funcctx);

			//text* debug_t = palloc(VARHDRSZ + 6);
			//SET_VARSIZE(debug_t, VARHDRSZ + 6);
			//memcpy(VARDATA_ANY(debug_t), "IHMOI0", 6);
		dt[0] = PointerGetDatum(pm->prods[funcctx->call_cntr]);
	}
	//elog(INFO, "K: %p %lf", pm->indexes, pm->indexes[0]);

	dt[1] = Float8GetDatumFast(pm->indexes[funcctx->call_cntr]);

	//Cria e retorna tupla
	ht = heap_form_tuple(funcctx->tuple_desc, dt, isnull);
	ret = HeapTupleGetDatum(ht);
	SRF_RETURN_NEXT(funcctx, ret);
}

BG_FUNCTION_INFO_V1(common_eci);

//groups cgroup[], start_year integer, end_year integer, hs_digit_pairs integer
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

	elog(INFO, "Vai sai");
	PG_RETURN_DATUM(pm->vecs_tuple);
}

void query_series(text* c1, text* c2, int start_yi, int end_yi, VarChar* prod)
{
	int len = VARSIZE_ANY_EXHDR(prod);
	mystring* query;

	if (len < 0 || len > 6 || len & 1)
		ereport(ERROR,
		(
			errcode(ERRCODE_FEATURE_NOT_SUPPORTED),
			errmsg("hs_code must be 0, 2, 4 or 6 digits")
		));

#define Q0 "SELECT exporter, product, year, sum(COALESCE((exp_val + imp_val) / 2, exp_val, imp_val))"\
	" FROM transaction WHERE (exporter = '"
#define Q1 "' OR exporter = '"
#define Q2 "')"
#define Q3 " AND year >= "
#define Q4 " AND year <= "
#define Q5 " AND left(product, 0) = '"
#define Q6 " GROUP BY exporter, product, year ORDER BY year, product\0"

	len += sizeof(Q0) + sizeof(Q1) + sizeof(Q2) + sizeof(Q3) + sizeof(Q4) +
		sizeof(Q5) + sizeof(Q6) + 2 * MAXINTSIZE - 5;
	//5 = 1(\') - 6 (sizeof \0s)

	PMYSTRING_INIT(query, VARHDRSZ + len);
	query->litcat(Q0);
	query->concat(VARDATA_ANY(c1), VARSIZE_ANY_EXHDR(c1));
	query->litcat(Q1);
	query->concat(VARDATA_ANY(c2), VARSIZE_ANY_EXHDR(c2));
	query->litcat(Q2);

	if (start_yi != 0)
	{
		query->litcat(Q3);
		query->concat(start_yi);
	}

	if (end_yi != 0)
	{
		query->litcat(Q4);
		query->concat(end_yi);
	}

	if (VARSIZE_ANY_EXHDR(prod) > 0)
	{
		query->litcat(Q5);
		query->data[query->len - 6] += VARSIZE_ANY_EXHDR(prod);
		query->concat(VARDATA_ANY(prod), VARSIZE_ANY_EXHDR(prod));
		query->concat('\'');
	}

	query->concat(Q6, sizeof(Q6));

	//elog(INFO, "query= '%s'", query->data);

	SPI_execute(query->data, true, 0);
	pfree(query);

#undef Q0
#undef Q1
#undef Q2
#undef Q3
#undef Q4
#undef Q5
#undef Q6
}

BG_FUNCTION_INFO_V1(euclidean_distance);

//country_1 text, country_2 text, start_year integer, end_year integer, hs_code varchar(6)
Datum euclidean_distance(PG_FUNCTION_ARGS)
{
	bool is_null;
	l_str p_compare(6);

	SPI_connect();

	query_series(PG_GETARG_TEXT_PP(0), PG_GETARG_TEXT_PP(1),
		PG_GETARG_INT32(2), PG_GETARG_INT32(3), PG_GETARG_VARCHAR_PP(4));

	double dist = 0, aux;

	//if (status <= 0 || SPI_tuptable == NULL)
	//	ERROR

	TupleDesc tupdesc = SPI_tuptable->tupdesc;

	if (!SPI_tuptable->numvals)
	{
		SPI_freetuptable(SPI_tuptable);
		SPI_finish();
		PG_RETURN_FLOAT8(0.0);
	}

	HeapTuple current, last;
	last = SPI_tuptable->vals[0];
	for (int i = 1; i < SPI_tuptable->numvals; i++)
	{
		current = SPI_tuptable->vals[i];

		if (!p_compare(SBI_getString(current, tupdesc, 2, &is_null), SBI_getString(last, tupdesc, 2, &is_null))
			||
			SBI_getInt(current, tupdesc, 3, &is_null) != SBI_getInt(last, tupdesc, 3, &is_null))
		{
			aux = SBI_getDouble(last, tupdesc, 4, &is_null);
			//elog(INFO, "dif: %lf", aux);
			dist += aux * aux;
			last = current;
			continue;
		}

		aux = SBI_getDouble(current, tupdesc, 4, &is_null) - SBI_getDouble(last, tupdesc, 4, &is_null);
		//elog(INFO, "eq: %lf", aux);
		dist += aux * aux;

		if (++i + 1 < SPI_tuptable->numvals)
		{
			last = SPI_tuptable->vals[i];
		}
		else if (i == SPI_tuptable->numvals)
			goto euclidean_distance_last_one_ok;
		else
		{
			last = SPI_tuptable->vals[i];
			break;
		}
	}

	aux = SBI_getDouble(last, tupdesc, 4, &is_null);
	//elog(INFO, "last: %lf", aux);
	dist += aux * aux;

euclidean_distance_last_one_ok:

	SPI_freetuptable(SPI_tuptable);
	SPI_finish();

	PG_RETURN_FLOAT8(sqrt(dist));
}