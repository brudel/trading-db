extern "C" {
#include "postgres.h"
#include "fmgr.h" //Module magic
#include "funcapi.h" //Return row
#include "utils/array.h" //Coisas de array
#include "executor/spi.h" //SPI
#include <lapacke.h> //AlgeLin

PG_MODULE_MAGIC;
}

#include <unordered_map>

//# Usar umas função hash descente 



#define BG_FUNCTION_INFO_V1(func) extern "C" { \
PG_FUNCTION_INFO_V1(common_eci); \
}

struct country_equal
{
	bool operator()(char* a, char* b) const
	{
		return *(short*)a == *(short*)b && a[2] == b[2];
	}
};

struct country_hash
{
	size_t operator()(char* str) const
	{
		return (size_t)((*(short*)str) << 16) | (size_t) str[2];
	}
};

struct product_equal
{
	bool operator()(char* a, char* b) const
	{
		return *(int*)a == *(int*)b && ((short*)a)[2] == ((short*)b)[2];
	}
};

struct product_hash
{
	size_t operator()(char* str) const
	{
		return (size_t)((*(int*)str) << 16) | (size_t) ((short*)str)[2];
	}
};

typedef struct perm_mem {
	ArrayIterator ar_it;
	double* ecis;
	int n;
} perm_mem;

#define talloc(size) palloc(VARHDRSZ + size)

#define SBI_getString(tuple, tupdesc, col, is_null) (char*) \
	VARDATA_ANY( DatumGetTextP( \
		SPI_getbinval(tuple, tupdesc, col, is_null) \
	))

typedef std::unordered_map<char*, int, country_hash, country_equal> country_map;
typedef std::unordered_map<char*, int, product_hash, product_equal> product_map;

double get_double_from_SPI(HeapTuple row, TupleDesc rowdesc, int colnumber)
{
	double ret;
	char* res;

	res = SPI_getvalue(row, rowdesc, colnumber);
	ret = atof(res);

	pfree(res);
	return ret;
}

void create_all_countrs_map(country_map *countrs_map)
{
	//# Usar variável global SPI_tuptable e desalocar depois
	int count = 0;
	bool is_null;

	int status = SPI_execute(
		"SELECT code "
		"FROM country ",
		true, 0);


	if (status > 0 && SPI_tuptable != NULL)
	{
		SPITupleTable *tuptable = SPI_tuptable;
		TupleDesc tupdesc = tuptable->tupdesc;

		for (int i = 0; i < tuptable->numvals; i++)
		{
			HeapTuple tuple = tuptable->vals[i];

			countrs_map->insert({SBI_getString(tuple, tupdesc, 1, &is_null), ++count});
		}
	}
}

void create_common_countrs_map(country_map *countrs_map, ArrayType* groups)
{
	int count = 0;
	Datum daux;
	HeapTupleHeader hth;
	ArrayType* members;
	bool is_null;
	text* taux;

	ArrayIterator itg = array_create_iterator(groups, 0, NULL);
	while (array_iterate(itg, &daux, &is_null))
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

void create_prods_map(product_map *prods_map, int hs_digits)
{
	int count = 0;
	bool is_null;

	int status = SPI_execute(
		"SELECT hs_code "
		"FROM product ",
		true, 0);

	if (status > 0 && SPI_tuptable != NULL)
	{
		TupleDesc tupdesc = SPI_tuptable->tupdesc;

		//elog(INFO, "create prods: %d", SPI_tuptable->numvals);
		for (int i = 0; i < SPI_tuptable->numvals; i++)
		{
			HeapTuple tuple = SPI_tuptable->vals[i];

			prods_map->insert({{SBI_getString(tuple, tupdesc, 1, &is_null), count++}});
		}
	}
	SPI_freetuptable(SPI_tuptable);
}

void calc_X(country_map* countrs_map, int n_groups, product_map* prods_map, int year, int hs_digits,	
	double** _X, double* Xp, double* Xc)
{
	double (*X)[prods_map->size()] = (void*)_X;
	bool is_null;

	//Primeira iteração do for em seguida, mas iniciando Xp
	for (int i = 0; i < prods_map->size(); ++i)
	{
		X[0][i] = 0;
		Xp[i] = 0;
	}
	Xc[0] = 0;

	//Inicia X e Xc
	for (int i = 1; i < n_groups; ++i)
	{
		//elog(INFO, "Bucket[%d] > %d", i, countrs_map->bucket_size(i));
		for (int j = 0; j < prods_map->size(); ++j)
			X[i][j] = 0;

		Xc[i] = 0;
	}

	//elog(INFO, "Total %d", countrs_map->bucket_count());

	int status = SPI_execute(
		"SELECT exporter, left(product, 2 * 3), sum(exp_val) "
		"FROM transaction "
		"WHERE year = 2015 "
		"GROUP BY exporter, left(product, 2 * 3)",
		true, 0);

	if (status <= 0 && SPI_tuptable == NULL)
		return;

	SPITupleTable *tuptable = SPI_tuptable;
	TupleDesc tupdesc = tuptable->tupdesc;

	for (int i = 0; i < tuptable->numvals; i++)
	{
		HeapTuple tuple = tuptable->vals[i];

		int c = (*countrs_map)[SBI_getString(tuple, tupdesc, 1, &is_null)],
			p = (*prods_map)[SBI_getString(tuple, tupdesc, 2, &is_null)];

		double value = atof(SPI_getvalue(tuple, tupdesc, 3));

		X[c][p] += value;
		//# Mais operações do que calcular em x[c][p] pronto, mas não precisa percorrer dnv
		Xc[c] += value;
		Xp[p] += value;

		//K[(*countrs_map)[SPI_getvalue(tuple, tupdesc, 1)]] += atof(SPI_getvalue(tuple, tupdesc, 3));
	}
	SPI_freetuptable(SPI_tuptable);
}

//https://www.sanfoundry.com/c-program-implement-interpolation-search-array-integers/
int interpolation_search(int* vec, int top, int value)
{
	//elog(INFO, "\tsearch in %d %d", top, value);
	int bottom = 0, mid;
	--top;

	//O condicional precisa ser disjunto em bottom == top
	//	para ser falso e evitar divisão por zero.
	while (vec[bottom] < value && value <= vec[top])
	{
		mid = bottom + (top - bottom) * (value - vec[bottom]) / (vec[top] - vec[bottom]);
		//elog(INFO, "\tsearch loop %d %d %d %d %d", bottom, mid, top, vec[bottom], vec[top]);

		if (value == vec[mid]){
			//elog(INFO, "\tsearch find");
			return mid;}

		if (value < vec[mid])
			top = mid - 1;

		else
			bottom = mid + 1;
	}
	//elog(INFO, "\tsearch out not");

	return value == vec[bottom] ? bottom : -1;
}

//# Se usar valor mínimo > 0 vai precisar ajustar Xc pra tirar os produtos excluídos
void filter_products(double* Xp, double** _X, int n_groups, product_map* prods_map)
{
	//elog(INFO, "Entro");
	double (*X)[prods_map->size()] = (void*)_X;
	int count = 0, idx, m_count, aux;
	int* eliminated = palloc(sizeof(*eliminated)* prods_map->size()), *moved;

	//Identifica
	for (int i = 0; i < prods_map->size(); ++i)
		if (Xp[i] == 0)
			eliminated[count++] = i;

	if (count == 0)
		return;
	//elog(INFO, "Identifico");


	//Cria moved[]
	moved = palloc(sizeof(*eliminated) * count);

	m_count = 0;
	aux = prods_map->size() - count;

	idx = count;
	while (eliminated[--idx] >= aux);

	for (; ++idx < count; ++aux)
		while (aux != eliminated[idx])
			moved[m_count++] = aux++;

	while (aux < prods_map->size())
		moved[m_count++] = aux++;

	aux = prods_map->size() - count;

	//Substituí
	for (auto it = prods_map->begin(); it != prods_map->end(); ++it)
	{
		//elog(INFO, "for in");
		while (idx = interpolation_search(eliminated, count, it->second) >= 0)
		{
			//elog(INFO, "while in");
			it = prods_map->erase(it);
			if (it == prods_map->end())
				goto end_loop;
			//elog(INFO, "while out");
		}

		//elog(INFO, "for mid");

		if (it->second >= aux)
		{
			//elog(INFO, "if in");
			idx = interpolation_search(moved, m_count, it->second);

			if (idx == -1)
				elog(INFO, "DEU MERDAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA AQUI Ó: %d", it->second);

			//elog(INFO, "\tinter");
			//elog(INFO, "\tatrib");

			Xp[eliminated[idx]] = Xp[it->second];

			for (int i = 0; i < n_groups; ++i)
				X[i][eliminated[idx]] = X[i][it->second];
			
			it->second = eliminated[idx];
		}
		//elog(INFO, "for out");
	}
end_loop:

	for (int i = 0; i < prods_map->size(); ++i)
		if (Xp[i] == 0)
			elog(INFO, "Deu merda");

	elog(INFO, "NULLs: %d", count);
	pfree(eliminated);
	pfree(moved);
}

void calc_M(double** _X, double* Xp, double* Xc, int n_groups, product_map* prods_map, int total,
	char**_M, double* Mc, double* Mp)
{
	char (*X)[prods_map->size()] = (void*) _X;
	char (*M)[prods_map->size()] = (void*) _M;

	for (int i = 0; i < prods_map->size(); ++i)
		Mp[i] = 0; //# Da pra iniciar com o primeiro valor do for do M[0][i]

	//elog(INFO, "for i %d j %d", n_groups, prods_map->size());

	for (int i = 0; i < n_groups; ++i)
	{
		Mc[i] = 0;
		for (int j = 0; j < prods_map->size(); ++j)
		{
			//elog(INFO, "M[%d][%d] = %lf * %lf / (%lf * %lf)", i, j, X[i][j], total, Xc[i],  Xp[j]);
			M[i][j] = (X[i][j] * total / (Xc[i] * Xp[j]) >= 1 ? 1 : 0);
			Mc[i] += M[i][j];
			Mp[j] += M[i][j];
		}
		//elog(INFO, "i: %d", i);
	}
}

void calc_Kc(double** W, int n_groups, double* K)
{
	double* avlr = palloc(sizeof(*avlr)*n_groups);
	double* avli = palloc(sizeof(*avli)*n_groups);
	double (*avtr)[n_groups] = palloc(sizeof(*avtr)*n_groups);

	int info = LAPACKE_dgeev(LAPACK_ROW_MAJOR, 'N', 'V', n_groups, (double*)W, n_groups, avlr, avli, NULL, n_groups, (double*)avtr, n_groups);
	elog(INFO, "info %d", info);

	int max = 1;

	//# Se o 1 sempre for o primeiro
	for (int i = 2; i < n_groups; ++i)
		if (avlr[i] > max)
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

void calc_eci(country_map* countrs_map, int n_groups, product_map* prods_map, int year, int hs_digits, perm_mem* pm)
{
	//# Mc e Mp podem virar int* 
	double *Xc, *Xp, *Mc, *Mp, *K, total = 0;

	double (*X)[prods_map->size()] = palloc(sizeof(*X)*n_groups);
	Xc = palloc(sizeof(*Xc)*n_groups);
	Xp = palloc(sizeof(*Xp)*prods_map->size());

	calc_X(countrs_map, n_groups, prods_map, year, hs_digits, (double**) X, Xp, Xc);
	filter_products(Xp, (double**) X, n_groups, prods_map);

	for (int i = 0; i < n_groups; ++i)
		total += Xc[i];

	char (*M)[prods_map->size()] = palloc(sizeof(*M)*n_groups);
	Mc = palloc(sizeof(*Mc)*n_groups);
	Mp = palloc(sizeof(*Mp)*prods_map->size());
	//# Conferir valores de Xc e Xp

	calc_M((double**) X, Xp, Xc, n_groups, prods_map, total, (char**) M, Mc, Mp);
	pfree(X);
	pfree(Xc);
	pfree(Xp);

	double (*W)[n_groups] = palloc(sizeof(*W)*n_groups);

	//Dava pra pular produtos que um país não tem especialidade
	for (int i = 0; i < n_groups; ++i)
	{
		for (int j = 0; j < n_groups; ++j)
		{
			W[i][j] = 0;
			for (int p = 0; p < prods_map->size(); ++p)
				//# Da pra fazer condicional (bool** M;) ou conversão
				W[i][j] += (M[i][p] & M[j][p]) / Mp[p];
			W[i][j] /= Mc[i];
		}
		//elog(INFO, "i: %d", i);
	}
	pfree(M);
	pfree(Mc);
	pfree(Mp);

	K = SPI_palloc(sizeof(*K)*n_groups);

	calc_Kc((double**) W, n_groups, K);

	pm->ecis = K;
}

BG_FUNCTION_INFO_V1(common_eci);

Datum common_eci(PG_FUNCTION_ARGS)
{
	FuncCallContext *funcctx;
	ArrayType* groups;
	ArrayMetaState *mstate = NULL;
	TupleDesc td0;
	HeapTuple ht;
	HeapTupleHeader hth;
	Datum dt[2], ret, daux;
	bool isnull[2] = {false, false}, isNullAux;
	perm_mem* pm;

	groups = PG_GETARG_ARRAYTYPE_P(0);
	//hth = DatumGetHeapTupleHeader();

	//Primeira chamada: cálculo dos valores
	if (SRF_IS_FIRSTCALL())
	{
		//Cria contexto de chamada
		funcctx = SRF_FIRSTCALL_INIT();

		//Identifica tipo de retorno
		if (get_call_result_type(fcinfo, NULL, &td0) != TYPEFUNC_COMPOSITE)
			ereport(ERROR,
			(
				errcode(ERRCODE_FEATURE_NOT_SUPPORTED),
				errmsg("function returning record called in context that cannot accept type record")
			));
		funcctx->tuple_desc = BlessTupleDesc(td0);

		//Cria memória do contexto
		pm = (perm_mem*) palloc(sizeof(perm_mem));
		funcctx->user_fctx = pm;

		pm->ar_it = array_create_iterator(groups, 0, mstate);

		//Calcula ECI
		SPI_connect();

		country_map* countrs_map = new country_map(241);
		product_map* prods_map = new product_map(5300);
		create_common_countrs_map(countrs_map, groups);
		create_prods_map(prods_map, PG_GETARG_INT32(2));
		elog(INFO, "main, created maps: %d %d", countrs_map->size(), prods_map->size());

		calc_eci(countrs_map, *ARR_DIMS(groups), prods_map, PG_GETARG_INT32(1), PG_GETARG_INT32(2), pm);

		delete prods_map, countrs_map;
		SPI_finish();
		elog(INFO, "final first");
	}

	//Recupera contexto da chamada
	funcctx = SRF_PERCALL_SETUP();
	pm = (perm_mem*) funcctx->user_fctx;

	//Encerra a última chamada
	if (funcctx->call_cntr == *ARR_DIMS(groups))
	{
		array_free_iterator(pm->ar_it);
		SRF_RETURN_DONE(funcctx);
	}

	//Cria resultados da chamada
	array_iterate(pm->ar_it, &daux, &isNullAux);
	hth = DatumGetHeapTupleHeader(daux);
	dt[0] = GetAttributeByNum(hth, 1, &isNullAux);
	dt[1] = Float8GetDatum(pm->ecis[funcctx->call_cntr]);
   
	//Cria e retorna tupla
	ht = heap_form_tuple(funcctx->tuple_desc, dt, isnull);
	ret = HeapTupleGetDatum(ht);
	SRF_RETURN_NEXT(funcctx, ret);
}