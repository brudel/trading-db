extern "C" {
#include "postgres.h" // Main Postgres stuff
#include "fmgr.h" // Postgres function interface
#include "funcapi.h" // Return row types
#include "utils/array.h" // Array stuff
#include "utils/lsyscache.h" // get_typlenbyvalalign
#include "executor/spi.h" // SPI: Queries interface
#include <lapacke.h> // Eigenvectors and eigenvalues
// Commented because do not improve execution time
//#include <omp.h> // Parallelization

// Postgres libraries version check
PG_MODULE_MAGIC;
}

#include <unordered_map>
#include "utils.h"

/*		Notes      */
/*  Continent prefixed functions lack a significant optimization as it query for
 * transactions could be grouped by continent rather than by country. But this would
 * implicate a laborious interface change.
 */

/*		Configuration macros      */
// Log time of principal subroutines
//#define PROFILE
#define COUNTRY_SIZE 3
#define PRODUCT_SIZE 6
// Value to be considered for the calculations
#define REAL_VALUE "COALESCE((exp_val + imp_val) / 2, exp_val, imp_val, 0)"

// Subroutines time profile implementation
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


// Declare function for Postgres with C ABI
#define BG_FUNCTION_INFO_V1(func) extern "C" { \
PG_FUNCTION_INFO_V1(func); \
}

// Get text from HeapTuple
#define SBI_getText(tuple, tupdesc, col, is_null) DatumGetTextPP( \
		SPI_getbinval(tuple, tupdesc, col, is_null) \
	)

// Get data from text from HeapTuple
#define SBI_getString(tuple, tupdesc, col, is_null) ((char*) \
	VARDATA_ANY(SBI_getText(tuple, tupdesc, col, is_null)))

// Get double from HeapTuple
#define SBI_getDouble(tuple, tupdesc, col, is_null) DatumGetFloat8( \
		SPI_getbinval(tuple, tupdesc, col, is_null) \
	)

// Get int from HeapTuple
#define SBI_getInt(tuple, tupdesc, col, is_null) DatumGetInt32( \
		SPI_getbinval(tuple, tupdesc, col, is_null) \
	)

// Initiate mstr mystring with n max length
#define PMYSTRING_INIT(mstr, n) mstr = (mystring*) palloc(n + sizeof(mystring)); \
mstr->len = 0;

// Get ArrayType number of elements
#define ARR_DIM(v) ArrayGetNItems(ARR_NDIM(v), ARR_DIMS(v))

// Test if HeapTupleHeader has null elements
#define HeapTupleHeaderHasNulls(hth) \
	(((hth)->t_infomask & HEAP_HASNULL) != 0)


// Store output parameters
typedef struct perm_mem {
	int n_indexes; // eci and pci
	text** prods; // pci and eci_pci
	text** cntrs; // eci and eci_pci
	double* indexes; // eci and pci
	Datum vecs_tuple; // eci_pci
} perm_mem;

// Map text to int with equal and hash function on t_aux
typedef std::unordered_map<text*, int, t_aux, t_aux> t_map;

// Index output type (function suffix)
enum index_t {
	ECI = 1,
	PCI,
	ECI_PCI
};


/*	  Populate map from country of a group to index of a group.
 * 'groups' are cgroup[] and contain the groups to populate 'countrs_map'. Also validate
 * groups vector.
 * Optionally construct a group name vector on 'pm'->cntrs. (if 'index' have ECI flag).
 */
void create_common_countrs_map(t_map* countrs_map, ArrayType* groups, perm_mem* pm, index_t index)
{
	int count = 0;
	Datum daux;
	HeapTupleHeader hth;
	ArrayType* members;
	bool is_null;
	text* taux;

	if (index & ECI)
		pm->cntrs = (text**) palloc(sizeof(*pm->cntrs) * ARR_DIM(groups));

	// Iterate through cgroups
	ArrayIterator itv = array_create_iterator(groups, 0, NULL);
	while (array_iterate(itv, &daux, &is_null))
	{
		hth = DatumGetHeapTupleHeader(daux);

		if (HeapTupleHeaderHasNulls(hth))
		{
			delete countrs_map;
			elog(ERROR, "groups has null elements");
		}

		members = DatumGetArrayTypeP(GetAttributeByNum(hth, 2, &is_null));

		if (ARR_HASNULL(members))
		{
			delete countrs_map;
			elog(ERROR, "cgroup members has null elements");
		}

		ArrayIterator itv = array_create_iterator(members, 0, NULL);

		if (index & ECI)
			pm->cntrs[count] = DatumGetTextPP(GetAttributeByNum(hth, 1, &is_null));

		// Iterate through cgroups members
		while (array_iterate(itv, &daux, &is_null))
		{
			taux = DatumGetTextPP(daux);
			if (VARSIZE_ANY_EXHDR(taux) != COUNTRY_SIZE)
			{
				delete countrs_map;
				elog(ERROR, "country code text must be of length %d", COUNTRY_SIZE);
			}

			// Add group member to map with index of the group
			if (!countrs_map->insert({taux, count}).second)
			{
				delete countrs_map;
				elog(ERROR, "country belongs to more than one group");
			}
		}

		++count;
	}
}

/*	  Populate map from country_group names to respective indexes.
 * 'groups' are text[] and contain the country_group names to populate 'countrs_map'.
 * Also validate groups vector.
 * Optionally construct a country_group name vector on 'pm'->cntrs. (if 'index' have ECI flag).
 */
void create_groups_countrs_map(t_map* countrs_map, ArrayType* groups, perm_mem* pm, index_t index)
{
	int count = 0;
	Datum daux;
	bool is_null;
	text* taux;

	if (index & ECI)
		pm->cntrs = (text**) palloc(sizeof(*pm->cntrs) * ARR_DIM(groups));

	// Add country_group name to map
	ArrayIterator itv = array_create_iterator(groups, 0, NULL);
	while (array_iterate(itv, &daux, &is_null))
	{
		taux = DatumGetTextPP(daux);

		if (VARSIZE_ANY_EXHDR(taux) != COUNTRY_SIZE)
		{
			delete countrs_map;
			elog(ERROR, "country code text must be of length %d", COUNTRY_SIZE);
		}

		if (index & ECI)
			pm->cntrs[count] = taux;

		if (!countrs_map->insert({taux, count++}).second)
		{
			delete countrs_map;
			elog(ERROR, "two groups with same name");
		}
	}
}

/*	  Populate map from product aggregations to respective indexes.
 * 'hs_digits' contains the aggregation level to be queried from product table to populate
 * 'countrs_map'.
 * Optionally construct a product code vector on 'pm'->prods. (if 'index' have PCI flag).
 */
int create_prods_map(t_map* prods_map, int hs_digits, perm_mem* pm, index_t index,
	t_map* countrs_map)
{
	int n;
	TupleDesc td;
	HeapTuple tuple;
	text* prod;
	bool is_null;

	/*		Query product aggregations      */
#define Q "SELECT DISTINCT left(hs_code, 0) FROM product"
	char* query = (char*) palloc(sizeof(Q));
	memcpy(query, Q, sizeof(Q));
	query[30] |= hs_digits; // OR char zero, 00110000, equivalent to sum
#undef Q

	int status = SPI_execute(query, true, 0);
	pfree(query);

	if (status <= 0 || SPI_tuptable == NULL)
	{
		delete countrs_map, prods_map;
		elog(ERROR, "can't successfully access needed data on database");
	}

	/*		Create Map      */
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

/*	  Query transactions for the given parameters.
 * If 'f_year' equals zero, the period it's only s_year. Else it's the open interval
 * between 's_year' and 'f_year'. 'c_groups' indicate if the function group from
 * country_group and 'hs_digits' the aggregation level of the products.
 */
void calc_X_query(int s_year, int f_year, int hs_digits, bool c_groups,
	t_map* countrs_map, t_map* prods_map)
{
	int status;
	mystring* query;

// Queried entity
#define QEXP "exporter"
#define QC_G "c_group"

/*		Parts of the query string      */
#define QSELECT "SELECT "
#define QSCOLUMNS ", left(product, 0), sum(" REAL_VALUE ") FROM transaction "
/*	Optimization possibility
 * Add c_group restriction on join clause to get only the ones passed as arguments,
 * would be needed to construct a c_group's list string similarly as done for countries
 * in euclidean distance functions query.
 */
#define QJOIN "JOIN country_group_belonging ON " \
	"(exporter = country AND year >= entry_year AND (exit_year IS NULL OR year < exit_year)) "
#define QWHERE "WHERE year "
#define QWEQ "= "
#define QWGE ">= "
#define QWFYEAR " AND year <= "
// GROUP BY makes it a lot faster than group only on matrix filling
#define QGROUP " GROUP BY "
#define QGPROD ", left(product, 0)\0"

	// Allocate query string
	PMYSTRING_INIT(query, sizeofl(QSELECT) + 2 * MAX(sizeofl(QEXP),
		sizeofl(QC_G)) + sizeofl(QSCOLUMNS) + sizeofl(QJOIN) + sizeofl(QWHERE)
		+ MAX(sizeofl(QWEQ) + MAXINTSIZE, sizeofl(QWGE) + sizeofl(QWFYEAR) +
		2 * MAXINTSIZE) + sizeofl(QGROUP) + sizeofl(QGPROD));

	query->litcat(QSELECT);

	// Queried entity
	if (c_groups)
		query->litcat(QC_G);
	else
		query->litcat(QEXP);

	query->litcat(QSCOLUMNS);
	query->data[query->len - 81] |= hs_digits;
	// OR char zero, 00110000, equivalent to sum

	if (c_groups)
		query->litcat(QJOIN);

	query->litcat(QWHERE);

	// Interval
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

	// Queried entity
	if (c_groups)
		query->litcat(QC_G);
	else
		query->litcat(QEXP);

	query->litcat(QGPROD);
	query->data[query->len - 3] |= hs_digits;
	// OR char zero, 00110000, equivalent to sum

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
	{
		delete countrs_map, prods_map;
		elog(ERROR, "can't successfully access needed data on database");
	}
}

/*	  Allocate and calculate X matrix and sum vectors.
 * If 'f_year' equals zero, the period it's only s_year. Else it's the open interval
 * between 's_year' and 'f_year'. 'c_groups' indicate if the function group from
 * country_group and 'hs_digits' the aggregation level of the products.
 */
int calc_X(ArrayType* groups, int s_year, int f_year, int hs_digits,
	double*** _X, double** _Xp, double** _Xc, perm_mem* pm, index_t index)
{
	double *Xc, *Xp;
	int c, p, n_prods, n_groups = ARR_DIM(groups);
	bool is_null, c_groups = groups->elemtype == TEXTOID;
	t_map* countrs_map;
	t_map* prods_map = new t_map(0, t_aux(hs_digits), t_aux(hs_digits));
	TupleDesc td;
	HeapTuple tuple;
	double value;

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
	n_prods = create_prods_map(prods_map, hs_digits, pm, index, countrs_map);
	tick("prods_map");

	/*		Allocate and initialize arrays      */
	// Points the argument references to the correct memory intervals
	double (*X)[n_prods] = (decltype(X)) SPI_palloc(sizeof(*X) * n_groups);
	*_X = (double**) X;
	*_Xc = Xc = (double*) SPI_palloc(sizeof(*Xc) * n_groups);
	*_Xp = Xp = (double*) SPI_palloc(sizeof(*Xp) * n_prods);


	// Splitted first iteration of the next double for to initialize Xp
	for (int i = 0; i < n_prods; ++i)
	{
		X[0][i] = 0;
		Xp[i] = 0;
	}
	Xc[0] = 0;

	// Initialize X and Xc
	for (int i = 1; i < n_groups; ++i)
	{
		for (int j = 0; j < n_prods; ++j)
			X[i][j] = 0;

		Xc[i] = 0;
	}
	tick("calc_X set");

	/*		Populate arrays      */
	calc_X_query(s_year, f_year, hs_digits, c_groups, countrs_map, prods_map);
	tick("SPI_execute");

	td = SPI_tuptable->tupdesc;

	if (!SPI_tuptable->numvals)
	{
		delete countrs_map, prods_map;
		elog(ERROR, "no transactions were found with these parameters");
	}

	for (int i = 0; i < SPI_tuptable->numvals; i++)
	{
		tuple = SPI_tuptable->vals[i];

		auto it = countrs_map->find(SBI_getText(tuple, td, 1, &is_null));
		if (it == countrs_map->end()) // Not of interest for the function
			continue;

		c = it->second;
		p = (*prods_map)[SBI_getText(tuple, td, 2, &is_null)];

		value = SBI_getDouble(tuple, td, 3, &is_null);

		X[c][p] += value;
		/*  For multimembered groups sum here imply more operations than iterate ahead
		 * through X, but avoid a whole through matrix iteration.
		 */
		Xc[c] += value;
		Xp[p] += value;
	}
	SPI_freetuptable(SPI_tuptable);

	SPI_finish();
	delete countrs_map, prods_map;

	return n_prods;
}

/*	  Remove groups with no non-zero transactions.
 * Iterate through 'Xc' for groups with no non-zero transactions, assuming no negative
 * transaction value. Correct X('_X') and 'Xc' (and 'cntrs' if 'index' flags ECI)
 * swapping the invalid groups with the last ones and decreasing n_groups.
 */
// If minimum value > 0, would need to correct Xp by subtracting excluded products
int filter_groups(double* Xc, double** _X, text** cntrs, int n_groups, int n_prods, index_t index)
{
	double (*X)[n_prods] = (decltype(X)) _X;

	// Verify if there are some valid group
	while (--n_groups >= 0 && Xc[n_groups] == 0);

	if (n_groups < 0)
		elog(ERROR, "no transactions were found with non-zero value for these parameters");
	n_groups++;

	// Iterate through Xc
	for (int i = 0; i < n_groups; ++i)
		if (Xc[i] == 0)
		{
			while (Xc[--n_groups] == 0); // Find last valid group to swap
			if (n_groups < i) // Ends loop if the are no more valid groups
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

/*	  Remove products with no non-zero transactions.
 * Iterate through 'Xp' for products with no non-zero transactions, assuming no negative
 * transaction value. Correct X and 'Xp' (and 'prods' if 'index' flags PCI) swapping the
 * invalid products with the last ones and decreasing n_prods.
 */
// If minimum value > 0, would need to correct Xc by subtracting excluded countries
int filter_products(double* Xp, double** _X, text** prods, int n_groups, int n_prods, index_t index)
{
	double (*X)[n_prods] = (decltype(X)) _X;

	// Iterate through Xp
	for (int i = 0; i < n_prods; ++i)
		if (Xp[i] == 0)
		{
			while (Xp[--n_prods] == 0); // Find last valid product to swap
			if (n_prods < i) // Ends loop if the are no more valid products
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

/*	  Calculate M matrix and sum vectors from X.
 * 'n_total_prods' is the leading dimension of X, with may be different from n_prods.
 */
void calc_M(double** _X, double* Xp, double* Xc, int n_groups, int n_prods, int n_total_prods,
	double X_total, bool**_M, double* Mc, double* Mp)
{
	double (*X)[n_total_prods] = (decltype(X)) _X;
	bool (*M)[n_prods] = (decltype(M)) _M;

	// Splitted first iteration of the next double for to initialize Mp
	Mc[0] = 0;
	for (int j = 0; j < n_prods; ++j)
	{
		M[0][j] = (X[0][j] * X_total / (Xc[0] * Xp[j]) >= 1 ? 1 : 0);
		Mc[0] += M[0][j];
		Mp[j] = M[0][j];
	}

	// Calculation
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

//	Calculate W matrix from M.
void calc_W(bool**_M, double* Mc, double* Mp, int n_groups, int n_prods, double** _W)
{
	bool (*M)[n_prods] = (decltype(M)) _M;
	double (*W)[n_groups] = (decltype(W)) _W;

	// #pragma omp parallel for
	/*	Optimization possibility
	 * Sum only products that each country is specialized in, loops would have to be
	 * swapped to country : product : country, rather than country : country : product,
	 * but that would spoil spatial locality unless W are transposed.
	 */
	for (int i = 0; i < n_groups; ++i)
		for (int j = 0; j < n_groups; ++j)
		{
			W[i][j] = 0;
			for (int p = 0; p < n_prods; ++p)
				/* This also can be than by a conditional sum of 1/Mp[p] to avoid float
				 * point division in exchange of a branch operation.
				 */
				W[i][j] += (M[i][p] & M[j][p]) / Mp[p];
			W[i][j] /= Mc[i];
		}
}

//	Calculate Kc vector from W.
void calc_Kc(double** W, int n_groups, double* Kc)
{
	double *avlr, *avli, mean = 0, sum_quad = 0, stdev;
	int info, max[2];

	// Matrix for the (right) eigenvectors
	double (*avtr)[n_groups] = (decltype(avtr)) palloc(sizeof(*avtr) * n_groups);
	avlr = (double*) palloc(sizeof(*avlr) * n_groups); // Real part of eigenvalues
	avli = (double*) palloc(sizeof(*avli) * n_groups); // Imaginary part of eigenvalues, should be zeros

	// Calculate eigenvalues and eigenvectors
	info = LAPACKE_dgeev(LAPACK_ROW_MAJOR, // matrix_layout
		'N', // Don't calculate left eigenvectors
		'V', //	Calculate (right) eigenvectors
		n_groups, // Order of the matrix
		(double*) W,
		n_groups, // Leading dimension of W
		avlr, // (OUT)
		avli, // (OUT)
		NULL, // (OUT) Matrix for the left eigenvectors
		n_groups, // Leading dimension previous matrix
		(double*) avtr, // (OUT)
		n_groups); // Leading dimension previous matrix

	if (info)
		elog(ERROR, "LAPACKE_dgeev returned error code %d", info);

	/*		Find the second greater eigenvalue     */
	// max[0] must be the greater one.
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

	// Move right eigenvector to Kc
	for (int i = 0; i < n_groups; ++i)
		Kc[i] = avtr[i][max[1]];

	pfree(avlr);
	pfree(avli);
	pfree(avtr);
}

//	Calculate Kp vector from Kc.
void calc_Kp(double* Kc, int n_groups, int n_prods, bool** _M, double* Mp, double* Kp)
{
	bool (*M)[n_prods] = (decltype(M)) _M;

	// Splitted first iteration of the next for to initialize Kp
	for (int j = 0; j < n_prods; ++j)
		Kp[j] = M[0][j] * Kc[0];

	// Sum specialized exporters
	for (int i = 1; i < n_groups; ++i)
		for (int j = 0; j < n_prods; ++j)
			Kp[j] += M[i][j] * Kc[i];

	// Make it a mean
	for (int j = 0; j < n_prods; ++j)
		Kp[j] /= Mp[j];
}

/*	Pack output vectors into a HeapTuple so it can be returned.
 * This function is called for functions with eci_pci suffix, so 'Kc' and 'pm'->cntrs
 * must be packed with Kp and 'pm'->prods so both can be returned.
 * 'call_td': TupleDesc of the return tuple.
 */
void pack_indexes(double* Kc, int n_groups, double* Kp, int n_prods, perm_mem* pm,
	TupleDesc call_td)
{
	short elmlen; // Type length
	bool elmbyval; // False if type is passed by reference
	char elmalign; // Type memory alignment
	Datum data[2], *elems;
	bool isnull[] = {false, false};
	ArrayType* arr; // To temporally store ECI array
	TupleDesc td;

	elems = (Datum*) palloc(sizeof(*elems) * MAX(n_groups, n_prods));
	/*		Pack ECI      */
	td = RelationNameGetTupleDesc("eciout");

	// Pack Kc with 'pm'->cntrs into a Postgres vector of tuples
	for (int i = 0; i < n_groups; ++i)
	{
		data[0] = PointerGetDatum(pm->cntrs[i]);
		data[1] = Float8GetDatumFast(Kc[i]);
		elems[i] = PointerGetDatum(heap_form_tuple(td, data, isnull)->t_data);
	}

	get_typlenbyvalalign(td->tdtypeid, &elmlen, &elmbyval, &elmalign);
	arr = construct_array(elems, n_groups, td->tdtypeid, elmlen, elmbyval, elmalign);

	/*		Pack PCI      */
	td = RelationNameGetTupleDesc("pciout");

	// Pack Kp with 'pm'->prods into a Postgres vector of tuples
	for (int i = 0; i < n_prods; ++i)
	{
		data[0] = PointerGetDatum(pm->prods[i]);
		data[1] = Float8GetDatumFast(Kp[i]);
		elems[i] = PointerGetDatum(heap_form_tuple(td, data, isnull)->t_data);
	}

	get_typlenbyvalalign(td->tdtypeid, &elmlen, &elmbyval, &elmalign);
	data[1] = PointerGetDatum(construct_array(elems, n_prods, td->tdtypeid,
		elmlen, elmbyval, elmalign));

	data[0] = PointerGetDatum(arr); // ECI array

	// Pack both vectors into a tuple
	pm->vecs_tuple = HeapTupleGetDatum(heap_form_tuple(call_td, data, isnull));
}

/*	Calculate appropriate index (ECI, PCI or both).
 * 'fcinfo': All function arguments.
 * 'pm': Indexes store struct.
 * 'index': index to calculate.
 * 'call_td': TupleDesc of the return tuple (for ECI_PCI).
 */
void calc_indexes(FunctionCallInfo fcinfo, perm_mem* pm, index_t index, TupleDesc call_td)
{
	// Mc and Mp can have int* type if converted to double on use time.
	double **X, *Xc, *Xp, X_total = 0, *Mc, *Mp, **W, *Kc, *Kp;
	bool** M;
	int n_total_prods, n_prods, n_groups = ARR_DIM(PG_GETARG_ARRAYTYPE_P(0));

	/*		Calculate X matrix and filter valid entities     */
	n_total_prods = calc_X(PG_GETARG_ARRAYTYPE_P(0), PG_GETARG_INT32(1), PG_GETARG_INT32(2),
		PG_GETARG_INT32(3) << 1, &X, &Xp, &Xc, pm, index);
	tick("calc_X");
	n_groups = filter_groups(Xc, X, pm->cntrs, n_groups, n_total_prods, index);
	tick("filter_groups");
	n_prods = filter_products(Xp, X, pm->prods, n_groups, n_total_prods, index);
	tick("filter_products");

	for (int i = 0; i < n_groups; ++i)
		X_total += Xc[i];

	/*		Calculate M matrix     */
	M = (bool**) palloc(sizeof(*M) * n_groups * n_prods);
	Mc = (double*) palloc(sizeof(*Mc) * n_groups);
	Mp = (double*) palloc(sizeof(*Mp) * n_prods);

	calc_M(X, Xp, Xc, n_groups, n_prods, n_total_prods, X_total, M, Mc, Mp);
	pfree(X);
	pfree(Xc);
	pfree(Xp);
	tick("calc_M");


	/*		Calculate W matrix     */
	W = (double**) palloc(sizeof(*W) * n_groups * n_groups);

	calc_W(M, Mc, Mp, n_groups, n_prods, W);
	pfree(Mc);
	if (index == ECI)
	{
		pfree(M);
		pfree(Mp);
	}
	tick("calc_W");

	/*		Calculate Kc vector     */
	Kc = (double*) palloc(sizeof(*Kc) * n_groups);

	calc_Kc(W, n_groups, Kc);
	pfree(W);
	tick("calc_Kc");

	// Return ECI
	if (index == ECI)
	{
		z_transform(Kc, n_groups);
		pm->indexes = Kc;

		pm->n_indexes = n_groups;
		return;
	}

	/*		Calculate Kp vector     */
	Kp = (double*) palloc(sizeof(*Kp) * n_prods);

	calc_Kp(Kc, n_groups, n_prods, M, Mp, Kp);
	tick("calc_Kp");
	pfree(M);
	pfree(Mp);

	// Return PCI
	if (index == PCI)
	{
		pfree(Kc);

		z_transform(Kp, n_prods);
		pm->indexes = Kp;

		pm->n_indexes = n_prods;
		return;
	}

	/*		Pack ECI and PCI to be returned     */
	z_transform(Kc, n_groups);
	z_transform(Kp, n_prods);
	pack_indexes(Kc, n_groups, Kp, n_prods, pm, call_td);
}

// Run in on first call of the functions, validate args and prepare output.
void common_index_init(FunctionCallInfo fcinfo, index_t index)
{
	FuncCallContext *funcctx;
	perm_mem* pm;
	TupleDesc td;
	ArrayType* groups = PG_GETARG_ARRAYTYPE_P(0);
	int hs_digits = PG_GETARG_INT32(3);
	MemoryContext original_context;

	/*		Validate args      */
	if (ARR_HASNULL(groups))
		elog(ERROR, "groups has null elements");

	if (ARR_DIM(groups) < 2)
		elog(ERROR, "groups must have at least 2 elements");

	if (hs_digits > 3 || hs_digits < 1)
		elog(ERROR, "hs_digit_pairs must be 1, 2 or 3");

	/*		Create multi-call environment      */
	funcctx = SRF_FIRSTCALL_INIT();
	original_context = MemoryContextSwitchTo(funcctx->multi_call_memory_ctx);

	// Identify return type
	if (get_call_result_type(fcinfo, NULL, &td) != TYPEFUNC_COMPOSITE)
		ereport(ERROR,
		(
			errcode(ERRCODE_FEATURE_NOT_SUPPORTED),
			errmsg("function returning record called in context that cannot accept type record")
		));
	funcctx->tuple_desc = BlessTupleDesc(td);

	// Create memory to return data
	pm = (perm_mem*) palloc(sizeof(perm_mem));
	funcctx->user_fctx = pm;

	/*		Calculate indexes      */

	initick();
	calc_indexes(fcinfo, pm, index, funcctx->tuple_desc);
	tick("calc_indexes");

	MemoryContextSwitchTo(original_context);
}

/*	Construct and return rows of a set.
 * Extract data to return from 'fcinfo's pm, according with index type passed as 'index',
 * build HeapTuple and return a row of a time.
 */
Datum return_table(FunctionCallInfo fcinfo, index_t index)
{
	FuncCallContext *funcctx;
	HeapTuple ht;
	HeapTupleHeader hth;
	Datum dt[2], ret, daux;
	bool isnull[2] = {false, false}, isNullAux;
	perm_mem* pm;

	// Retrieve multi-call environment
	funcctx = SRF_PERCALL_SETUP();
	pm = (perm_mem*) funcctx->user_fctx;

	// Ends last call
	if (funcctx->call_cntr == pm->n_indexes)
		SRF_RETURN_DONE(funcctx);

	// Get data to return
	if (index == ECI)
		dt[0] = PointerGetDatum(pm->cntrs[funcctx->call_cntr]);
	else
		dt[0] = PointerGetDatum(pm->prods[funcctx->call_cntr]);

	dt[1] = Float8GetDatumFast(pm->indexes[funcctx->call_cntr]);

	// Build and return row
	ht = heap_form_tuple(funcctx->tuple_desc, dt, isnull);
	ret = HeapTupleGetDatum(ht);
	SRF_RETURN_NEXT(funcctx, ret);
}

/*	  Arguments for the next three functions
 * groups (cgroup[] or text[] for groups prefixed),
 * start_year integer,
 * end_year integer,
 * hs_digit_pairs integer
 */

BG_FUNCTION_INFO_V1(common_eci);

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

/*	  Add next country from Postgres array iterator 'itv' on 'query' mytring and set
 * flags 'mask' on 'countrs_map' map.
 */
void add_country(mystring* query, ArrayIterator itv, t_map* countrs_map, int mask)
{
	Datum daux;
	bool is_null;
	text* country;

	array_iterate(itv, &daux, &is_null);
	country = DatumGetTextPP(daux);

	if (VARSIZE_ANY_EXHDR(country) != COUNTRY_SIZE)
	{
		delete countrs_map;
		elog(ERROR, "country code text must be of length %d", COUNTRY_SIZE);
	}

	query->concat(VARDATA_ANY(country), COUNTRY_SIZE);

	(*countrs_map)[country] |= mask;
	// It can be done because C++ initialize integers with 0
}

// Construct a list string of the countries of 'g1' and 'g2' on 'query' and a map of it.
t_map* add_countries(mystring* query, ArrayType* g1, ArrayType* g2)
{
	int n;
	Datum daux;
	bool is_null;
	text* country;
	ArrayIterator itv;
	t_map* countrs_map;

	// Validate input vectors
	if (ARR_HASNULL(g1) || ARR_HASNULL(g2) || !(ARR_NDIM(g1) && ARR_NDIM(g2)))
		elog(ERROR, "country array is empty or has null elements");
	
	countrs_map = new t_map((ARR_DIM(g1) + ARR_DIM(g2)) * 1.3, t_aux(COUNTRY_SIZE),
		t_aux(COUNTRY_SIZE));

#define QSEP "', '"

	// Iterate through g1
	itv = array_create_iterator(g1, 0, NULL);
	n = ARR_DIM(g1);
	for (int i = 0; i < n; ++i)
	{
		add_country(query, itv, countrs_map, 1);
		query->litcat(QSEP);
	}

	// Iterate through g2
	itv = array_create_iterator(g2, 0, NULL);
	n = ARR_DIM(g2) - 1;
	for (int i = 0; i < n; ++i)
	{
		add_country(query, itv, countrs_map, 2);
		query->litcat(QSEP);
	}

	// Adds g2 last element
	add_country(query, itv, countrs_map, 2);

	return countrs_map;
}

/*	  Query transactions for the given series.
 * Search for transactions of countries that belong to 'g1' or 'g2', whose product belong
 * to 'prod' hs group, this is, products that start with prod digits. Optional limit by
 * 'start_yi' and 'end_yi' as each of then are set, including then.
 */
t_map* query_series(ArrayType* g1, ArrayType* g2, int start_yi, int end_yi, VarChar* prod)
{
	int len = VARSIZE_ANY_EXHDR(prod);
	mystring* query;
	t_map* countrs_map;
	int status;

	// Validate prod argument
	if (len < 0 || len > 6 || len & 1)
		elog(ERROR, "hs_code must have 0, 2, 4 or 6 digits");

/*		Parts of the query string      */
#define QSELECT "SELECT exporter, product, year, sum(" REAL_VALUE ")"\
	" FROM transaction WHERE exporter IN ('"
#define QFCLOSE "')"
#define QWSYEAR " AND year >= "
#define QWFYEAR " AND year <= "
#define QWPROD " AND left(product, 0) = '"
#define QGROUP " GROUP BY exporter, product, year ORDER BY year, product\0"

	// Calculate max query string length
	len = sizeofl(QSELECT) + sizeofl(QFCLOSE) + sizeofl(QWSYEAR) + sizeofl(QWFYEAR) +
		sizeofl(QWPROD) + sizeofl(QGROUP) + 2 * MAXINTSIZE + COUNTRY_SIZE +
		(sizeofl(QSEP) + COUNTRY_SIZE) * (ARR_NDIM(g1) + ARR_NDIM(g2) - 1);

	PMYSTRING_INIT(query, len);
	query->litcat(QSELECT);

	// Add countries to query of
	countrs_map = add_countries(query, g1, g2);

	query->litcat(QFCLOSE);

	// Bottom limit
	if (start_yi != 0)
	{
		query->litcat(QWSYEAR);
		query->concat(start_yi);
	}

	// Upper limit
	if (end_yi != 0)
	{
		query->litcat(QWFYEAR);
		query->concat(end_yi);
	}

	// Filter by product hs group
	if (VARSIZE_ANY_EXHDR(prod) > 0)
	{
		query->litcat(QWPROD);
		query->data[query->len - 6] |= VARSIZE_ANY_EXHDR(prod);
		// OR char zero, 00110000, equivalent to sum
		query->concat(VARDATA_ANY(prod), VARSIZE_ANY_EXHDR(prod));
		query->concat('\'');
	}

	query->litcat(QGROUP);

#undef QSEP
#undef QSELECT
#undef QFCLOSE
#undef QWSYEAR
#undef QWFYEAR
#undef QWPROD
#undef QGROUP

	status = SPI_execute(query->data, true, 0);
	pfree(query);

	if (status <= 0 || SPI_tuptable == NULL)
	{
		delete countrs_map;
		elog(ERROR, "can't successfully access needed data on database");
	}

	return countrs_map;
}

BG_FUNCTION_INFO_V1(euclidean_distance);

/*	  euclidean_distance arguments
 * country_1 text(optionally [])
 * country_2 text(optionally [])
 * start_year integer
 * end_year integer
 * hs_code varchar(6)
 */
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

	// Query series
	countrs_map = query_series(PG_GETARG_ARRAYTYPE_P(0), PG_GETARG_ARRAYTYPE_P(1),
		PG_GETARG_INT32(2), PG_GETARG_INT32(3), PG_GETARG_VARCHAR_PP(4));

	td = SPI_tuptable->tupdesc;

	if (!SPI_tuptable->numvals)
	{
		SPI_freetuptable(SPI_tuptable);
		delete countrs_map;
		SPI_finish();
		PG_RETURN_FLOAT8(0.0);
	}

	// Calculate distance²
	for (int i = 0; i < SPI_tuptable->numvals;)
	{
		current = first = SPI_tuptable->vals[i]; // First value of a dimension

		g1 = g2 = 0;

		// Accumulate dimension values on g1 and g2
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
				SBI_getText(first, td, 2, &is_null)) // Same product
			&& SBI_getInt(current, td, 3, &is_null) == SBI_getInt(first, td, 3, &is_null));
				//Same year

		aux = g1 - g2;

		dist += aux * aux;
	}

	SPI_freetuptable(SPI_tuptable);
	SPI_finish();
	delete countrs_map;

	PG_RETURN_FLOAT8(sqrt(dist));
}