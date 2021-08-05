extern "C" {
#include "common/hashfn.h"
}

#include <cstring>


#define MAX(X, Y) ((X) > (Y) ? (X) : (Y))
#define MAXINTSIZE 11
// For hardcoded strings
#define litcat(lit) concat(lit, sizeof(lit) - 1)


 /* Working state for array_iterate() */
 typedef struct ArrayIteratorData
 {
	/* basic info about the array, set up during array_create_iterator() */
	ArrayType  *arr;            /* array we're iterating through */
	bits8      *nullbitmap;     /* its null bitmap, if any */
	int         nitems;         /* total number of elements in array */
	int16       typlen;         /* element type's length */
	bool        typbyval;       /* element type's byval property */
	char        typalign;       /* element type's align property */

	/* information about the requested slice size */
	int         slice_ndim;     /* slice dimension, or 0 if not slicing */
	int         slice_len;      /* number of elements per slice */
	int        *slice_dims;     /* slice dims array */
	int        *slice_lbound;   /* slice lbound array */
	Datum      *slice_values;   /* workspace of length slice_len */
	bool       *slice_nulls;    /* workspace of length slice_len */

	/* current position information, updated on each iteration */
	char       *data_ptr;       /* our current position in the array */
	int         current_item;   /* the item # we're at in the array */
 } ArrayIteratorData;


int itoa(int num, char* str)
{
	char buffer[10];
	int i = 0, j = 0;

	if (num < 0)
	{
		str[j++] = '-';
		num = -num;
	}

	do {
		buffer[i++] = '0' + num % 10;
	} while ( num /= 10 );

	while(i)
		str[j++] = buffer[--i];

	str[j] = '\0';

	return j;
}

struct mystring
{
	int len;
	char data[];

	void concat(mystring* str)
	{
		memcpy(data + len, str->data, str->len);
		len += str->len;
	}

	void concat(char* cstr, int clen)
	{
		memcpy(data + len, cstr, clen);
		len += clen;
	}

	void concat(char c)
	{
		data[len++] = c;
	}

	void concat(int d)
	{
		len += itoa(d, data + len);
	}

	char* end()
	{
		return data + len;
	}
};

struct t_aux
{
	uint fix_len;

	t_aux() : fix_len(0){}
	t_aux(uint n) : fix_len(n){}

	bool operator()(text* a, text* b) const
	{
		int len;

		if (fix_len)
		{
			if (fix_len & 4 && *(int*)VARDATA_ANY(a) != *(int*)VARDATA_ANY(b))
				return false;

			if (fix_len & 2 && ((short*)VARDATA_ANY(a))[(fix_len & 4) >> 1]
					!= ((short*)VARDATA_ANY(b))[(fix_len & 4) >> 1])
				return false;

			if (fix_len & 1 && VARDATA_ANY(a)[fix_len & 6] != VARDATA_ANY(b)[fix_len & 6])
				return false;

			return true;
		}

		if (len = VARSIZE_ANY_EXHDR(a) != VARSIZE_ANY_EXHDR(b))
			return false;

		return !memcmp(VARDATA_ANY(a), VARDATA_ANY(b), len);
	}

	size_t operator()(text* str) const
	{
		size_t hash = 0;

		if (fix_len)
		{
			if (fix_len & 4)
				hash |= *(int*)VARDATA_ANY(str);

			if (fix_len & 2)
				hash = (hash << 16) | ((short*)VARDATA_ANY(str))[(fix_len & 4) >> 1];

			if (fix_len & 1)
				hash = (hash << 8) | VARDATA_ANY(str)[fix_len & 6];

			return hash_bytes_uint32_extended(hash, 0);
		}
		else
			return hash_bytes_extended((unsigned char*) VARDATA_ANY(str), VARSIZE_ANY_EXHDR(str), 0);
	}
};

void z_transform(double* vec, int n)
{
	double mean = 0, sum_quad = 0, stdev;

	for (int i = 0; i < n; ++i)
	{
		mean += vec[i];
		sum_quad += vec[i] * vec[i];
	}

	mean = mean / n;
	stdev = sqrt((sum_quad - mean * mean * n) / n);

	for (int i = 0; i < n; ++i)
		vec[i] = (vec[i] - mean) / stdev;
}