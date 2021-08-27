extern "C" {
#include "common/hashfn.h"
}


/*		Util macros      */
#define MAX(X, Y) ((X) > (Y) ? (X) : (Y)) // Greater element
#define MAXINTSIZE 11 // Max possible size of a int in a decimal string
#define sizeofl(lit) sizeof(lit) -1 // Size of a string literal without '\0'
// Concatenate literal string into a struct mystring
#define litcat(lit) concat(lit, sizeofl(lit))


// Write a decimal into a char* address
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

// char* with length stored and some util functions
struct mystring
{
	int len;
	char data[];

	// Concatenate another mystring in the end of if this
	void concat(mystring* str)
	{
		memcpy(data + len, str->data, str->len);
		len += str->len;
	}

	// Concatenate a char* with in the end of if this, with length given
	void concat(char* cstr, int clen)
	{
		memcpy(data + len, cstr, clen);
		len += clen;
	}

	// Contatenate a char in the end of if this
	void concat(char c)
	{
		data[len++] = c;
	}

	// Contatenate a decial number in the end of if this
	void concat(int d)
	{
		len += itoa(d, data + len);
	}

	// Address after last character
	char* end()
	{
		return data + len;
	}
};

/*	  Struct to hash and compare text
 * 'fix_len' may be an integer not greater than 8 to fixed length texts, otherwise each
 * size is get of text itself.
 */
struct t_aux
{
	uint fix_len;

	t_aux() : fix_len(0){}
	t_aux(uint n) : fix_len(n){}

	// Equal
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

	// hash
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

			return hash * 31; //hash_bytes_uint32_extended 31        9973    11400714819323198485
		}
		else
			return hash_bytes_extended((unsigned char*) VARDATA_ANY(str), VARSIZE_ANY_EXHDR(str), 0) * 31;
	}
};

// Apply a z transform into 'vec' with 'n' elements
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