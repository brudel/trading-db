#include <cstring>


#define MAX(X, Y) ((X) < (Y) ? (X) : (Y))
#define MAXINTSIZE 11
// For hardcoded strings
#define litcat(lit) concat(lit, sizeof(lit) - 1)


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

struct l_str
{
	int n;

	l_str(int _n) : n(_n){}

	bool operator()(char* a, char* b) const
	{
		if (n & 4 && *(int*)a != *(int*)b)
			return false;

		if (n & 2 && ((short*)a)[(n & 4) >> 1] != ((short*)b)[(n & 4) >> 1])
			return false;

		if (n & 1 && a[n & 6] != b[n & 6])
			return false;

		return true;
	}

	size_t operator()(char* str) const
	{
		size_t hash = 0;

		if (n & 4)
			hash |= *(int*)str;

		if (n & 2)
			hash = (hash << 16) | ((short*)str)[(n & 4) >> 1];

		if (n & 1)
			hash = (hash << 8) | str[n & 6];

		return hash * 9973; //31	9973	11400714819323198485
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