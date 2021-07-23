
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

struct l_str
{
	int n;

	l_str(int _n) : n(_n){}

	bool operator()(char* a, char* b) const
	{
		if (n & 4 && *(int*)a != *(int*)b)
			{if(n==6)elog(INFO, "out0: '%.6s' vs '%.6s'", a, b);return false;}

		if (n & 2 && ((short*)a)[(n & 4) >> 1] != ((short*)b)[(n & 4) >> 1])
			{if(n==6)elog(INFO, "out1: '%.6s' vs '%.6s'", a, b);return false;}

		if (n & 1 && a[n & 6] != b[n & 6])
			{if(n==6)elog(INFO, "out2: '%.6s' vs '%.6s'", a, b);return false;}

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

		return hash * 31;//11400714819323198485;
	}
};