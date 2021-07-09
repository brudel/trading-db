ts := sudo docker exec -it ts0
tsparams := -U postgres
Alocalparams := timescale_teste
pg11 := /usr/include/postgresql/11/server
pg13 := "/home/brudel/Área de trabalho/Faculdade/TCC/Operações/pglib"
so_file := "/home/brudel/Área de trabalho/Faculdade/TCC/Docker/tsdocker/trading.so"
static := -llapacke -lblas -lgfortran -lm
#static := /usr/lib/x86_64-linux-gnu/liblapacke.a #/usr/lib/x86_64-linux-gnu/libblas.a /usr/lib/x86_64-linux-gnu/liblapack.a
#-Wl,-Bstatic -Wl,-Bdynamic

teste: ${so_file}
	${ts} psql ${tsparams} ${localparams} -c "SELECT * FROM common_eci(array(SELECT (left(code, 1), array_agg(code))::cgroup FROM country GROUP BY left(code, 1) ORDER BY left(code, 1)), 2015, 3);"
#"SELECT * FROM common_eci(array(SELECT (code, ARRAY[code])::cgroup FROM country ORDER BY code), 2015, 3);"
#"SELECT * FROM common_eci(array(SELECT (left(code, 1), array_agg(code))::cgroup FROM country GROUP BY left(code, 1) ORDER BY left(code, 1)), 2015, 3);"
#"SELECT * FROM calculate_eci(ARRAY[('oie','{a,b}'), ('euro','{a,b,m,kk}'),  ('g2','{a,b,m}')]::cgroup[], 2015, 3);"

${so_file}: trading.o
	g++ -fPIC -shared trading.o -I${pg13} -o ${so_file} ${static}

trading.o: trading.cpp Makefile
	g++ -fPIC -shared trading.cpp -I${pg13} -c -fpermissive -fno-exceptions

redef:
	${ts} psql ${tsparams} ${localparams} -c \
	"DROP FUNCTION IF EXISTS common_eci;"
	${ts} psql ${tsparams} ${localparams} -c \
	"CREATE OR REPLACE FUNCTION common_eci(groups cgroup[], year integer, hs_digits integer)\
\
    RETURNS TABLE\
        (\"group\" text, eci double precision)\
\
    AS\
        '/tsdata/trading.so', 'common_eci'\
\
    LANGUAGE\
        C STRICT;"

conn:
	${ts} psql ${tsparams}

cteste: teste.c
	gcc teste.c -o exe
	./exe

vcteste: teste.c
	gcc teste.c -o exe -g
	valgrind ./exe