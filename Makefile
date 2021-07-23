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
	${ts} psql ${tsparams} ${localparams} -c "SELECT * FROM common_eci(array(SELECT (left(code, 1), array_agg(code))::cgroup FROM country GROUP BY left(code, 1) ORDER BY left(code, 1)), 2015, 2015, 3);"
#"SELECT euclidean_distance('chn', 'uni')"
#"SELECT * FROM common_eci(array(SELECT (exporter, ARRAY[exporter])::cgroup FROM transaction WHERE year = 2015 GROUP BY exporter ORDER BY exporter), 2015);"
#"SELECT * FROM common_eci(array(SELECT (left(code, 1), array_agg(code))::cgroup FROM country GROUP BY left(code, 1) ORDER BY left(code, 1)), 2015);"

${so_file}: trading.o
	g++ -fPIC -shared trading.o -I${pg13} -o ${so_file} ${static}

trading.o: trading.cpp Makefile text.h
	g++ -fPIC -shared trading.cpp -I${pg13} -c -fpermissive -fno-exceptions

conn:
	${ts} psql ${tsparams}

cteste: teste.cpp
	g++ teste.cpp -o exe
	./exe

vcteste: teste.cpp
	g++ teste.cpp -o exe -g
	valgrind ./exe