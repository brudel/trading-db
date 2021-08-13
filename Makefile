ts := sudo docker exec -it ts0
tsparams := -U postgres
psql := ${ts} psql ${tsparams}
pg11 := /usr/include/postgresql/11/server
pg13 := /home/brudel/Área\ de\ trabalho/Faculdade/TCC/Operações/pglib
so_file := /home/brudel/Área\ de\ trabalho/Faculdade/TCC/Docker/tsdocker/trading.so
static := -llapacke -lblas -lgfortran -lm
#static := /usr/lib/x86_64-linux-gnu/liblapacke.a #/usr/lib/x86_64-linux-gnu/libblas.a /usr/lib/x86_64-linux-gnu/liblapack.a
#-Wl,-Bstatic -Wl,-Bdynamic

teste: ${so_file}
	${psql} -c "SELECT * FROM common_eci(ARRAY(SELECT (left(code, 1), array_agg(code))::cgroup FROM country GROUP BY left(code, 1) ORDER BY left(code, 1)), 1998, hs_digit_pairs => 1);"
#"SELECT * FROM countries_eci(ARRAY(SELECT code FROM country ORDER BY code), 1998, hs_digit_pairs => 1);"
#"SELECT * FROM common_eci(ARRAY(SELECT (left(code, 1), array_agg(code))::cgroup FROM country GROUP BY left(code, 1) ORDER BY left(code, 1)), 1998, hs_digit_pairs => 1);"

####	Testes de argumentos    ####
#	Erro de argumento
#ARRAY[]::cgroup[]
#ARRAY[NULL, NULL]::cgroup[]
#ARRAY[('chn', ARRAY['chn'])]::cgroup[]
#ARRAY[(NULL, ARRAY['chn']), ('niu', ARRAY['niu'])]::cgroup[]
#ARRAY[('chn', NULL), ('niu', ARRAY['niu'])]::cgroup[]
#ARRAY[('chn', ARRAY[NULL]), ('niu', ARRAY['niu'])]::cgroup[]
#ARRAY[('chn', ARRAY['1234']), ('niu', ARRAY['niu'])]::cgroup[]

#	Erro esperado de execução
#ARRAY[('chn', ARRAY[]::text[]), ('niu', ARRAY[]::text[])]::cgroup[]
#ARRAY[('chn', ARRAY['chn']), ('niu', ARRAY[]::text[])]::cgroup[]
#ARRAY[('chn', ARRAY[]::text[]), ('niu', ARRAY['niu'])]::cgroup[]

#	Testes válidos
#ARRAY[('chn', ARRAY['chn']), ('niu', ARRAY['niu'])]::cgroup[]
#ARRAY[ARRAY['chn', 'niu'], ARRAY['bra', 'ind']]
#ARRAY[('chn', ARRAY['chn']), ('niu', ARRAY['niu'], ('oie', ARRAY[]::text[]))]::cgroup[] #Remoção de país sem transação no ano
#ARRAY(SELECT code FROM country ORDER BY code)

${so_file}: trading.so
	cp trading.so ${so_file}

trading.so: trading.o
	g++ -shared trading.o -o trading.so ${static}

trading.o: trading.cpp Makefile utils.h
	g++ -fPIC trading.cpp -I${pg13} -c -fno-exceptions -Wno-write-strings

conn:
	${psql}

cteste: teste.cpp
	g++ teste.cpp -o exe
	./exe

vcteste: teste.cpp
	g++ teste.cpp -o exe -g
	valgrind ./exe