test-db := td0
pglib := pglib
ts := sudo docker exec -it ${test-db}
tsparams := -U postgres
psql := ${ts} psql ${tsparams}
so_flag := .so_flag
static := -llapacke -lblas -lgfortran -lm

test: test-update
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

start-test:
	make -C Docker start-test

install-reqs:
	sudo apt-get install git g++ make docker.io liblapacke-dev gfortran

test-update: ${so_flag}

${so_flag}: trading.so
	touch ${so_flag}
	sudo docker cp trading.so ${test-db}:/usr/local/lib/postgresql/

Docker/pgdata:
	make -C Docker test

trading.so: trading.o
	g++ -shared trading.o -o trading.so ${static}

trading.o: trading.cpp Makefile utils.h ${pglib}
	g++ -fPIC trading.cpp -I${pglib} -c -fno-exceptions -Wno-write-strings

${pglib}:
	sudo docker run -d --name td-aux postgres:13-alpine
	sudo docker start td-aux
	sudo docker cp td-aux:/usr/local/include/postgresql/server/ ${pglib}
	sudo docker stop td-aux
	sudo docker rm td-aux

conn:
	${psql}

cteste: teste.cpp
	g++ teste.cpp -o exe -I${pglib} ${static}
	./exe

vcteste: teste.cpp
	g++ teste.cpp -o exe -g -I${pglib} ${static}
	valgrind ./exe

clean:
	sudo rm -rf ${pglib} trading.so trading.o ${so_flag}