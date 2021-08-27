test-db := td0

build: .build

.build: Dockerfile trading.so 002_set_trading_db.sql
	sudo docker build -t brudel/trading-db .
	touch ../.so_flag # Updated in cointainer with the image
	touch .build

trading.so: ../trading.so
	cp ../trading.so .

../%: FORCE
	make -C ../ $*

FORCE:

002_set_trading_db.sql: ../set_trading_db.sql
	cp ../set_trading_db.sql 002_set_trading_db.sql

test: build
	mkdir pgdata
	sudo docker run -d --name ${test-db} -v "$$PWD""/pgdata:/var/lib/postgresql/data/" brudel/trading-db
	make start-test
	sudo chown -R 70:$$USER pgdata

start-test:
	sudo docker start ${test-db}
	sleep 3
	sudo chmod -R g+rwx pgdata

test-conn:
	sudo docker exec -it ${test-db} psql -U postgres

clean-docker: test-conn
	sudo docker image rm brudel/trading-db

clean-test:
	sudo docker stop ${test-db}
	sudo docker rm ${test-db}

clean-files:
	rm trading.so 002_set_trading_db.sql .build

clean-all: clean-docker clean-files

push:
	sudo docker push brudel/trading-db