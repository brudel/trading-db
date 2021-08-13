test-db := td0

build: Dockerfile trading.so 002_set_trading_db.sql
	sudo docker build -t brudel/trading-db .

trading.so: ../trading.so
	cp ../trading.so .

../%:
	make -C ../ $*

002_set_trading_db.sql: ../set_trading_db.sql
	cp ../set_trading_db.sql 002_set_trading_db.sql

test: build
	mkdir -p pgdata
	sudo chown 70:$$USER pgdata
	sudo docker run -d --name ${test-db} -v "$$PWD""/pgdata:/var/lib/postgresql/data/" brudel/trading-db
	sudo docker start ${test-db}
	sleep 3
	sudo chmod g+rwx pgdata

test-conn:
	sudo docker exec -it ${test-db} psql -U postgres

clean-img:
	sudo docker image rm brudel/trading-db

clean-test:
	sudo docker stop ${test-db}
	sudo docker rm ${test-db}

clean-docker:
	make clean-test
	make clean-img

clean-files:
	rm trading.so 002_set_trading_db.sql

push:
	sudo docker push brudel/trading-db