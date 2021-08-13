build: Dockerfile trading.so 002_set_trading_db.sql
	sudo docker build -t brudel/trading-db .

trading.so: ../trading.so
	cp ../trading.so .

../%:
	make -C ../ $*

002_set_trading_db.sql: ../set_trading_db.sql
	cp ../set_trading_db.sql 002_set_trading_db.sql

test: build
	sudo docker run -d --name ts1 brudel/trading-db
	sudo docker start ts1

test-conn:
	sudo docker exec -it ts1 psql -U postgres

clean-img:
	sudo docker image rm brudel/trading-db

clean-test:
	sudo docker stop ts1
	sudo docker rm ts1

clean-docker:
	make clean-test
	make clean-img

clean-files:
	rm trading.so 002_set_trading_db.sql

push:
	sudo docker push brudel/trading-db