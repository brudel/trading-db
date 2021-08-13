# syntax = docker/dockerfile:1
FROM timescale/timescaledb:latest-pg13
ENV POSTGRES_PASSWORD secret
COPY 002_set_trading_db.sql /docker-entrypoint-initdb.d/
COPY trading.so /usr/local/lib/postgresql/
RUN apk add lapack