
-------------------------------- International trade database create script  --------------------------------

/*
 * Author: Bruno Del Monde <brudel@usp.br>
 * Usage: psql <connection-uri> < import.sql
 */

/* ---- Sumary ----
 *
 * -- Set script configuration --
 *
 * -- Create relations --
 ** Create tables
 ** Create indexes
 ** Clean
 ** Convert to timescale
 ** Fill continents
 *
 * -- Create functions --
 ** Create types
 ** Common functions
 ** Countries functions
 ** Continents functions
 ** Groups functions
 ** Euclidian Distance
 */


---------------- Set script configuration ----------------

\set ON_ERROR_STOP on
SET statement_timeout = 0;
SET lock_timeout = 0;
SET idle_in_transaction_session_timeout = 0;
SET client_encoding = 'UTF8';
SET client_min_messages = warning;
SET row_security = off;
SET work_mem to '256MB';
SET maintenance_work_mem to '4048MB';


---------------- Migrate to final tables ----------------

-- Create tables --

CREATE TABLE IF NOT EXISTS continent
(
    code char(2) PRIMARY KEY,
    name text
);

CREATE TABLE IF NOT EXISTS country
(
    code            char(3) PRIMARY KEY,
    name            text,
    continent       char(2) REFERENCES continent
);

CREATE TABLE IF NOT EXISTS country_annual_measure
(
    country     char(3) REFERENCES country,
    year        integer,
    measure     text,
    value       double precision,
    PRIMARY KEY (country, year, measure)
);

CREATE TABLE IF NOT EXISTS country_group 
(
    name text PRIMARY KEY
);

CREATE TABLE IF NOT EXISTS country_group_belonging
(
    c_group     text REFERENCES country_group,
    country     char(3) REFERENCES country,
    entry_year  integer NOT NULL,
    exit_year   integer,
    PRIMARY KEY (c_group, country, entry_year)
);

CREATE TABLE IF NOT EXISTS product
(
    hs_code         text PRIMARY KEY,
    description     text
);

CREATE TABLE IF NOT EXISTS transaction --# Add constraints latter
(
    year                int,
    exporter            char(3) REFERENCES country(code),
    importer            char(3) REFERENCES country(code),
    product             text REFERENCES product(hs_code),
    exp_val             double precision,
    imp_val             double precision,
    PRIMARY KEY (exporter, importer, product, year)
        WITH (fillfactor = 100)
);

-- Create indexes --

CREATE INDEX
    ON transaction(importer, product)

    WITH
        (fillfactor = 100);

CREATE INDEX
    ON transaction(product, exporter)

    WITH
        (fillfactor = 100);

-- Clean --

TRUNCATE continent CASCADE;
TRUNCATE country CASCADE;
TRUNCATE country_annual_measure CASCADE;
TRUNCATE country_group CASCADE;
TRUNCATE country_group_belonging CASCADE;
TRUNCATE product CASCADE;
TRUNCATE transaction CASCADE;

-- Convert to timescale --

SELECT create_hypertable
(
    'country_annual_measure',
    'year',
    chunk_time_interval => 1,
    create_default_indexes => false,
    if_not_exists => true
);

SELECT create_hypertable
(
    'transaction',
    'year',
    chunk_time_interval => 1,
    create_default_indexes => false,
    if_not_exists => true
);

-- Fill continents

INSERT INTO continent
    (
        code, name
    )
    VALUES
        ('sa', 'South America'),
        ('eu', 'Europe'),
        ('oc', 'Oceania'),
        ('na', 'North America'),
        ('af', 'Africa'),
        ('as', 'Asia');


---------------- Create functions ----------------

-- Create types

CREATE TYPE cgroup
    AS (
        name    text,
        members text[]
    );

CREATE TYPE eciout
    AS (
        "group" text,
        eci     double precision
    );

CREATE TYPE pciout
    AS (
        product text,
        pci     double precision
    );

-- Common functions

CREATE OR REPLACE FUNCTION common_eci
    (
        groups          cgroup[],
        start_year      integer,
        end_year        integer = 0,
        hs_digit_pairs  integer = 3
    )

    RETURNS TABLE
        ("group" text, eci double precision)

    LANGUAGE C
    STRICT
    STABLE
    COST 500000

    AS
        'trading', 'common_eci';

CREATE OR REPLACE FUNCTION common_pci
    (
        groups          cgroup[],
        start_year      integer,
        end_year        integer = 0,
        hs_digit_pairs  integer = 3
    )

    RETURNS TABLE
        ("product" text, pci double precision)

    LANGUAGE C
    STRICT
    STABLE
    COST 500000

    AS
        'trading', 'common_pci';

CREATE OR REPLACE FUNCTION common_eci_pci
    (
        groups          cgroup[],
        start_year      integer,
        end_year        integer = 0,
        hs_digit_pairs  integer = 3
    )

    RETURNS TABLE
        (eci eciout[], pci pciout[])

    LANGUAGE C
    STRICT
    STABLE
    COST 500000

    AS
        'trading', 'common_eci_pci';

-- Countries functions

CREATE OR REPLACE FUNCTION countries_eci
    (
        countries       text[],
        start_year      integer,
        end_year        integer = 0,
        hs_digit_pairs  integer = 3
    )

    RETURNS TABLE
        (country text, eci double precision)

    LANGUAGE SQL
    STRICT
    STABLE
    COST 500000

    AS
        'SELECT common_eci
            (
                ARRAY
                (
                    SELECT
                        (cntr, ARRAY[cntr])::cgroup

                        FROM
                            unnest($1) AS cntr
                ), $2, $3, $4
            )';

CREATE OR REPLACE FUNCTION countries_pci
    (
        countries       text[],
        start_year      integer,
        end_year        integer = 0,
        hs_digit_pairs  integer = 3
    )

    RETURNS TABLE
        (product text, pci double precision)

    LANGUAGE SQL
    STRICT
    STABLE
    COST 500000

    AS
        'SELECT common_pci
            (
                ARRAY
                (
                    SELECT
                        (cntr, ARRAY[cntr])::cgroup

                        FROM
                            unnest($1) AS cntr
                ), $2, $3, $4
            )';

CREATE OR REPLACE FUNCTION countries_eci_pci
    (
        countries       text[],
        start_year      integer,
        end_year        integer = 0,
        hs_digit_pairs  integer = 3
    )

    RETURNS TABLE
        (eci eciout[], pci pciout[])

    LANGUAGE SQL
    STRICT
    STABLE
    COST 500000

    AS
        'SELECT common_eci_pci
            (
                ARRAY
                (
                    SELECT
                        (cntr, ARRAY[cntr])::cgroup

                        FROM
                            unnest($1) AS cntr
                ), $2, $3, $4
            )';

-- Continents functions

CREATE OR REPLACE FUNCTION continents_eci
    (
        continents      text[],
        start_year      integer,
        end_year        integer = 0,
        hs_digit_pairs  integer = 3
    )

    RETURNS TABLE
        (continent text, eci double precision)

    LANGUAGE SQL
    STRICT
    STABLE
    COST 500000

    AS
        'SELECT common_eci
            (
                ARRAY
                (
                    SELECT
                        (cntnnt.code, array_agg(country.code))::cgroup

                        FROM
                            unnest($1) AS cntnnt (code)
                            JOIN country
                                ON
                                    (cntnnt.code = country.continent)

                        GROUP BY
                             cntnnt.code
                ), $2, $3, $4
            )';

CREATE OR REPLACE FUNCTION continents_pci
    (
        continents      text[],
        start_year      integer,
        end_year        integer = 0,
        hs_digit_pairs  integer = 3
    )

    RETURNS TABLE
        (product text, pci double precision)

    LANGUAGE SQL
    STRICT
    STABLE
    COST 500000

    AS
        'SELECT common_pci
            (
                ARRAY
                (
                    SELECT
                        (cntnnt.code, array_agg(country.code))::cgroup

                        FROM
                            unnest($1) AS cntnnt (code)
                            JOIN country
                                ON
                                    (cntnnt.code = country.continent)

                        GROUP BY
                             cntnnt.code
                ), $2, $3, $4
            )';

CREATE OR REPLACE FUNCTION continents_eci_pci
    (
        continents      text[],
        start_year      integer,
        end_year        integer = 0,
        hs_digit_pairs  integer = 3
    )

    RETURNS TABLE
        (eci eciout[], pci pciout[])

    LANGUAGE SQL
    STRICT
    STABLE
    COST 500000

    AS
        'SELECT common_eci_pci
            (
                ARRAY
                (
                    SELECT
                        (cntnnt.code, array_agg(country.code))::cgroup

                        FROM
                            unnest($1) AS cntnnt (code)
                            JOIN country
                                ON
                                    (cntnnt.code = country.continent)

                        GROUP BY
                             cntnnt.code
                ), $2, $3, $4
            )';

-- Groups functions

CREATE OR REPLACE FUNCTION groups_eci
    (
        groups          text[],
        start_year      integer,
        end_year        integer = 0,
        hs_digit_pairs  integer = 3
    )

    RETURNS TABLE
        ("group" text, eci double precision)

    LANGUAGE C
    STRICT
    STABLE
    COST 500000

    AS
        'trading', 'common_eci';

CREATE OR REPLACE FUNCTION groups_pci
    (
        groups          text[],
        start_year      integer,
        end_year        integer = 0,
        hs_digit_pairs  integer = 3
    )

    RETURNS TABLE
        (product text, pci double precision)

    LANGUAGE C
    STRICT
    STABLE
    COST 500000

    AS
        'trading', 'common_pci';

CREATE OR REPLACE FUNCTION groups_eci_pci
    (
        groups          text[],
        start_year      integer,
        end_year        integer = 0,
        hs_digit_pairs  integer = 3
    )

    RETURNS TABLE
        (eci eciout[], pci pciout[])

    LANGUAGE C
    STRICT
    STABLE
    COST 500000

    AS
        'trading', 'common_eci_pci';

CREATE OR REPLACE FUNCTION euclidean_distance
    (
            country_1   text[],
            country_2   text[],
            start_year  integer = 0,
            end_year    integer = 0,
            hs_code     varchar(6) = ''
    )

    RETURNS double precision

    LANGUAGE C 
    STRICT
    STABLE

    AS
        'trading', 'euclidean_distance';

-- Euclidian Distance

CREATE OR REPLACE FUNCTION euclidean_distance
    (
        country_1   text,
        country_2   text[],
        start_year  integer = 0,
        end_year    integer = 0,
        hs_code     varchar(6) = ''
    )

    RETURNS double precision

    LANGUAGE SQL
    STRICT
    STABLE

    AS
        'SELECT euclidean_distance(ARRAY[$1], $2, $3, $4, $5)';

CREATE OR REPLACE FUNCTION euclidean_distance
    (
        country_1   text[],
        country_2   text,
        start_year  integer = 0,
        end_year    integer = 0,
        hs_code     varchar(6) = ''
    )

    RETURNS double precision

    LANGUAGE SQL
    STRICT
    STABLE

    AS
        'SELECT euclidean_distance($1, ARRAY[$2], $3, $4, $5)';

CREATE OR REPLACE FUNCTION euclidean_distance
    (
        country_1   text,
        country_2   text,
        start_year  integer = 0,
        end_year    integer = 0,
        hs_code     varchar(6) = ''
    )

    RETURNS double precision

    LANGUAGE SQL
    STRICT
    STABLE

    AS
        'SELECT euclidean_distance(ARRAY[$1], ARRAY[$2], $3, $4, $5)';