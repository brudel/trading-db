# Trading-DB

This repository contains material related to the thesis [*Armazenamento e Processamento de Dados de Comercio Internacional Usando TimescaleDB*](https://github.com/brudel/trading-db/blob/main/Monografia.pdf), including the SQL script to create the data base, the code for the implemented operations and the configuration of docker image that encapsulate the base.

Note: Arrays passed as input to functions, can have any dimensionality.

## Make

`trading-db$ make install-reqs`: Install compilation requisites.

`trading-db/Docker$ make build`: Build the shared library and brudel/trading-db image.

`trading-db/Docker$ make test`: Create and start test container.

`trading-db$ make start-test`: Start test container.

`trading-db$ make conn`: Open psql client at test container.

`trading-db$ make`: Execute query specified at Makefile at test container.


## Schema

### continent
Column | Type | Constraints 
------ | ---- | ----------
code | char(2) | PRIMARY KEY
name | text |

<br/>

### country
Column | Type | Constraints
------ | ---- | ----------
code | char(3) | PRIMARY KEY
name | text | 
continent | char(2) | REFERENCES continent

<br/>

### country_annual_measure
Column | Type | Constraints
------ | ---- | ----------
country | char(3) | PRIMARY KEY, REFERENCES country
year | integer | PRIMARY KEY
measure | text | PRIMARY KEY
value | double precision |

<br/>

### country_group 
Column | Type | Constraints
------ | ---- | ----------
name | text | PRIMARY KEY

<br/>

### country_group_belonging
Column | Type | Constraints
------ | ---- | ----------
c_group | text | PRIMARY KEY, REFERENCES country_group
country | char(3) | PRIMARY KEY, REFERENCES country
entry_year | integer | PRIMARY KEY
exit_year | integer |

<br/>

### product
Column | Type | Constraints
------ | ---- | ----------
hs_code | text | PRIMARY KEY
description | text |

<br/>

### transaction
Column | Type | Constraints
------ | ---- | ----------
product | text | PRIMARY KEY, REFERENCES product
importer | char(3) | PRIMARY KEY, REFERENCES country
exporter | char(3) | PRIMARY KEY, REFERENCES country
year | int | PRIMARY KEY
exp_val | double precision |
imp_val | double precision |

<br/>


## Complexity indexes API

Set of function to economic complexity index (ECI) and product complexity index (PCI) calculation.

### Functions:

* `countries_eci`
* `countries_pci`
* `countries_eci_pci`

<br/>

* `continents_eci`
* `continents_pci`
* `continents_eci_pci`

<br/>

* `groups_eci`
* `groups_pci`
* `groups_eci_pci`

<br/>

* `common_eci`
* `common_pci`
* `common_eci_pci`

### Arguments

#### Required

Column | Type | Constraints
---- | ---- | ---------
groups | text[] | Array of countries or aggregations.
start_year | integer | Year to calculate the index, or first year of the index calculation interval.

#### Optional

Name | Type | Default | Description
---- | ---- | ------------ | ---------
end_year | integer | 0 | If defined, last year of the index calculation interval.
hs_digit_pairs | integer | 3 | Number of HS digit pairs to be considered. Must be 1, 2 or 3.

### Prefixes

`countries`: The argument `groups` contains a set of codes of the country to be considered at indexes calculation.

`continents`: The argument `groups` contains a set of codes of the continents to be considered at indexes calculation, with the transactional values of the countries belonging to each continent aggregated.

`groups`: The argument `groups` contains a set of group of `country_group` names to be considered at index calculation. The calculation takes in account the entry and exit of countries to the group during the time interval and a same country may be in more than one group.

`common`: Functions with this prefix exceptionally have `group` argument with `cgroup[]` type. Here, group is a  representation of an arbitrary named group, with the respective countries list. One country can't belong to more than one group and each group members are static for the calculation interval. `cgroup` is defined as:

`cgroup = (name text, members text[])`

### Return and suffixes

`eci`: Produces a tuple set with format `(name text, eci double precision)`, `name` is an country or continent code or an group name.

`pci`: Produces a tuple set with format `(product text, pci double precision)`.

`eci_pci`: Produces a tuple with format `((name text, eci double precision)[], (product text, pci double precision))`, with contains two vectors, with are the ECI and PCI respective result vectors. The functions take virtually the same time no matter the suffix, therefore functions with `eci_pci` suffix are indicated to optimize situations where both measures are wanted.

### Sample Usage

All countries ECI in 1998:
~~~ SQL
SELECT * FROM countries_eci(array(SELECT code FROM country), 1998);
~~~

Continent aggregated PCI, temporally aggregated between 2000 and 2006:
~~~ SQL
SELECT * FROM continents_pci(array(SELECT code FROM continent), 2000, 2006);
~~~

ECI and PCI aggregated by country code first letter and two digit pairs HS product code in 2015:
~~~ SQL
SELECT * FROM
	common_eci_pci(
	array(
		SELECT
			(left(code, 1), array_agg(code))::cgroup
			FROM
				country
			GROUP BY
				left(code, 1)
	),
	2015,
	hs_digit_pairs => 2
	);
~~~


## euclidean_distance function

Calculate the distance between two transaction temporal series with different product and interval aggregations. The function returns a `double precision` value.

### Arguments

#### Required

Column | Type | Constraints
---- | ---- | ---------
country_1 | text or text[] | Country 1.
country_2 | text or text[] | Country 2.

#### Optional

Name | Type | Default | Description
---- | ---- | ------------ | ---------
start_year | integer | 0 | Define a inferior limit for the interval (limit included).
end_year | integer | 0 | Define a superior limit for the interval (limit included).
hs_code | varchar(6) | '' | Product or by HS aggregation.

### Sample usage

Euclidean distance between two country series:
~~~ SQL
SELECT euclidean_distance('chn', 'eua');
~~~

Euclidean distance between a country and two other countries, from the beginning of time up to 1980:
~~~ SQL
SELECT euclidean_distance('chn', ARRAY['bra', 'chl'], end_year => 1980);
~~~

Euclidean distance between two country series, from 1980 up to 1999, for products inside coffee, tea, mate e spices category:
~~~ SQL
SELECT euclidean_distance('chn', 'eua', 1980, 1999, '09');
~~~
