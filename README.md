# Armazenamento e Processamento de Dados de Comercio Internacional Usando TimescaleDB

Esse repositório contém material relacionado ao trabalho Armazenamento e Processamento de Dados de Comercio Internacional Usando TimescaleDB, inclundo o script SQL para a criação da bade de dados, o código das operações implementadas e a configuração da imagem docker que encapsula a base.

Observação: As funções que recebem vetores, podem receber vetores de qualquer dimencionalidade.

## Esquema

### continent
 Coluna | Tipo | Restrições
 ------ | ---- | ----------
 code | char(2) | PRIMARY KEY
 name | text |

<br/>

### country
 Coluna | Tipo | Restrições
 ------ | ---- | ----------
 code | char(3) | PRIMARY KEY
 name | text | 
 continent | char(2) | REFERENCES continent

<br/>

### country_annual_measure
 Coluna | Tipo | Restrições
 ------ | ---- | ----------
 country | char(3) | PRIMARY KEY, REFERENCES country
 year | integer | PRIMARY KEY
 measure | text | PRIMARY KEY
 value | double precision |

<br/>

### country_group 
 Coluna | Tipo | Restrições
 ------ | ---- | ----------
 name | text | PRIMARY KEY

<br/>

### country_group_belonging
 Coluna | Tipo | Restrições
 ------ | ---- | ----------
 c_group | text | PRIMARY KEY, REFERENCES country_group
 country | char(3) | PRIMARY KEY, REFERENCES country
 entry_year | integer | PRIMARY KEY
 exit_year | integer |

<br/>

### product
 Coluna | Tipo | Restrições
 ------ | ---- | ----------
 hs_code | text | PRIMARY KEY
 description | text |

<br/>

### transaction
 Coluna | Tipo | Restrições
 ------ | ---- | ----------
 product | text | PRIMARY KEY, REFERENCES product
 importer | char(3) | PRIMARY KEY, REFERENCES country
 exporter | char(3) | PRIMARY KEY, REFERENCES country
 year | int | PRIMARY KEY
 exp_val | double precision |
 imp_val | double precision |

<br/>


## API índices de complexidade

Conjunto de funções para cálculo de índices de complexidade econômica (ECI) e de produto (PCI).

### Funções:

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

### Argumentos

#### Obrigatórios

Nome | Tipo | Descrição
---- | ---- | ---------
groups | text[] | Vetor com os países ou agregações.
start_year | integer | Ano a ser calculado o índice, ou primeiro ano do intervalo a ser calculado.

#### Opcionais

Nome | Tipo | Valor padrão | Descrição
---- | ---- | ------------ | ---------
end_year | integer | 0 | Se definido, último ano do intervalo a ser calculado o índice.
hs_digit_pairs | integer | 3 | Define quantos pares de dígitos do HS serão considerados. Deve ser 1, 2 ou 3.

### Prefixos

`countries`: O argumento `groups` contém uma lista de códigos de países a serem considerados para o cálculo dos índices.

`continents`: O argumento `groups` contém uma lista de códigos de continentes a serem considerados para o cálculo dos índices, sendo os valores transacionais dos países pertencentes a cada continente agregados.

`groups`: O argumento `groups` contém uma lista de nomes de grupos de países pertencentes a `country_group` a serem considerados para o cálculo dos índices. O cálculo considera a entrada e saída de países de um grupo durante o intervalo de tempo e um mesmo país pode estar em múltiplos grupos.

`common`: Funções com esse prefixo excepcionalmente tem o argumento `groups` com o tipo `cgroup[]`. Sendo esse tipo a representação arbitrária de um grupo nomeado com a lista dos respectivos países membros. Um país não pode estar em mais de um grupo e os membros de cada grupo são estáticos para o intervalo de cálculo. `cgroup` é definido como:

`cgroup = (name text, members text[])`

### Retorno e sufixos

`eci`: Produz um conjunto de tuplas no formato `(name text, eci double precision)`, sendo `name` o código do país ou continente ou o nome do grupo.

`pci`: Produz um conjunto de tuplas no formato `(product text, pci double precision)`.

`eci_pci`: Produz uma tupla no formato `((name text, eci double precision)[], (product text, pci double precision))`, contendo dois vetores, sendo eles os respectivos resultados do cálculo de ECI e PCI. As funções demoram virtualmente o mesmo tempo independentemente do sufixo, portanto funções com esse sufixo são indicadas para otimizar situações em que ambas da medidas são desejadas.

### Exemplos

ECI de todos os países em 1998:
~~~ SQL
SELECT * FROM countries_eci(array(SELECT code FROM country), 1998);
~~~

PCI, considerando agregação geográfica em continentes e agregação temporal de 2000 a 2006:
~~~ SQL
SELECT * FROM continents_pci(array(SELECT code FROM continent), 2000, 2006);
~~~

ECI e PCI, com os países agregadados pela letra inicial do código e produtos agregados por dois pares de dígitos do HS em 2015:
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


## Função euclidean_distance

Calcula a distancia euclidiana entre duas séries de transações com diferentes agregações de produto e intervalo variável. A função retorna um valor `double precision`.

### Argumentos

#### Obrigatórios

Nome | Tipo | Descrição
---- | ---- | ---------
country_1 | text ou text[] | País 1.
country_2 | text ou text[] | País 2.

#### Opcionais

Nome | Tipo | Valor padrão | Descrição
---- | ---- | ------------ | ---------
start_year | integer | 0 | Define limite inferior para o intervalo (limite incluso).
end_year | integer | 0 | Define limite superior para o intervalo (limite incluso).
hs_code | varchar(6) | '' | Produto ou agregação seguindo o HS.

### Exemplos

Distancia euclidiana entre as séries de dois países:
~~~ SQL
SELECT euclidean_distance('chn', 'eua');
~~~

Distancia euclidiana entre as séries da um país e outros dois países do começo dos tempos até 1980:
~~~ SQL
SELECT euclidean_distance('chn', ARRAY['bra', 'chl'], end_year => 1980);
~~~

Distancia euclidiana entre as séries de dois países entre 1980 e 1999, para produtos que se enquadram na categoria café, chá, mate e pimentas:
~~~ SQL
SELECT euclidean_distance('chn', 'eua', 1980, 1999, '09');
~~~