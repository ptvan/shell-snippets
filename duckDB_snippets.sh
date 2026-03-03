# read in multiple CSVs of the same format as a single transient database
# eg. file1.csv, file2.csv all with the same columns
# NOTE: this _DOES NOT_ create a new table
SELECT * FROM 'file*.csv';

# attach a persistent database on disk (will be created if it doesn't exist)
ATTACH 'persistent_duck.db';

# copy in-memory database into attached database
COPY FROM DATABASE memory TO persistent_duck;

# create a new persistent database
duckdb ~/working/FHIR-sandbox/synthea_nov2021.db

# launch a web-based interactive UI on `localhost:4213`
duckdb -ui

# list databases
.databases

# list tables
.tables

# read in and execute SQL from external file
.read import_data.sql

# show current database
SELECT current_database();

# switch databases
USE other_database;

# detach database(s):
# NOTE: use only name, no .db suffix
DETACH persistent_duck;

# read in a CSV containing Synthea FHIR data (https://synthea.mitre.org/downloads/, '1K Sample Synthetic Patient Records, CSV'
# autodetecting column types and creating a new table
CREATE TABLE allergies AS SELECT * FROM read_csv('~/working/FHIR-sandbox/synthea_nov2021_CSV/allergies.csv');

# show columns
DESCRIBE allergies;

# order by one column and select a number of rows
SELECT * FROM allergies ORDER BY PATIENT LIMIT 3;

# !DuckDB-specific SQL: exclude specific columns
SELECT * EXCLUDE ('REACTION2', 'DESCRIPTION2', 'SEVERITY2') FROM allergies;

# !DuckDB-specific SQL: select columns by regex pattern
SELECT COLUMNS('SEVERITY*') FROM allergies;

# delete table
DROP TABLE allergies;

## handling Parquet files
# describe the columns of a 8GB Parquet file with ~41M rows:
# (from https://datasets-documentation.s3.eu-west-3.amazonaws.com/amazon_reviews/amazon_reviews_2015.snappy.parquet)
DESCRIBE SELECT * FROM 'amazon_reviews_2015.snappy.parquet';

# lower-level column information
SELECT * FROM parquet_schema('amazon_reviews_2015.snappy.parquet');

# built-in functions
SELECT MAX(total_votes) FROM 'amazon_reviews_2015.snappy.parquet';

## aggregation 
# breaking down star rating from the same dataset above is instantaneous
duckdb -c "SELECT star_rating, COUNT(star_rating) FROM 'amazon_reviews_2015.snappy.parquet' GROUP BY star_rating;"

## views
# list views
FROM duckdb_views WHERE NOT internal;

# duckdb views are virtual, results are re-calculated on demand
CREATE VIEW by_star_verified AS SELECT * FROM 'amazon_reviews_2015.snappy.parquet' WHERE verified_purchase = true ORDER BY star_rating DESC;

# querying views is the same as querying an actual table
SELECT product_id, star_rating, helpful_votes, total_votes FROM by_star_verified ORDER BY helpful_votes DESC;

# drop views
DROP VIEW by_star_verified;

# converting from CSV to Parquet, autodetecting columns and renaming fields
COPY (SELECT DATE AS ENCOUNTER_DATE, PATIENT AS PATIENT_ID, ENCOUNTER AS ENCOUNTER_ID 
	    FROM read_csv('~/working/FHIR-sandbox/synthea_nov2021_CSV/observations.csv', AUTO_DETECT=TRUE))
  TO '~/working/FHIR-sandbox/observations_recoded.parquet' (COMPRESSION zstd);