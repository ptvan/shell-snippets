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
# autodetecting column types
CREATE TABLE allergies 
AS SELECT * FROM read_csv('~/working/FHIR-sandbox/synthea_nov2021/allergies.csv');

# show columns
DESCRIBE allergies;

# order by one column and select a number of rows
SELECT * FROM allergies ORDER BY PATIENT LIMIT 3;

# delete table
DROP TABLE allergies;

## handling Parquet files
# describe the columns of a 8GB Parquet file with ~41M rows:
# (from https://datasets-documentation.s3.eu-west-3.amazonaws.com/amazon_reviews/amazon_reviews_2015.snappy.parquet)
DESCRIBE SELECT * FROM 'amazon_reviews_2015.snappy.parquet';

# aggregation is very fast
SELECT star_rating, COUNT(star_rating) FROM 'amazon_reviews_2015.snappy.parquet' GROUP BY star_rating;
