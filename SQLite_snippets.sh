# performance tweaks
## Allows concurrent reads while a write operation is active.
PRAGMA journal_mode = WAL;

## Speed up writes by relying on the OS to flush data to disk
PRAGMA synchronous = NORMAL;

## Enforce foreign key constraints (disabled by default in SQLite!)
PRAGMA foreign_keys = ON;

## Store temporary tables and indices in RAM instead of disk
PRAGMA temp_store = MEMORY;

# perform operations without starting an interactive session
sqlite3 my_database.db ".tables"
sqlite3 my_database.db "SELECT * FROM intervals WHERE name IN ('IGK','IGH','IGL');"


