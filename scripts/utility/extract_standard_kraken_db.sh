#!/bin/bash

# Set up paths
SCRATCH_DIR="/ibex/scratch/$USER"
DB_DESTINATION="$SCRATCH_DIR/kraken2_standard_db"
SOURCE_TAR="/ibex/reference/KSL/kraken2/Refseq/k2_standard_20210517.tar.gz"

# Check if the compressed file exists
if [ ! -f "$SOURCE_TAR" ]; then
    echo "Standard database tarball not found at $SOURCE_TAR"
    echo "Trying alternative location..."
    SOURCE_TAR="/ibex/reference/KSL/kraken2/Refseq/k2_standard_8gb_20210517.tar.gz"
    
    if [ ! -f "$SOURCE_TAR" ]; then
        echo "8GB standard database tarball also not found."
        echo "Please check the correct paths on your system."
        exit 1
    fi
fi

echo "Found Kraken2 standard database at: $SOURCE_TAR"
echo "Will extract to: $DB_DESTINATION"

# Create destination directory
mkdir -p $DB_DESTINATION

# Extract the database
echo "Extracting Kraken2 standard database..."
echo "This might take some time depending on the database size..."
tar -xzvf $SOURCE_TAR -C $DB_DESTINATION

# Check if extraction was successful
if [ $? -eq 0 ]; then
    echo "Database extraction successful"
    
    # List files in the database directory
    echo "Files in the database directory:"
    ls -la $DB_DESTINATION
    
    # Check for essential Kraken2 database files
    if [ -f "$DB_DESTINATION/hash.k2d" ] && [ -f "$DB_DESTINATION/opts.k2d" ] && [ -f "$DB_DESTINATION/taxo.k2d" ]; then
        echo "✓ Found all required Kraken2 database files"
        echo "Database is ready to use"
    else
        echo "⚠ Some required Kraken2 database files are missing"
        echo "The extracted database might have a different structure"
        echo "Check inside the extracted directory for the actual database location"
    fi
else
    echo "Failed to extract database"
    exit 1
fi

# If the extraction created a subdirectory, find where the actual database files are
DB_FILES_DIR=$(find $DB_DESTINATION -name "hash.k2d" -exec dirname {} \; | head -1)
if [ -n "$DB_FILES_DIR" ] && [ "$DB_FILES_DIR" != "$DB_DESTINATION" ]; then
    echo "Database files found in subdirectory: $DB_FILES_DIR"
    echo "This is the path you should use in your script"
    FINAL_DB_PATH=$DB_FILES_DIR
else
    FINAL_DB_PATH=$DB_DESTINATION
fi

echo "Final Kraken2 database path to use: $FINAL_DB_PATH"
echo "You can update your script to use this database path"

echo "export KRAKEN_DB_PATH=$FINAL_DB_PATH" > $SCRATCH_DIR/kraken_db_path.sh
echo "Database path saved to $SCRATCH_DIR/kraken_db_path.sh"
