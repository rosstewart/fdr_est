#!/bin/bash

# Check if a directory argument is provided
if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <subdirectory> [--verbose]"
    exit 1
fi

# Assign the directory from the first argument
directory="$1"
verbose=false

# Check if the verbose flag is set
if [ "$#" -eq 2 ] && [ "$2" == "--verbose" ]; then
    verbose=true
fi

total_count=0

# Find all .mgf files in the specified directory (not subdirectories)
for file in "$directory"/*.mgf; do
    # Check if the file exists (this prevents the case where no .mgf files are found)
    if [ -f "$file" ]; then
        count=$(grep -c "BEGIN IONS" "$file")
        total_count=$((total_count + count))

        # Print count for each file if verbose is true
        if $verbose; then
            echo "File: $file - BEGIN IONS count: $count"
        fi
    else
        echo "No .mgf files found in '$directory'."
    fi
done

# Output the total count
echo "Total BEGIN IONS entries in '$directory': $total_count"


# Check if a directory argument is provided
#if [ "$#" -lt 1 ]; then
#    echo "Usage: $0 <subdirectory> [--verbose]"
#    exit 1
#fi

# Assign the directory from the first argument
#directory="$1"
#verbose=false

# Check if the verbose flag is set
#if [ "$#" -eq 2 ] && [ "$2" == "--verbose" ]; then
#    verbose=true
#fi

#total_count=0

# Find all .mgf files in the specified directory and its subdirectories
#while IFS= read -r file; do
#    count=$(grep -c "BEGIN IONS" "$file")
#    total_count=$((total_count + count))
#    
#    # Print count for each file if verbose is true
#    if $verbose; then
#        echo "File: $file - BEGIN IONS count: $count"
#    fi
#done < <(find "$directory" -type f -name '*.mgf')

# Output the total count
#echo "Total BEGIN IONS entries in '$directory': $total_count"

