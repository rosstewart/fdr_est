#!/bin/bash

# List of files to rename
files=(
    "Synthetic scores 1.txt"
    "Synthetic scores 4.txt"
    "Synthetic scores 8.txt"
    "Synthetic scores 10.txt"
    "Synthetic scores 5.txt"
    "Synthetic scores 9.txt"
    "Synthetic scores 2.txt"
    "Synthetic scores 6.txt"
    "Synthetic scores 3.txt"
    "Synthetic scores 7.txt"
)

# Loop through each file and rename it
for file in "${files[@]}"; do
    # Extract the number from the filename
    num=$(echo "$file" | grep -oP '\d+')

    # Define the new filename
    new_filename="synthetic_${num}.txt"

    # Rename the file
    mv "$file" "$new_filename"

    echo "Renamed '$file' to '$new_filename'"
done

