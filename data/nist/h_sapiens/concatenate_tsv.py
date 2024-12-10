import sys

# Check if two file paths are provided
if len(sys.argv) != 3:
    print("Usage: python concatenate_tsv.py <file1.tsv> <file2.tsv>")
    sys.exit(1)

# Assign file paths from command-line arguments
file1 = sys.argv[1]
file2 = sys.argv[2]
output_file = f'{sys.argv[1].replace("1","").replace(".tsv","")}_nod.tsv'

# Open the output file and concatenate the contents
with open(output_file, 'w') as outfile:
    # Write contents of the first file
    with open(file1, 'r') as f1:
        for line in f1:
            outfile.write(line)
    
    # Write contents of the second file, skipping the header
    with open(file2, 'r') as f2:
        next(f2)  # Skip the first line (header) of the second file
        for line in f2:
            outfile.write(line)

print(f"Files {file1} and {file2} have been concatenated into {output_file}")

