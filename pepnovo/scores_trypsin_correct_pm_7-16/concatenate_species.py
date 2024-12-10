import glob
import sys
import os

def sorting_key(filename):
    # Split the filename by '.' to separate the name from the extension
    # Then split the name by '_' and get the last element
    return int(filename.split('.')[0].split('_')[-1])

def concatenate_files(input_pattern, output_file):
    # Get a sorted list of all files matching the input pattern
    file_list = sorted(glob.glob(input_pattern),key=sorting_key)
    
    # To keep track of the last index for the '>>' lines
    current_index = 0
    
    with open(output_file, 'w') as outfile:
        for filename in file_list:
            print('writing',filename)
            with open(filename, 'r') as infile:
                for line in infile:
                    if line.startswith('>>'):
                        # Increment the third element
                        parts = line.split()
                        if len(parts) >= 3:
                            parts[2] = str(current_index)  # Update the third element
                            current_index += 1  # Increment for the next '>>' line
                        outfile.write(' '.join(parts) + '\n')  # Write the modified line
                    else:
                        # Write all other lines unchanged
                        outfile.write(line)

if __name__ == "__main__":
    # Expect two arguments: input pattern and output file
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_pattern> <output_file>")
        sys.exit(1)

    input_pattern = sys.argv[1] + '_*.txt'
    output_file = sys.argv[2]

    concatenate_files(input_pattern, output_file)
    print(f"Files concatenated into {output_file}")

