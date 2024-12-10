import os

# Input and output file paths
input_file = 'h_sapiens.mgf'
output_file1 = 'h_sapiens_part1.mgf'
output_file2 = 'h_sapiens_part2.mgf'

# Function to count BEGIN IONS blocks
def count_blocks(filename):
    with open(filename, 'r') as f:
        return sum(1 for line in f if line.strip() == 'BEGIN IONS')

# Split the file after a certain number of blocks
def split_mgf(filename, output1, output2, block_count):
    with open(filename, 'r') as infile, open(output1, 'w') as out1, open(output2, 'w') as out2:
        current_block = 0
        out = out1
        
        for line in infile:
            # Switch output files after the desired number of blocks
            if line.strip() == 'BEGIN IONS':
                current_block += 1
                if current_block > block_count:
                    out = out2
            
            out.write(line)

# Count total BEGIN IONS blocks
total_blocks = count_blocks(input_file)

# Split the file after half of the blocks
split_mgf(input_file, output_file1, output_file2, total_blocks // 2)

print(f"File has been split into {output_file1} and {output_file2}")

