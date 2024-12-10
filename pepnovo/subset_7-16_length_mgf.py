import os
import sys
import re

def filter_mgf_files(input_dir, output_dir):
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Regular expression to capture sequence in the TITLE line
    title_pattern = re.compile(r'^TITLE=(\w+)/')

    for filename in os.listdir(input_dir):
        if filename.endswith('.mgf'):
            input_path = os.path.join(input_dir, filename)
            output_path = os.path.join(output_dir, filename)

            with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
                write_entry = False
                entry = []  # Buffer for each spectrum entry

                for line in infile:
                    # Start of a new spectrum entry
                    if line.startswith('BEGIN IONS'):
                        entry = [line]
                        write_entry = False  # Reset the flag for each entry

                    elif line.startswith('TITLE='):
                        entry.append(line)
                        match = title_pattern.search(line)
                        if match:
                            sequence = match.group(1)
                            if 7 <= len(sequence) <= 16:
                                write_entry = True  # Mark for writing if sequence length is valid

                    # End of a spectrum entry
                    elif line.startswith('END IONS'):
                        entry.append(line)
                        if write_entry:
                            outfile.writelines(entry)  # Write the entire entry if flagged

                    # Middle of a spectrum entry
                    else:
                        entry.append(line)

    print(f"Filtered .mgf files saved to {output_dir}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_directory> <output_directory>")
        sys.exit(1)

    input_dir = sys.argv[1]
    output_dir = sys.argv[2]
    filter_mgf_files(input_dir, output_dir)

