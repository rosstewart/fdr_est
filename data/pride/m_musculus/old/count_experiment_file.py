import sys

def count_begin_ions(filename, target_file_number):
    total_lines = sum(1 for _ in open(filename, 'r'))  # Count total lines for progress tracking
    print(f"Total lines in {filename}: {total_lines}")
    
    count = 0  # Initialize the count for matching entries
    processed_lines = 0  # Keep track of processed lines for progress

    counts = {}
    scans = []

    with open(filename, 'r') as f:
        line = f.readline()
        while line:
            # Check if the line contains "BEGIN IONS"
            if line.startswith("BEGIN IONS"):
                next_line = f.readline()  # Read the next line after "BEGIN IONS"
                processed_lines += 1
                if not next_line.startswith('TITLE=id=PXD013092;20190112_LL_A.mzML;'):
                    print(next_line)
                    raise Exception

                f.readline()
                charge_line = f.readline()
                if not charge_line.startswith('CHARGE=2+'):
                    print('charge error',charge_line)
                    raise Exception

                #scan = next_line[len('TITLE=id=PXD013092;20190112_LL_A.mzML;scan='):].split()[0]
                #if scan not in scans:
                #    scans.append(scan)
                #else:
                #    print(scan,'twice error')
                #    raise Exception

                id_ = next_line.split('=')[-1].strip()
                if id_ not in counts:
                    counts[id_] = 1
                else:
                    counts[id_] += 1
                
                # Extract the file number from the TITLE line (after "file=")
                #print(next_line[-6-len(f"{target_file_number}"):].strip())
                #print(f"file={target_file_number}")
                #if next_line[-6-len(f"{target_file_number}"):].strip() == f"file={target_file_number}":
                #    count += 1

                # Print progress at every 1% completion
                if processed_lines % 50000 == 0:
                    print(f"{100*processed_lines//487000}% complete")

            # Move to the next line
            line = f.readline()
    for id_,count_ in counts.items():
        print(f"Number of 'BEGIN IONS' with file={id_}: {count_}")

# Check if the script is being run with the correct number of arguments
if len(sys.argv) != 2:
    print("Usage: python count_ions.py <filename>")
else:
    # Get the filename and file number from the command line arguments
    filename = sys.argv[1]
    file_number = -1

    # Call the function to count "BEGIN IONS"
    count_begin_ions(filename, file_number)


