import sys
import os

def split_mgf(input_file, output_dir, spectra_per_file=10000):
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Get the base name of the input file without the extension
    base_name = os.path.splitext(os.path.basename(input_file))[0]
    
    # Open the .mgf file
    with open(input_file, 'r') as file:
        file_count = 1
        spectra_count = 0
        output_lines = []
        
        for line in file:
            # Collect lines for the current spectrum
            output_lines.append(line)
            
            if line.strip() == 'END IONS':
                spectra_count += 1
                
                # Check if we have reached the limit for this split file
                if spectra_count >= spectra_per_file:
                    # Write to the new file with the modified name
                    output_file = os.path.join(output_dir, f"{base_name}_{file_count}.mgf")
                    with open(output_file, 'w') as out:
                        out.writelines(output_lines)
                    
                    # Reset for the next file
                    file_count += 1
                    spectra_count = 0
                    output_lines = []

        # Write any remaining spectra to the final file
        if output_lines:
            output_file = os.path.join(output_dir, f"{base_name}_{file_count}.mgf")
            with open(output_file, 'w') as out:
                out.writelines(output_lines)

if __name__ == "__main__":
    # Check that two arguments are provided
    if len(sys.argv) != 3:
        print("Usage: python split_mgf.py <input_file.mgf> <output_directory>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_dir = sys.argv[2]
    split_mgf(input_file, output_dir)

