import itertools
import sys

def read_file_and_get_lines_and_numbers(filename):
    lines_and_numbers = []
    with open(filename, 'r') as file:
        for line in file:
            # Split the line and get the last element as a number
            num = int(line.split()[-1])
            lines_and_numbers.append((line.strip(), num))  # Store both the line and the number
    return lines_and_numbers

def has_consistent_tags(combo_lines):
    """Returns True if all lines in the combination have the same tag (either 'SDS' or 'Tricine')."""
    tags = {"SDS", "Tricine"}
    first_tag = None
    for line in combo_lines:
        if any(tag in line for tag in tags):
            current_tag = "SDS" if "SDS" in line else "Tricine"
            if first_tag is None:
                first_tag = current_tag
            elif current_tag != first_tag:
                return False
    return True

def find_combination_closest_to_target(lines_and_numbers, target):
    closest_sum = float('inf')  # Initialize with infinity (so any sum will be less than this initially)
    best_combination = []

    # Extract numbers for combination logic
    numbers = [num for _, num in lines_and_numbers]

    # Iterate over combinations with a maximum of 3 numbers
    for r in range(1, 4):  # r goes from 1 to 3 (inclusive)
        for combo_indices in itertools.combinations(range(len(numbers)), r):
            combo_lines = [lines_and_numbers[i][0] for i in combo_indices]
            if not has_consistent_tags(combo_lines):
                continue  # Skip this combination if tags don't match
            
            current_sum = sum(numbers[i] for i in combo_indices)
            if current_sum >= target and current_sum < closest_sum:
                closest_sum = current_sum
                best_combination = combo_indices
            # Early exit if exact match is found
            if closest_sum == target:
                return best_combination, closest_sum

    return best_combination, closest_sum

if __name__ == "__main__":
    filename = sys.argv[1]  # Take the filename as a command-line argument
    target = 37087

    # Read the file and get the list of lines and corresponding numbers
    lines_and_numbers = read_file_and_get_lines_and_numbers(filename)

    # Find the best combination that adds up to target or greater, ensuring consistent tags
    best_combination, closest_sum = find_combination_closest_to_target(lines_and_numbers, target)

    # Output the results
    if best_combination:
        print(f"Best combination sum: {closest_sum}")
        print("Lines corresponding to the best combination:")
        for index in best_combination:
            print(lines_and_numbers[index][0])  # Print the original line
    else:
        print("No valid combination found with consistent tags.")

