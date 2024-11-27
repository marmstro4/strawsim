
def remove_third_column_duplicates(input_file, output_file):
    try:
        with open(input_file, 'r') as infile:
            lines = infile.readlines()

        # Dictionary to store unique third column values
        unique_third_col = {}
        for line in lines:
            columns = line.strip().split(',')
            if len(columns) != 3:
                print(f"Skipping malformed line: {line.strip()}")
                continue

            third_col = columns[2]  # Third column value
            if third_col not in unique_third_col:
                unique_third_col[third_col] = line

        # Write only lines with unique third column values
        with open(output_file, 'w') as outfile:
            outfile.writelines(unique_third_col.values())

        print(f"Lines with unique third column values written to '{output_file}'.")
    except FileNotFoundError:
        print(f"Error: File '{input_file}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

# Usage example
if __name__ == "__main__":
    input_filename = "build/rows.csv"  # Replace with your input file
    output_filename = "build/vertex.csv"  # Replace with your desired output file
    remove_third_column_duplicates(input_filename, output_filename)
