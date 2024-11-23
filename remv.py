import pandas as pd

def remove_duplicates(input_file="build/output.csv", output_file="build/vertex.csv"):
    try:
        # Read the input file
        data = pd.read_csv(input_file, header=None, names=["Value"])

        # Drop duplicate values
        unique_data = data.drop_duplicates()

        # Save the unique values to a new file
        unique_data.to_csv(output_file, index=False, header=False)

        print(f"Duplicates removed. Unique values saved to '{output_file}'.")
    except Exception as e:
        print(f"An error occurred: {e}")

# Call the function
remove_duplicates()
