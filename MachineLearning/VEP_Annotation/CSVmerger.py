#this script merges all the CSV files in a folder into a single CSV file
import glob
import csv
import os
import sys

def merge_csv_files(folder, output_file='merged.csv'):
    """
    Merges all CSV files in the specified folder into a single CSV file.
    The output file will be named 'merged.csv' as default and the produced file is placed in the batch folder.
    """
    csv_files = sorted(glob.glob(os.path.join(folder, '*.csv')))
    output_file = os.path.join(folder, output_file)
    # Exclude the output file if it already exists in the folder
    csv_files = [f for f in csv_files if os.path.abspath(f) != os.path.abspath(output_file)]

    if not csv_files:
        print("No CSV files found!")
        sys.exit(1)

    with open(csv_files[0], 'r', newline='', encoding='utf-8') as firstfile:
        reader = csv.reader(firstfile)
        try:
            header = next(reader)
        except StopIteration:
            print(f"First file {csv_files[0]} is empty.")
            sys.exit(1)

    with open(output_file, 'w', newline='', encoding='utf-8') as outfile:
        writer = csv.writer(outfile)
        writer.writerow(header)
        for filename in csv_files:
            with open(filename, 'r', newline='', encoding='utf-8') as infile:
                reader = csv.reader(infile)
                try:
                    file_header = next(reader)
                except StopIteration:
                    print(f"File {filename} is empty, skipping.")
                    continue
                if file_header != header:
                    print(f"Header mismatch in file: {filename}\nExpected: {header}\nFound:    {file_header}")
                    sys.exit(1)
                for row in reader:
                    writer.writerow(row)

    print(f"Merged {len(csv_files)} files into {output_file}")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        target_folder = sys.argv[1]
    else:
        target_folder = '.'
    merge_csv_files(target_folder)
    print("CSV files merged successfully.")