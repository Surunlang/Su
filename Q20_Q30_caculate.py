import os
import re
import argparse

def extract_q20_q30_after_filtering(file_path):
    q20_values = []
    q30_values = []
    with open(file_path, 'r') as file:
        read_section = None
        for line in file:
            if "Read1 after filtering:" in line:
                read_section = "Read1"
                continue
            elif "Read2 after filtering:" in line:
                read_section = "Read2"
                continue
            elif read_section and "Q20 bases" in line:
                q20_match = re.search(r"Q20 bases: \d+\((\d+\.\d+)%\)", line)
                if q20_match:
                    q20_value = float(q20_match.group(1))
                    q20_values.append(q20_value)
            elif read_section and "Q30 bases" in line:
                q30_match = re.search(r"Q30 bases: \d+\((\d+\.\d+)%\)", line)
                if q30_match:
                    q30_value = float(q30_match.group(1))
                    q30_values.append(q30_value)
            elif read_section and "Filtering result:" in line:
                read_section = None
    return q20_values, q30_values

def calculate_average_q20_q30(log_dir):
    all_q20_values = []
    all_q30_values = []
    
    for file_name in os.listdir(log_dir):
        if file_name.endswith(".log"):  # Adjust the file extension if necessary
            file_path = os.path.join(log_dir, file_name)
            q20_values, q30_values = extract_q20_q30_after_filtering(file_path)
            print(f"Processing {file_name}: Q20={q20_values}, Q30={q30_values}")  # Debug info
            all_q20_values.extend(q20_values)
            all_q30_values.extend(q30_values)

    if all_q20_values and all_q30_values:
        average_q20 = sum(all_q20_values) / len(all_q20_values)
        average_q30 = sum(all_q30_values) / len(all_q30_values)
        return average_q20, average_q30
    else:
        return None, None

def main():
    parser = argparse.ArgumentParser(description="Calculate average Q20 and Q30 values from fastp log files.")
    parser.add_argument('-i', '--input_dir', required=True, help="Directory containing fastp log files")
    
    args = parser.parse_args()
    log_directory = args.input_dir

    average_q20, average_q30 = calculate_average_q20_q30(log_directory)

    if average_q20 is not None and average_q30 is not None:
        print(f"Average Q20: {average_q20:.2f}%")
        print(f"Average Q30: {average_q30:.2f}%")
    else:
        print("No Q20 or Q30 values found in the log files.")

if __name__ == "__main__":
    main()

