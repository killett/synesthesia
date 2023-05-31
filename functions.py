import csv
import os

def load_cie_functions():
    cie = {}

    # Override the usual map output
    cie['options'] = {'output_choice': 2}

    cie['xy'] = {
        'titles': ["CIE 1931 color matching functions"],
        'x_units': ["nm"],
        'y_units': ["power"],
        'x_values': [[] for _ in range(3)],
        'y_values': [[] for _ in range(3)],
        'legends': [['x', 'y', 'z']]
    }

    file = os.path.join(".", "sealevel_spectra", "ciexyz31_1_trimmed_420nm_690nm.csv")

    # Open input file for reading
    try:
        with open(file, 'r') as in_fp:
            reader = csv.reader(in_fp)
            for row in reader:
                temp_long = float(row[0])  # Read wavelengths, record for all 3 functions
                for i in range(3):
                    cie['xy']['x_values'][i].append(temp_long)

                # Read x, y, z
                for i, val in enumerate(row[1:], start=0):
                    cie['xy']['y_values'][i].append(float(val))

    except IOError:
        print(f"The CIE file, {file}, failed to open.")

    return cie

if __name__ == "__main__":
    cie = load_cie_functions()
    print(f"{cie = }")

