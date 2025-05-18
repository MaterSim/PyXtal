from pyxtal.symmetry import Group
import csv

def csv_to_rst(csv_path, rst_path):
    with open(csv_path, 'r') as csvfile, open(rst_path, 'w') as rstfile:
        # Write RST header
        rstfile.write("Wyckoff Positions\n")
        rstfile.write("=================\n\n")

        # Read CSV
        reader = csv.reader(csvfile)
        headers = next(reader)

        # Write table header
        rstfile.write(".. list-table::\n")
        rstfile.write("   :header-rows: 1\n")
        rstfile.write("   :widths: auto\n\n")

        # Write headers
        rstfile.write("   * - " + "\n     - ".join(headers) + "\n")

        # Write data rows
        for row in reader:
            rstfile.write("   * - " + "\n     - ".join(row) + "\n")



# Open a CSV file for writing
with open('doc/wyckoff_positions.csv', 'w', newline='') as csvfile:
    # Create CSV writer
    csvwriter = csv.writer(csvfile)

    # Write header
    csvwriter.writerow(['Space Group Number', 'Space Group Symbol', 'Wyckoff Label', 'Site Symmetry'])

    # Iterate through space groups
    for g in range(1, 231):
        spg = Group(g)
        for wp in spg:
            wp.get_site_symmetry()
            # Write data row
            csvwriter.writerow([spg.number, spg.symbol, wp.get_label(), wp.site_symm])

# Convert CSV to RST
csv_to_rst('doc/wyckoff_positions.csv', 'doc/wyckoff_positions.rst')
