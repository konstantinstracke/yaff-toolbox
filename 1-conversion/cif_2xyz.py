#!/usr/bin/env python

from ase.io import read, write

# Input/output files
import sys
input_file = sys.argv[1] if len(sys.argv) > 1 else "init.cif"
output_file = input_file[:-4]+".xyz" if len(sys.argv) == 2 else sys.argv[2] if len(sys.argv) > 2 else "output.xyz"

try:
    # Read CIF using ASE
    atoms = read(input_file)

    # Prepare XYZ second line with cell vectors
    cell = atoms.cell
    cell_flat = cell.reshape(-1)  # Flatten to 9 values
    comment_line = "cell= " + ' '.join(f"{val:.6f}" for val in cell_flat)

    # Write XYZ with cell info as comment
    write(output_file, atoms, format='xyz', comment=comment_line)

    print(f"✅ Converted {input_file} → {output_file}")
except Exception as e:
    print(f"❌ Error: {e}")
