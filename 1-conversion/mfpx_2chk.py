from yaff import System, angstrom
import numpy as np
import sys

def preprocess_mfpx(infile, outfile):
    """
    Remove neighbor columns from an .mfpx file.
    Keeps only:
    id, element, x, y, z, ffatype, group
    """
    with open(infile, 'r') as fin, open(outfile, 'w') as fout:
        lines = fin.readlines()

        # Copy header
        fout.write(lines[0])
        fout.write(lines[1])

        for line in lines[2:]:
            if not line.strip():
                continue
            parts = line.split()
            trimmed = parts[:7]
            fout.write("{:>5s} {:<2s} {:>12s} {:>12s} {:>12s} {:<15s} {:<10s}\n"
                       .format(*trimmed))
# ---- main script ----

if len(sys.argv) != 3:
    prog = sys.argv[0]
    print(f"Usage: python {prog} init.mfpx init.xyz")
    print("Converts a Yaff .mfpx AND .xyz file to a .chk file.\nNeeds rvecs hardcoded.")
    sys.exit(1)

init_mfpx = sys.argv[1]
xyz_file = sys.argv[2]
clean_mfpx = init_mfpx[:-5] + "_clean.mfpx"

preprocess_mfpx(init_mfpx, clean_mfpx)

atomtypes = np.genfromtxt(clean_mfpx, skip_header=2, dtype=str)[:, -2:]
ffatypes = [atom[0] + "_" + atom[1] for atom in atomtypes]

a = [20.92565778, 0., 0.]
b = [1.28132699e-15, 2.09256578e+01, 0.]
c = [1.28132699e-15, 1.28132699e-15, 2.09256578e+01]
rvecs = np.array([a, b, c])


print(xyz_file)

system = System.from_file(xyz_file, ffatypes=ffatypes, rvecs=rvecs * angstrom)
print(system)

system.detect_bonds()
system.to_file(clean_mfpx[:-5] + '.chk')
