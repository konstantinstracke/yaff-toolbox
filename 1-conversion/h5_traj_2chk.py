import numpy as np
import h5py
from yaff import *
import sys
import molmod
from yaff.analysis.utils import *
from yaff.analysis.basic import *
from yaff.log import log
from yaff.pes.dlist import DeltaList
from yaff.pes.iclist import InternalCoordinateList
from yaff.sampling.utils import cell_lower

def map_integers_to_strings(integer_list, string_list):
    # Create a new list to hold the mapped strings
    mapped_strings = []
    for i in range(len(integer_list)):
        # Get the integer value
        integer = integer_list[i]
        # Get the corresponding string
        string = string_list[integer]
        # Append the string to the new list
        mapped_strings.append(string)
    return mapped_strings

# Must provide at least the h5 file
if len(sys.argv) < 2:
    prog = sys.argv[0]
    print(f"Usage: python {prog} init.h5 indexes")
    print("Converts a Yaff .h5 file to a .cif file for indexes or all frames if no indexes are given.\nindexes is a list of integers")
    sys.exit(1)

hfile = sys.argv[1]
f = h5py.File(hfile, "r")

# Optional indexes
if len(sys.argv) > 2:
    indexs = [int(i) for i in sys.argv[2:]]  # from argv[2] onward
else:
    indexs = list(range(f['trajectory/pos'].shape[0]))

print("Processing frames:", indexs)


system = System.from_hdf5(f)
for i in indexs:
    ffatypes_t = np.array(f['system/ffatypes'], dtype=str)
    ffatype_ids = np.array(f['system/ffatype_ids'])
    ffatypes_d = map_integers_to_strings(ffatype_ids, ['Al', 'O', 'O', 'C', 'C', 'C', 'H', 'H'])
    pos =  np.array(f['trajectory/pos'])
    rvecs = np.array(f['trajectory/cell'])
    from molmod.periodic import periodic
    numbers = np.array(f['system/numbers'])
    s = System(numbers, pos[i], ffatypes = ffatypes_d, rvecs=rvecs[i])
    outfile = str(hfile[:-3])+'_'+str(i)+'.chk'
    s.to_file(outfile)

print(len(pos))