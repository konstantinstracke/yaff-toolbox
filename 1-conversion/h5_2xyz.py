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

if len(sys.argv) != 2:
    prog = sys.argv[0]
    print(f"Usage: python {prog} init.h5")
    print("Converts a Yaff .h5 file to a .xyz file.")
    sys.exit(1)

hfile = sys.argv[1]
f = h5py.File(hfile)

system = System.from_hdf5(f)
system.to_file(str(hfile)+'.xyz')