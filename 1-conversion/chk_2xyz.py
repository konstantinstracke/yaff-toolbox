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
    print("Usage : python chk_2xyz.py <struct_file>")
    sys.exit(1)

struct_fname = sys.argv[1]

system = System.from_file(struct_fname)
system.to_file(struct_fname[:-4]+'.xyz')