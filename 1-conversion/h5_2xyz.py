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

hfile = sys.argv[1]
f = h5py.File(hfile)

system = System.from_hdf5(f)
system.to_file(str(hfile)+'.xyz')