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


fn_chk = sys.argv[1]
system = system = System.from_file(fn_chk).supercell(2,1,2)
#system.cell.nvec = 3
#x = np.array([1,2,2])
#system2 = system.supercell(x)
system.to_file(fn_chk[:-4]+'_super'+'.chk')