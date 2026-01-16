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
    print(f"Usage: python {prog} init.chk")
    print("Converts a Yaff .chk file to a .xyz file.")
    sys.exit(1)
    
struct_fname = sys.argv[1]
system = System.from_file(struct_fname)
fn_output = struct_fname[:-4]+".xyz" if len(sys.argv) == 2 else sys.argv[2] if len(sys.argv) > 2 else "output.xyz"
system.to_file(fn_output)