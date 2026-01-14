#! /usr/bin/env python 

from molmod.units import *
import h5py
from yaff import *
import numpy as np
import os

from yaff.system import System
from yaff.pes.ff import ForceField
from yaff.pes.eos import PREOS
from yaff.log import log
from yaff.external.lammpsio import get_lammps_ffatypes

from molmod.units import kcalmol, angstrom, kelvin, amu, pascal, deg
from molmod.constants import boltzmann
from molmod.periodic import periodic
from molmod.units import *
import sys
            
fn_xyz = sys.argv[1]
fn_chk = '/Users/konstantinstracke/work/3-MIL_volume/3-TI/4-cluster/j-Dihedral_CV/6-FF_dihedral/3-lets_go/2-fit/0-test/1-analys_output/mil53al.chk'
struct_fname = sys.argv[1].split('.')[0]

system_og = System.from_file(fn_chk)
rvecss =  np.array([[16.228734, 0.0,  0.0],
                        [ 0.0, 13.776714,  0.0],
                        [0.0, 00.0, 13.390357]])

system = System.from_file(fn_xyz)
system.ffatypes = system_og.ffatypes
system.ffatype_ids = system_og.ffatype_ids
system.bonds = system_og.bonds
#system.rvecs = rvecss
system.cell = system_og.cell

system.to_file(struct_fname+'.chk')
