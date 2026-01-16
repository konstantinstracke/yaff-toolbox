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
import sys

def dump_cif(host, fn, ffatypes=None):
    '''Write a CIF file

       **Arguments:**

       host
            A System instance

       fn
            The name of the new CIF file.

       **Optional arguments:**

       ffatypes
            A NumPy array containing atomtypes. If not give, each atom will be
            assigned a unique atom type
    '''
    if host.cell.nvec != 3:
        raise TypeError('The CIF format only supports 3D periodic systems.')
    symbols = [periodic[host.numbers[i]].symbol for i in range(host.natom)]
    if ffatypes == None:
        ffatypes = [symbols[i]+str(i+1) for i in range(host.natom)]
    assert len(ffatypes)==host.natom
    # Conversion to fractional coordinates
    frac = np.einsum('ab,ib->ia',host.cell.gvecs,host.pos)
    with open(fn, 'w') as f:
        f.write('data_\n')
        f.write('_symmetry_space_group_name_H-M       \'P1\'\n')
        f.write('_audit_creation_method            \'Yaff\'\n')
        f.write('_symmetry_Int_Tables_number       1\n')
        f.write('_symmetry_cell_setting            triclinic\n')
        f.write('loop_\n')
        f.write('_symmetry_equiv_pos_as_xyz\n')
        f.write('  x,y,z\n')
        lengths, angles = host.cell.parameters
        f.write('_cell_length_a     %12.6f\n' % (lengths[0]/angstrom))
        f.write('_cell_length_b     %12.6f\n' % (lengths[1]/angstrom))
        f.write('_cell_length_c     %12.6f\n' % (lengths[2]/angstrom))
        f.write('_cell_angle_alpha  %12.6f\n' % (angles[0]/deg))
        f.write('_cell_angle_beta   %12.6f\n' % (angles[1]/deg))
        f.write('_cell_angle_gamma  %12.6f\n' % (angles[2]/deg))
        f.write('loop_\n')
        f.write('_atom_site_label\n')
        f.write('_atom_site_type_symbol\n')
        f.write('_atom_site_fract_x\n')
        f.write('_atom_site_fract_y\n')
        f.write('_atom_site_fract_z\n')
        for i in range(host.natom):
            f.write('%10s %3s % 12.6f % 12.6f % 12.6f\n' %
             (ffatypes[i], symbols[i], frac[i,0], frac[i,1], frac[i,2]))

if len(sys.argv) != 2:
    prog = sys.argv[0]
    print(f"Usage: python {prog} init.chk")
    print("Converts a Yaff .chk file to a .cif file.")
    sys.exit(1)

fn_input = sys.argv[1]
fn_chk = fn_input
output_file = fn_input[:-4]+".cif" if len(sys.argv) == 2 else sys.argv[2] if len(sys.argv) > 2 else "output.cif"

system = System.from_file(fn_chk)
system.detect_bonds()

dump_cif(system, output_file, ffatypes=None)
