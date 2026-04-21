#! /usr/bin/env python
import sys
import os
import numpy as np
from molmod.units import *
from molmod.units import kcalmol, angstrom, kelvin, amu, pascal, deg
from molmod.constants import boltzmann
from molmod.periodic import periodic
from yaff import *
from yaff.system import System
from yaff.pes.ff import ForceField
from yaff.log import log

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
            
# --- input arguments ---
if len(sys.argv) < 4:
    print("Usage: python run_yaff_opt.py structure.chk pars.txt geo|cell [output.chk|output.cif|output.xyz]")
    sys.exit(1)

fn_chk   = sys.argv[1]
fn_ff    = sys.argv[2]
opt_type = sys.argv[3]          # 'geo' or 'cell'
fn_out   = sys.argv[4] if len(sys.argv) >= 5 else fn_chk[:-4] + '.chk'

# --- setup ---
system = System.from_file(fn_chk)
ff = ForceField.generate(system, fn_ff, rcut=12*angstrom, smooth_ei=False, gcut_scale=1.5, alpha_scale=3.2)

if opt_type == 'geo':
    dof = CartesianDOF(ff, gpos_rms=1e-8, dpos_rms=1e-6)
elif opt_type == 'cell':
    dof = StrainCellDOF(ff)
else:
    raise ValueError("opt_type must be 'geo' or 'cell', got: " + opt_type)

# --- Hooks ---
screen = OptScreenLog(step=1)
opt = CGOptimizer(dof, hooks=[screen])

# --- Run ---
try:
    opt.run()
except KeyboardInterrupt:
    log("INTERRUPTED")
    log.hline()

# --- Save output ---
if fn_out.endswith('.cif'):
    dump_cif(system, fn_out, ffatypes=None)
elif fn_out.endswith('.xyz'):
    system.to_file(fn_out)
else:
    system.to_file(fn_out)

print("Output written to:", fn_out)