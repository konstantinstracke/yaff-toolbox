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

def dump_cif(host, pos, rvecs, fn, ffatypes=None):
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
    symbols = ffatypes#[periodic[host.numbers[i]].symbol for i in range(host.natom)]
    if ffatypes == None:
        ffatypes = [symbols[i]+str(i+1) for i in range(host.natom)]
    #print(ffatypes, host.natom)
    assert len(ffatypes)==host.natom
    #print(host.cell.gvecs, host.pos, pos, rvecs)
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

import sys
import h5py

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
    s = System(ffatype_ids, pos[i], ffatypes = ffatypes_d, rvecs=rvecs[i])
    outfile = str(hfile[:-3])+'_'+str(i)+'.cif'
    dump_cif(s, pos[i], rvecs[i], outfile, ffatypes=ffatypes_d)

print(len(pos))