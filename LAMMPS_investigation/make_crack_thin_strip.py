#
# Copyright 2014, 2020 James Kermode (Warwick U.)
#           2016 Punit Patel (Warwick U.)
#
# matscipy - Materials science with Python at the atomic-scale
# https://github.com/libAtoms/matscipy
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

# ======================================================================
# matscipy - Python materials science tools
# https://github.com/libAtoms/matscipy
#
# Copyright (2014) James Kermode, King's College London
#                  Lars Pastewka, Karlsruhe Institute of Technology
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# ======================================================================

"""
Script to generate a crack slab, and apply initial strain ramp

James Kermode <james.kermode@kcl.ac.uk>
August 2014
"""

import numpy as np

from ase.lattice.cubic import Diamond
from ase.constraints import FixAtoms, StrainFilter
from ase.optimize import FIRE
import ase.io
import ase.io.lammpsdata
import ase.units as units

from matscipy.elasticity import (measure_triclinic_elastic_constants,
                                 rotate_elastic_constants,
                                 youngs_modulus, poisson_ratio)

from matscipy.fracture_mechanics.crack import (print_crack_system,
                                               G_to_strain,
                                               thin_strip_displacement_y,
                                               find_tip_stress_field)
import sys
sys.path.insert(0, '.')
import params
from matscipy import parameter
import sys



# ***** Find eqm. lattice constant ******

a0 = parameter('a0')
# ******* Find elastic constants *******
C = parameter('C')
# Get 6x6 matrix of elastic constants C_ij
print('Elastic constants (GPa):')
print((C / units.GPa).round(0))
print('')

E = youngs_modulus(C, parameter('cleavage_plane'))
print('Young\'s modulus %.1f GPa' % (E / units.GPa))
nu = poisson_ratio(C, parameter('cleavage_plane'), parameter('crack_direction'))
print('Poisson ratio % .3f\n' % nu)

# **** Setup crack slab unit cell ******

directions = [parameter('crack_direction'),
              parameter('cleavage_plane'),
              parameter('crack_front')]
print_crack_system(directions)

# now, we build system aligned with requested crystallographic orientation
unit_slab = Diamond(directions=directions,
                    size=(1, 1, 1),
                    symbol='C',
                    pbc=True,
                    latticeconstant=a0)



print('Unit slab with %d atoms per unit cell:' % len(unit_slab))
print(unit_slab.cell)
print('')

# center vertically half way along the vertical bond between atoms 0 and 1
unit_slab.positions[:, 1] += (unit_slab.positions[1, 1] -
                              unit_slab.positions[0, 1]) / 2.0

# map positions back into unit cell
unit_slab.set_scaled_positions(unit_slab.get_scaled_positions())

# ***** Setup crack slab supercell *****

# Now we will build the full crack slab system,
# approximately matching requested width and height
nx = int(parameter('width') / unit_slab.cell[0, 0])
ny = int(parameter('height') / unit_slab.cell[1, 1])

# make sure ny is even so slab is centered on a bond
if ny % 2 == 1:
    ny += 1

# make a supercell of unit_slab
crack_slab = unit_slab * (nx, ny, 1)
crack_slab.positions += [1,0.1,0]
crack_slab.wrap()

# open up the cell along x and y by introducing some vaccum
crack_slab.center(parameter('vacuum'), axis=0)
crack_slab.center(parameter('vacuum'), axis=1)

# centre the slab on the origin
xmean = crack_slab.positions[:, 0].mean()
ymean = crack_slab.positions[:, 1].mean()
crack_slab.positions[:, 0] -= xmean
crack_slab.positions[:, 1] -= ymean

orig_width = (crack_slab.positions[:, 0].max() -
              crack_slab.positions[:, 0].min())
orig_height = (crack_slab.positions[:, 1].max() -
               crack_slab.positions[:, 1].min())

print(('Made slab with %d atoms, original width and height: %.1f x %.1f A^2' %
       (len(crack_slab), orig_width, orig_height)))

top = crack_slab.positions[:, 1].max()
bottom = crack_slab.positions[:, 1].min()
left = crack_slab.positions[:, 0].min()
right = crack_slab.positions[:, 0].max()

# fix atoms in the top and bottom rows
fixed_mask = ((abs(crack_slab.positions[:, 1] - top) < 1.0) |
              (abs(crack_slab.positions[:, 1] - bottom) < 1.0))
const = FixAtoms(mask=fixed_mask)
crack_slab.set_constraint(const)
print('Fixed %d atoms\n' % fixed_mask.sum())

# Save all calculated materials properties inside the Atoms object
crack_slab.info['nneightol'] = 1.3 # nearest neighbour tolerance
crack_slab.info['LatticeConstant'] = a0
crack_slab.info['C11'] = C[0, 0]
crack_slab.info['C12'] = C[0, 1]
crack_slab.info['C44'] = C[3, 3]
crack_slab.info['YoungsModulus'] = E
crack_slab.info['PoissonRatio_yx'] = nu
crack_slab.info['SurfaceEnergy'] = parameter('surface_energy')
crack_slab.info['OrigWidth'] = orig_width
crack_slab.info['OrigHeight'] = orig_height
crack_slab.info['CrackDirection'] = parameter('crack_direction')
crack_slab.info['CleavagePlane'] = parameter('cleavage_plane')
crack_slab.info['CrackFront'] = parameter('crack_front')
crack_slab.info['cell_origin'] = -np.diag(crack_slab.cell)/2.0

crack_slab.set_array('fixed_mask', fixed_mask)
ase.io.write('slab.xyz', crack_slab, format='extxyz')

# ****** Apply initial strain ramp *****
initial_G = (parameter('initial_K_ase')**2)/E
print('E is', E)
strain = G_to_strain(initial_G, E, nu, orig_height)

print(strain)
crack_slab.positions[:, 1] += thin_strip_displacement_y(
                                 crack_slab.positions[:, 0],
                                 crack_slab.positions[:, 1],
                                 strain,
                                 left + parameter('crack_seed_length'),
                                 left + parameter('crack_seed_length') +
                                        parameter('strain_ramp_length'))


print('Applied initial load: strain=%.4f, G=%.2f J/m^2' %
      (strain, initial_G / (units.J / units.m**2)))


# ***** Relaxation of crack slab  *****
calc = parameter('calc')

# optionally, relax the slab, keeping top and bottom rows fixed
if parameter('relax_slab',False):
    print('Relaxing slab...')
    crack_slab.set_calculator(calc)
    opt = FIRE(crack_slab)
    opt.run(fmax=0.005)

# Find initial position of crack tip
#crack_pos = find_tip_stress_field(crack_slab, calc=calc)
#print('Found crack tip at position %s' % crack_pos)

#crack_slab.info['strain'] = strain
#crack_slab.info['G'] = initial_G
#crack_slab.info['CrackPos'] = crack_pos
#crack_slab.wrap()
# ******** Save output file **********

# undo crack shift
crack_slab.positions[:, 0] += xmean
crack_slab.positions[:, 1] += ymean

# Save results in extended XYZ format, including extra properties and info
print('Writing crack slab to file "crack.xyz"')
ase.io.write('crack.xyz', crack_slab, format='extxyz')
ase.io.lammpsdata.write_lammps_data('crack.lj',crack_slab)