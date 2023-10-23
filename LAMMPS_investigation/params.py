#
# Copyright 2016 James Kermode (Warwick U.)
#           2014 Lars Pastewka (U. Freiburg)
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
import numpy as np
from ase.lattice.cubic import Diamond
import ase.io
from ase.constraints import ExpCellFilter
from ase.optimize import LBFGS
from matscipy.elasticity import fit_elastic_constants
import ase.units as units

from ase.calculators.lammpslib import LAMMPSlib
import sys

################################################################################
###################THIS IS THE VALUE YOU HAVE TO CHANGE#########################

initial_K = 3.8 #supply this in MPa m^(1/2)

################################################################################
################################################################################

cmds = ['pair_style hybrid/overlay pace table spline 10000',
        'pair_coeff * * pace c_v3.yace C', 'pair_coeff * * table d2.table D2 9.0']
calc = LAMMPSlib(amendments=["mass 1 12.101"],lmpcmds = cmds,log_file='test.log',keep_alive=True)
#set calc as lammpslib object.
#calc = atomistica.TersoffScr()

# Fundamental material properties
el              = 'C'
a0_init         = 3.56

print('optimising lattice parameter')
unit_cell = Diamond(size=[1,1,1],symbol=el,latticeconstant=a0_init,pbc=(1,1,1))
unit_cell.set_calculator(calc)
ecf = ExpCellFilter(unit_cell)
uc_optimise = LBFGS(ecf)
uc_optimise.run(fmax=0.0001)
a0 = unit_cell.get_cell()[0,0] #get the optimised lattice parameter
print('optimised lattice parameter:',a0)

#surface_energy  = 5.0  * 10    # GPa*A = 0.1 J/m^2
elastic_symmetry = 'triclinic'

#fit elastic constants
C, C_err = fit_elastic_constants(unit_cell,
                                    symmetry=elastic_symmetry)

# Crack system

cleavage_plane   = np.array([ 0,-1,1 ])
crack_front     = np.array([ -2,1,1 ])
crack_direction = np.array([1,1,1])


print('Finding surface energy')
#figure out what the potential thinks the surface energy is:
from matscipy.surface import find_surface_energy
surface_energy = find_surface_energy(el,calc,a0,'diamond110',unit='0.1J/m^2')

print('SURFACE ENERGY:', surface_energy)

#supply this in angstrom
vacuum = 4
width = 300
height = 100
crack_seed_length = width/4
strain_ramp_length = width/4

initial_K_ase = (initial_K/1000) * (units.GPa*np.sqrt(units.m))
print('kase',initial_K_ase)
relax_slab=False