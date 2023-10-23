import ase
import numpy as np
from matscipy.fracture_mechanics.crack import find_tip_coordination
from matscipy.fracture_mechanics.clusters import set_groups
mid_point_traj = ase.io.read('dump.lammpstrj',index='5')
final_traj = ase.io.read('dump.lammpstrj',index='10')

def get_tip_position(traj):
    n = [300,100,1]
    mid_point_traj.set_pbc([False,False,True])
    set_groups(traj,n,20,10)
    bond_atoms = find_tip_coordination(traj,bondlength=1.9,bulk_nn=4)
    ase.io.write(f'traj.xyz',traj)
    tip_pos = (traj.get_positions()[bond_atoms,:][:,0])
    #print(tip_pos[0],tip_pos[1])
    assert tip_pos[0]-tip_pos[1] < 1
    return np.mean(tip_pos)


mid_point_tip = get_tip_position(mid_point_traj)
final_tip = get_tip_position(final_traj)
crack_speed = (final_tip-mid_point_tip)/(5*0.1)

print(f'crack speed measured as {crack_speed} A/ps')
print(f'which is {crack_speed*100} m/s')
