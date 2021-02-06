#!/usr/bin/python3
import os
import numpy as np
import scipy.io
from io import StringIO
import natsort

# Script to extract SHARC Trajectories into .mat arrays for scattering calcs
# Assumes $SHARC/geo.py has generated Geo.out to analyse internal coords
# Assumes trajectory directories are named "TRAJ_XXXXX" - Default convention.
# OUTPUT - Gives a geometry array with dimensions -
# (n_atom, xyz, run_time/dt, number_trajs)
# diss_marker - dimensions (1, number_trajs) - 1 or 0 for dissociation or not
# diss_times - times of when diss threshold is reached - 0 for non-diss trajs
# cslabels - 0, 2 or 3 - depending on if no-dissociation, atom 2 (1st Sulfur
# dissociated) or atom 3 (2nd Sulfur dissociated)
# dissociated_multiplicities - 0 for non-diss trajs, 1 for singlet, 3 for
# triplet diss trajs
# final_t - dimensions (1, number_trajs) - stores final times of all trajs
# DISS_CHECK.txt - text file listing what trajs dissociated after time t

# Options to EDIT
n_atom = 3
diss_limit = float(3.4)
run_time = 1000
dt = 0.5
#######################
cwd = os.getcwd()
directories_list_full = [d for d in os.listdir(os.getcwd())
                         if os.path.isdir(d) and "TRAJ" in d]
directories_list = [d for d in directories_list_full if
                    os.path.isfile(d+"/output.xyz")]
directories_list = [d for d in directories_list if
                    os.path.getsize(d+"/output.xyz") != 0]
directories_list = natsort.natsorted(directories_list)  # order directories
number_trajs = len(directories_list)
number_geoms = int(run_time/dt)
geom_array_dim = (n_atom, 3, number_trajs, number_geoms+1)
geom_array = np.zeros(geom_array_dim)  # Contains xyz coords over time and traj
diss_marker = np.zeros(number_trajs)  # 1 for dissociated and 0 for not
final_t = np.zeros(number_trajs)  # Final times for all trajs
diss_time = np.zeros(number_trajs)
S_values = np.zeros(number_trajs)
CS_labels = np.zeros(number_trajs)
diss_counter = 0
with open('DISS_CHECK.txt', 'a') as out_file:
    out_file.write('Geometry Array has dimensions %s (atom no.), %s (x,y,z)\
%s (traj number), %s (time step)\n' % (n_atom, '3', number_trajs, number_geoms))
for traj_index in range(0, number_trajs):
    traj_dir = cwd+"/"+directories_list[traj_index]+"/"
    traj_num = directories_list[traj_index].rsplit('TRAJ_')[1]
    read = True
    with open('DISS_CHECK.txt', 'a') as out_file:
        with open('%sGeo.out' % (traj_dir), 'r') as geo_file:
            for idx, line in enumerate(geo_file):
                if read == True and idx > 1:
                    data = np.genfromtxt(StringIO(line))
                    if data[1] >= diss_limit or data[2] >= diss_limit:
                        if data[1] >= diss_limit:
                            CS_labels[traj_index] = 2
                        elif data[2] >= diss_limit:
                            CS_labels[traj_index] = 3
                        with open('%soutput.lis' % (traj_dir), 'r') as out_f:
                            for line_out in out_f:
                                pass
                            last_line_outf = line_out
                            spin_data = np.genfromtxt(StringIO(last_line_outf))
                        if spin_data[9] >= 1:
                            S_values[traj_index] = 3
                            spin = 'TRIPLET, S=3'
                        elif spin_data[9] < 1:
                            S_values[traj_index] = 1
                            spin = 'SINGLET, S=1'
                        diss_counter += 1
                        diss_marker[traj_index] = 1
                        diss_time[traj_index] = data[0]
                        time = str(data[0])
                        out_file.write('TRAJ %s (Index %s) DISSOCIATED after %s fs. %s product\n'
                                        % (traj_num, traj_index+1, time, spin))
                        read = False
                    pass
            last_line = line
            data = np.genfromtxt(StringIO(last_line))
        if data[1] >= diss_limit or data[2] >= diss_limit:  # Second column of Geo.out - C-S1 len
            final_t[traj_index] = data[0]
        else:
            diss_marker[traj_index] = 0
            final_t[traj_index] = data[0]
            out_file.write('TRAJ %s (Index %s) DID NOT dissociate.\n'
                           % (traj_num, traj_index+1))
            pass
        with open('%soutput.xyz' % (traj_dir), 'r') as xyz_file:
            lines = xyz_file.readlines()
            count = 0
            repeating_geom_lines = n_atom+2  # In output.xyz geometry
    # block at t=t+dt repeates at every n_atom+2 lines
            for i in range(0, len(lines), repeating_geom_lines):
                geom_section = lines[i:i+repeating_geom_lines]
                geom = geom_section[2:repeating_geom_lines]
                geom = [i[1:] for i in geom]  # Remove atom labels
                x = ''.join(geom)
                b = StringIO(x)
                geom_array[:, :, traj_index, count] = np.loadtxt(b)
                count += 1
        scipy.io.savemat('geometries.mat',
                         mdict={'geometries': geom_array})
        scipy.io.savemat('final_times.mat',  mdict={'final_t': final_t})
        scipy.io.savemat('diss_times.mat', mdict={'diss_times': diss_time})
        scipy.io.savemat('dissociation_marker.mat',
                         mdict={'diss_marker': diss_marker})
        scipy.io.savemat('dissociated_multiplicities.mat',
                         mdict={'dissociation_S': S_values})
        scipy.io.savemat('cslabels.mat',
                         mdict={'cslabels': CS_labels})
print(np.size(geom_array))
print('%s Total trajectories extracted - %s dissociated, %s did not dissociate\
, %s out of %s trajectories were unusuable' % (str(number_trajs),
str(diss_counter), str(number_trajs-diss_counter),
str(len(directories_list_full)-number_trajs), str(len(directories_list_full))))
