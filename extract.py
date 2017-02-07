#!/usr/bin/env python

import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


if len(sys.argv) < 3:
    print('Please provide a filename (with no extension), as follows:')
    print '> ./plotdata_v2.0.py FILENAME CHEMPOT'
    exit(0)
else:
    filename = str(sys.argv[1])
    chempot = float(sys.argv[2])
    print 'Reading', filename + '.data and', filename + '.par'

# Read the .par file to extract some parameters.
with open(filename + '.par', 'r') as f:
    lines = f.read().split('\n')
    theta = float(lines[10][7:])
    phi = float(lines[11][5:])
    if theta == 0 and phi == 0:
        spin = 'Sz'
    elif theta != 0 and phi == 0:
        spin = 'Sx'
    elif theta != 0 and phi != 0:
        spin = 'Sy'
    Temperature = lines[17][4:]
    str_list_energies = (lines[20][12:-1]).split(',')
    energies = [float(str_list_energies[i]) for i in
                range(len(str_list_energies))]
    index = energies.index(chempot)
    str_list_upos = (lines[22][8:-1]).split(',')
    upos = [float(str_list_upos[i]) for i in range(len(str_list_upos))]
    shape = [int(x) for x in lines[26][15:-1].split(',')]

# Read the array from disk. Note that this returns a 2D array.
data = np.loadtxt(filename + '.data')

# However, going back to 3D is easy if we know the original shape of the 
# array. The original shape is at the top of the .data file.
# The shape is (9, len(upos), len(energies)).
data = data.reshape((shape[0], shape[1], shape[2]))

# Store the data into individual matrices
g_ul_p = data[0, :, :]
g_ul_m = data[1, :, :]
g_ur_p = data[2, :, :]
g_ur_m = data[3, :, :]
g_rl_p = data[4, :, :]
g_rl_m = data[5, :, :]
v_p = data[6, :, :]
v_m = data[7, :, :]
dv = data[8, :, :]
dos_p = data[9, :, :]
dos_m = data[10, :, :]
ldos_p = data[11, :, :]
ldos_m = data[12, :, :]
rdos_p = data[13, :, :]
rdos_m = data[14, :, :]
udos_p = data[15, :, :]
udos_m = data[16, :, :]

# Store data in lists as a function of ulead position
g_ul_p_vs_upos = [g_ul_p[i][index] for i in range(len(g_ul_p))]
g_ul_m_vs_upos = [g_ul_m[i][index] for i in range(len(g_ul_m))]
g_ur_p_vs_upos = [g_ur_p[i][index] for i in range(len(g_ur_p))]
g_ur_m_vs_upos = [g_ur_m[i][index] for i in range(len(g_ur_m))]
g_rl_p_vs_upos = [g_rl_p[i][index] for i in range(len(g_rl_p))]
g_rl_m_vs_upos = [g_rl_m[i][index] for i in range(len(g_rl_m))]
v_p_vs_upos = [v_p[i][index] for i in range(len(v_p))]
v_m_vs_upos = [v_m[i][index] for i in range(len(v_m))]
dv_vs_upos = [dv[i][index] for i in range(len(dv))]

filename = filename + '_E' + str(chempot) + '.csv'
with open(filename, 'w') as f:
    f.write('chempot = ' + str(chempot) + '\n')
    f.write('g_ul_p, g_ul_m, g_ur_p, g_ur_m, g_rl_p, g_rl_m, v_p, v_m, dv \n')
    for i in range(len(g_ul_p_vs_upos)):
        f.write(str(g_ul_p_vs_upos[i]) + ', ' + str(g_ul_m_vs_upos[i])+ ', '
                + str(g_ur_p_vs_upos[i])+ ', ' + str(g_ur_m_vs_upos[i])+ ', '
                + str(g_rl_p_vs_upos[i])+ ', ' + str(g_rl_m_vs_upos[i])+ ', '
                + str(v_p_vs_upos[i])+ ', ' + str(v_m_vs_upos[i])+ ', '
                + str(dv_vs_upos[i]) + '\n')
f.close()