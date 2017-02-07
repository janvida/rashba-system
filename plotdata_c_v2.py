#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys

# Y10
energies_list = [-3.0, -2.9, -2.8, -2.7, -2.6, -2.5176380902, -2.5,
                 -2.4, -2.3, -2.2, -2.1, -2.0, -1.9, -1.8,
                 -1.7, -1.6, -1.5, -1.4823619098, -1.4, -1.3, -1.2,
                 -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5857864376,
                 -0.5, -0.4, -0.3, -0.2679491924, -0.2, -0.1,
                 -0.0681483474, 0.0, 0.0681483474, 0.1, 0.2,
                 0.2679491924, 0.3, 0.4, 0.5, 0.5857864376, 0.6,
                 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4,
                 1.4823619098, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,
                 2.1, 2.2, 2.3, 2.4, 2.5, 2.5176380902, 2.6, 2.7,
                 2.8, 2.9, 3.0]
# Y02
# energies_list = [-4.0, -3.9, -3.8, -3.7, -3.6, -3.5, -3.4142135624,
#                  -3.4, -3.3, -3.2, -3.1, -3.0, -2.9, -2.8, -2.7,
#                  -2.6, -2.5, -2.4, -2.3, -2.2, -2.1, -2.0, -1.9,
#                  -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0,
#                  -0.9, -0.8, -0.7, -0.6, -0.5857864376, -0.5, -0.4,
#                  -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5,
#                  0.5857864376, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2,
#                  1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3,
#                  2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3,
#                  3.4, 3.4142135624, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0]
# # Y04
# energies_list = [-4.0, -3.9, -3.8, -3.7320508076, -3.7, -3.6, -3.5,
#                  -3.4, -3.3, -3.2, -3.1, -3.0, -2.9, -2.8, -2.7, -2.6,
#                  -2.5, -2.4, -2.3, -2.2, -2.1, -2.0, -1.9, -1.8, -1.7,
#                  -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8,
#                  -0.7, -0.6, -0.5, -0.4, -0.3, -0.2679491924, -0.2,
#                  -0.1, 0.0, 0.1, 0.2, 0.2679491924, 0.3, 0.4, 0.5, 0.6,
#                  0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7,
#                  1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8,
#                  2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7,
#                  3.7320508076, 3.8, 3.9, 4.0]

if len(sys.argv) < 2:
    print('Please provide a filename (with no extension), as follows:')
    print '> ./plotdata_v2.0.py FILENAME'
    exit(0)
else:
    filename = str(sys.argv[1])
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
    # index = energies.index(chempot)
    str_list_upos = (lines[22][8:-1]).split(',')
    upos = [float(str_list_upos[i]) for i in range(len(str_list_upos))]
    shape = [int(x) for x in lines[26][15:-1].split(',')]

# Get the list of indices corresponding to the energies of interest
index_list = []
for energy in energies_list:
    index_list.append(energies.index(energy))

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

pp = PdfPages(filename + '_c.pdf')

for index in index_list:

    # Plot things as a function of ulead position
    g_ul_p_vs_upos = [g_ul_p[i][index] for i in range(len(g_ul_p))]
    g_ul_m_vs_upos = [g_ul_m[i][index] for i in range(len(g_ul_m))]
    g_ur_p_vs_upos = [g_ur_p[i][index] for i in range(len(g_ur_p))]
    g_ur_m_vs_upos = [g_ur_m[i][index] for i in range(len(g_ur_m))]
    g_rl_p_vs_upos = [g_rl_p[i][index] for i in range(len(g_rl_p))]
    g_rl_m_vs_upos = [g_rl_m[i][index] for i in range(len(g_rl_m))]
    v_p_vs_upos = [v_p[i][index] for i in range(len(v_p))]
    v_m_vs_upos = [v_m[i][index] for i in range(len(v_m))]
    dv_vs_upos = [dv[i][index] for i in range(len(dv))]


    # Make a 3x1 subplot of conductances.
    fig1 = plt.figure()
    fig1.suptitle('Conductances (' + spin + ', $\mu$ = ' + str(energies[index]) + 't, T = ' + Temperature + 't)', fontsize=16, fontweight='bold')
    sub1 = fig1.add_subplot(311)
    sub1.grid(True)
    sub1.set_ylabel('$G_{UL} \ [e^2/h]$')
    sub1.set_ylim((0, 1))
    sub1.plot(upos, g_ul_p_vs_upos, 'b-', label='p')
    sub1.plot(upos, g_ul_m_vs_upos, 'r-', label='m')
    # plt.legend(bbox_to_anchor=(0., 1.16, 1., .102), loc=3,ncol=3,
    #            mode="expand", borderaxespad=0.)
    plt.legend(shadow=True, fancybox=True)

    sub2 = fig1.add_subplot(312)
    sub2.grid(True)
    sub2.set_ylabel('$G_{UR} \ [e^2/h]$')
    sub2.set_ylim((0, 3))
    sub2.plot(upos, g_ur_p_vs_upos, 'b-')
    sub2.plot(upos, g_ur_m_vs_upos, 'r-')

    sub3 = fig1.add_subplot(313)
    sub3.grid(True)
    sub3.set_ylabel('$G_{RL} \ [e^2/h]$')
    sub3.set_xlabel('Upper Lead Position [a]')
    sub3.set_ylim((0, 15))
    sub3.plot(upos, g_rl_p_vs_upos, 'b-')
    sub3.plot(upos, g_rl_m_vs_upos, 'r-')

    # Fine-tune figure; make subplots close to each other and hide x ticks for
    # all but bottom plot.
    plt.setp([a.get_xticklabels() for a in fig1.axes[:-1]], visible=False)
    plt.subplots_adjust(left=0.11, bottom=0.09, right=0.96, top=0.92,
                        wspace=None, hspace=0.15)
    pp.savefig(fig1)
    plt.close('all')


print "Saving plots to " + filename + "_c.pdf"


pp.close()


