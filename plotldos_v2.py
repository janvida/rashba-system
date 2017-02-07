#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys

import kwant
import square_v2 as sq

if len(sys.argv) < 3:
    print('Please provide a filename (with no extension),'
          + ' chemical potential, and upper lead position as follows:')
    print '> ./plotldos_v2.0.py FILENAME CHEMPOT'
    exit(0)
else:
    filename = str(sys.argv[1])
    print 'Reading', filename + '.dos and', filename + '.par'
    chempot = float(sys.argv[2])

# Read the .par file to extract some parameters.
with open(filename + '.par', 'r') as f:
    lines = f.read().split('\n')
    L = int(lines[3][4:])
    W = int(lines[4][4:])
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
    shape = [3, len(energies), (L+1)*(W+1)]

# Read the array from disk. Note that this returns a 2D array.
data = np.loadtxt(filename + '.dos')

# However, going back to 3D is easy if we know the original shape of the 
# array. The original shape is at the top of the .data file.
# The shape is (9, len(upos), len(energies)).
data = data.reshape((shape[0], shape[1], shape[2]))

# Store the data into individual matrices
ldos_sys = data[0, :, :]
ldos_sys_p = data[1, :, :]
ldos_sys_m = data[2, :, :]

sys_p, leads_p = sq.make_system(L=L, W=W, srashba=0.3, lshift=4, lferro=4,
                                     theta=0, phi=0,
                                     uleftedge=1, uwidth=2)
sys_p = sys_p.finalized()
fig1 = kwant.plotter.map(sys_p, ldos_sys_p[index], num_lead_cells=3)
fig2 = kwant.plotter.map(sys_p, ldos_sys_m[index], num_lead_cells=3)
fig3 = kwant.plotter.map(sys_p, ldos_sys[index], num_lead_cells=3)

print "Saving plots to " + filename + "ldos.pdf"


pp = PdfPages(filename + '_ldos.pdf')
pp.savefig(fig1)
pp.savefig(fig2)
pp.savefig(fig3)
pp.close()