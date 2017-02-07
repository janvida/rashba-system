#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys

chempot = -2.5179
u = 100

if len(sys.argv) < 4:
    print('Please provide a filename (with no extension),'
          + ' chemical potential, and upper lead position as follows:')
    print '> ./plotdata_v2.0.py FILENAME CHEMPOT UPOS'
    exit(0)
else:
    filename = str(sys.argv[1])
    print 'Reading', filename + '.data and', filename + '.par'
    chempot = float(sys.argv[2])
    u = int(sys.argv[3])

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

# Plot as a function of energy, to see the effect of temperature smoothing
u = upos.index(u)
g_ul_p_vs_energy = g_ul_p[u][:]
g_ul_m_vs_energy = g_ul_m[u][:]
g_ur_p_vs_energy = g_ur_p[u][:]
g_ur_m_vs_energy = g_ur_m[u][:]
g_rl_p_vs_energy = g_rl_p[u][:]
g_rl_m_vs_energy = g_rl_m[u][:]
v_p_vs_energy = v_p[u][:]
v_m_vs_energy = v_m[u][:]
dv_vs_energy = dv[u][:]
dos_p_vs_energy = dos_p[u][:]
dos_m_vs_energy = dos_m[u][:]
ldos_p_vs_energy = ldos_p[u][:]
ldos_m_vs_energy = ldos_m[u][:]
rdos_p_vs_energy = rdos_p[u][:]
rdos_m_vs_energy = rdos_m[u][:]
udos_p_vs_energy = udos_p[u][:]
udos_m_vs_energy = udos_m[u][:]


# Make a 3x1 subplot of conductances.
fig3 = plt.figure()
fig3.suptitle('Conductances (' + spin + ', T = ' + Temperature + 't, upos = ' + str(upos[u]) + 'a)', fontsize=16, fontweight='bold')
sub1 = fig3.add_subplot(311)
sub1.grid(True)
sub1.set_ylabel('$G_{UL} \ [e^2/h]$')
sub1.set_ylim((0, 1))
sub1.plot(energies, g_ul_p_vs_energy, 'b-', label='p')
sub1.plot(energies, g_ul_m_vs_energy, 'r-', label='m')
# plt.legend(bbox_to_anchor=(0., 1.16, 1., .102), loc=3,ncol=3,
#            mode="expand", borderaxespad=0.)
plt.legend(shadow=True, fancybox=True)

sub2 = fig3.add_subplot(312)
sub2.grid(True)
sub2.set_ylabel('$G_{UR} \ [e^2/h]$')
sub2.set_ylim((0, 3))
sub2.plot(energies, g_ur_p_vs_energy, 'b-')
sub2.plot(energies, g_ur_m_vs_energy, 'r-')

sub3 = fig3.add_subplot(313)
sub3.grid(True)
sub3.set_ylabel('$G_{RL} \ [e^2/h]$')
sub3.set_xlabel('Energy [t]')
sub3.set_ylim((0, 10))
sub3.plot(energies, g_rl_p_vs_energy, 'b-')
sub3.plot(energies, g_rl_m_vs_energy, 'r-')

# Fine-tune figure; make subplots close to each other
# and hide x ticks for all but bottom plot.
plt.setp([a.get_xticklabels() for a in fig3.axes[:-1]], visible=False)
plt.subplots_adjust(left=0.11, bottom=0.09, right=0.96, top=0.92,
                    wspace=None, hspace=0.15)
plt.show()

# Make a 3x1 subplot with voltages, dV, and full system DoS.
fig4 = plt.figure()
fig4.suptitle('Voltages (' + spin + ', T = ' + Temperature + 't, upos = ' + str(upos[u]) + 'a)', fontsize=16, fontweight='bold')
sub1 = fig4.add_subplot(311)
sub1.grid(True)
sub1.set_ylabel('$V_U/V_L$')
sub1.set_ylim((0, 1))
sub1.plot(energies, v_p_vs_energy, 'b-', label='p')
sub1.plot(energies, v_m_vs_energy, 'r-', label='m')
# plt.legend(bbox_to_anchor=(0., 1.1, 1., .102), loc=3,
#            ncol=3, mode="expand", borderaxespad=0.)
plt.legend(shadow=True, fancybox=True)

sub2 = fig4.add_subplot(312)
sub2.grid(True)
sub2.set_ylabel('Difference in Voltages (dV)')
sub2.set_xlabel('Energy [t]')
sub2.set_ylim((-0.7, 0.7))
sub2.plot(energies, dv_vs_energy, 'k-', label='dV')
plt.legend(shadow=True, fancybox=True)

sub3 = fig4.add_subplot(313)
sub3.grid(True)
sub3.set_ylabel('DoS [t^-1 a^-2]')
sub3.set_xlabel('Energy [t]')
sub3.set_ylim((0, 1))
sub3.plot(energies, dos_p_vs_energy, 'b-', label='p')
sub3.plot(energies, dos_m_vs_energy, 'r-', label='m')
plt.legend(shadow=True, fancybox=True)


# Fine-tune figure; make subplots close to each other
# and hide x ticks for all but bottom plot.
plt.setp([a.get_xticklabels() for a in fig4.axes[:-1]], visible=False)
plt.subplots_adjust(left=0.11, bottom=0.09, right=0.96, top=0.92,
                    wspace=None, hspace=0.15)
plt.show()


# Make a 3x1 subplot with dos in the leads.
fig5 = plt.figure()
fig5.suptitle('DoS (' + spin + ', T = ' + Temperature + 't, upos = ' 
              + str(upos[u]) + 'a)', fontsize=16, fontweight='bold')
sub1 = fig5.add_subplot(311)
sub1.grid(True)
sub1.set_ylabel('Left DoS [t^-1 a^-2]')
sub1.set_ylim((0, 1.5))
sub1.plot(energies, ldos_p_vs_energy, 'b-', label='p')
sub1.plot(energies, ldos_m_vs_energy, 'r-', label='m')
# plt.legend(bbox_to_anchor=(0., 1.1, 1., .102), loc=3,
#            ncol=3, mode="expand", borderaxespad=0.)
plt.legend(shadow=True, fancybox=True)

sub2 = fig5.add_subplot(312)
sub2.grid(True)
sub2.set_ylabel('Right DoS [t^-1 a^-2]')
sub2.set_xlabel('Energy [t]')
sub2.set_ylim((0, 1.5))
sub2.plot(energies, rdos_p_vs_energy, 'b-', label='p')
sub2.plot(energies, rdos_m_vs_energy, 'r-', label='m')
plt.legend(shadow=True, fancybox=True)

sub3 = fig5.add_subplot(313)
sub3.grid(True)
sub3.set_ylabel('Upper DoS [t^-1 a^-2]')
sub3.set_xlabel('Energy [t]')
sub3.set_ylim((0, 1.5))
sub3.plot(energies, udos_p_vs_energy, 'b-', label='p')
sub3.plot(energies, udos_m_vs_energy, 'r-', label='m')
plt.legend(shadow=True, fancybox=True)


# Fine-tune figure; make subplots close to each other
# and hide x ticks for all but bottom plot.
plt.setp([a.get_xticklabels() for a in fig5.axes[:-1]], visible=False)
plt.subplots_adjust(left=0.11, bottom=0.09, right=0.96, top=0.92,
                    wspace=None, hspace=0.15)
plt.show()




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
plt.show()

# Make a 2x1 subplot with voltages and dV.
fig2 = plt.figure()
fig2.suptitle('Voltages (' + spin + ', $\mu$ = ' + str(energies[index]) + 't, T = ' + Temperature + 't)', fontsize=16, fontweight='bold')
sub1 = fig2.add_subplot(211)
sub1.grid(True)
sub1.set_ylabel('$V_U/V_L$')
sub1.set_ylim((0, 1))
sub1.plot(upos, v_p_vs_upos, 'b-', label='p')
sub1.plot(upos, v_m_vs_upos, 'r-', label='m')
# plt.legend(bbox_to_anchor=(0., 1.1, 1., .102), loc=3,
#            ncol=3, mode="expand", borderaxespad=0.)
plt.legend(shadow=True, fancybox=True)

sub2 = fig2.add_subplot(212)
sub2.grid(True)
sub2.set_ylabel('Difference in Voltages (dV)')
sub2.set_xlabel('Upper Lead Position [a]')
sub2.set_ylim((-0.7, 0.7))
sub2.plot(upos, dv_vs_upos, 'k-', label='dV')
plt.legend(shadow=True, fancybox=True)

# Fine-tune figure; make subplots close to each other and hide x ticks for
# all but bottom plot.
plt.setp([a.get_xticklabels() for a in fig2.axes[:-1]], visible=False)
plt.subplots_adjust(left=0.11, bottom=0.09, right=0.96, top=0.92,
                    wspace=None, hspace=0.15)
plt.show()

print "Saving plots to " + filename + ".pdf"


pp = PdfPages(filename + '.pdf')
pp.savefig(fig1)
pp.savefig(fig2)
pp.savefig(fig3)
pp.savefig(fig4)
pp.savefig(fig5)
pp.close()


