#!/usr/bin/env python

'''so that 1/2 == 0.5, and not 0'''
from __future__ import division
import sys
import datetime
import time
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import tinyarray
import numpy as np
from math import sin, cos, cosh, pi, isnan
import operator as op

import kwant

# define Pauli-matrices for convenience
s_0 = tinyarray.array([[1, 0], [0, 1]])
s_x = tinyarray.array([[0, 1], [1, 0]])
s_y = tinyarray.array([[0, -1j], [1j, 0]])
s_z = tinyarray.array([[1, 0], [0, -1]])


def divide_lists(list1, list2):
    """Combine two lists according to the specified operation (
    add, sub, div).
    """
    new_list = []
    if len(list1) != len(list2):
        print("Lengths of lists are not the same.")
        return 0
    else:
        for i in range(len(list1)):
            if list1[i] == 0 or list2[i] == 0:
                new_list.append(0)
            else:
                new_list.append(list1[i]/list2[i])
        return new_list


def floating_voltage(g_ul, g_ur):
    """Compute the voltage required to make U float while L is
    ferromagnetic.
    """
    volt = map(op.truediv, g_ul, map(op.add, g_ul, g_ur))
    for i in range(len(volt)):
        if isnan(volt[i]) is True:
            volt[i] = 0
    return volt


def bandstructure(flead, momenta, show=False):
    """Compute the band structure for each transverse
    subband.
    """
    bands = kwant.physics.Bands(flead)
    energies = [bands(k) for k in momenta]
    if show is True:
        plt.figure()
        plt.title("Band Structure")
        plt.plot(momenta, energies)
        plt.xlabel("Momentum [a^-1]")
        plt.ylabel("Energy [t]")
        plt.grid(True)
        plt.show()
    return energies


def band_minima(band):
    band = np.array(band).transpose()
    minima = []
    for i in range(len(band)):
        x = round(min(band[i]), 10)
        if x not in minima:
            minima.append(x)
            minima.append(-x)
    return minima


def conductance(smatrix, energies, to_lead, from_lead, show=False):
    """Calculate and plot conductance. Return conductance."""
    cond = []
    for i in range(len(smatrix)):
        cond.append(smatrix[i].transmission(to_lead, from_lead))
    if show is True:
        plt.figure()
        plt.plot(energies, cond)
        plt.xlabel("Energy [t]")
        plt.ylabel("Conductance [e^2/h]")
        plt.title("Conductance from Left to Right")
        plt.grid(True)
        plt.show()
    return cond


def density_of_states(sys, energies, L, W, show=False):
    """Calculate the density of states in system"""
    dos = []
    for energy in energies:
        dos.append(np.sum(kwant.ldos(sys, energy)))
    dos = [x / (W*L) for x in dos]
    if show is True:
        plt.figure()
        plt.title("Density of States")
        plt.plot(energies, dos)
        plt.xlabel("Energies [t]")
        plt.ylabel("Density of States [t^-1 a^-2]")
        plt.grid(True)
        plt.show()
    return dos


def fermi_scale(bands, momenta, energies, stop=False):
    ''' This finds all of the Fermi momenta that contribute to transport
        at a given energy (chemical potential) and stores them in a list.
        Then, there is a list of these lists corresponding to each
        different energy.
        This function returns a list of lists. 

        Arguments:
        bands -- list (of numpy arrays) returned by bandstructure()
        momenta -- list used by bandstructure()
        energies -- list of energies to calculate fermi momenta

        Returns: list (of lists)
    '''
    bands = np.array(bands)
    bands = bands.transpose()

    fermi_momentum = []
    fermi_wavelength = []
    for energy in energies:
        fermi_momentum_E = [[energy]]
        fermi_wavelength_E = [[energy]]
        for band_number in range(bands.shape[0]):
            band = bands[band_number]
            band_prev = band[0]
            final_count = float('nan')
            for i in range(1, int(len(momenta)/2)):
                count = i
                band_next = band[i]
                if (band_prev <= energy <= band_next) or (band_prev >= energy >= band_next):
                    final_count = count
                    fermo = 0.5*(momenta[count] + momenta[count - 1])
                    fermwl = 2 * pi / fermo
                if isnan(final_count):
                    fermo = float('nan')
                    fermwl = float('nan')
                band_prev = band_next
            fermi_momentum_E.append(fermo)
            fermi_wavelength_E.append(fermwl)
        fermi_momentum.append(fermi_momentum_E)
        fermi_wavelength.append(fermi_wavelength_E)
    return fermi_momentum, fermi_wavelength
    if stop is True:
        exit()


def local_dos(sys, energies, L, W, show=False):
    """Calculate the local density of states of the system. Returns
    the total ldos of the system, spin up ldos, and spin down
    ldos as numpy arrays."""
    ldos = np.empty((0, (L+1)*(W+1)))
    ldos_sys_p = np.empty((0, (L+1)*(W+1)))
    ldos_sys_m = np.empty((0, (L+1)*(W+1)))
    for energy in energies:
        ldos_e = kwant.ldos(sys, energy)
        # Extract the individual spin up and spin down ldos components.
        ldos_sys_p_e = ldos_e.reshape(-1, 2)[:, 0]
        ldos_sys_m_e = ldos_e.reshape(-1, 2)[:, 1]
        # Calculate ldos per site, by summing spin up and
        # spin down components.
        ldos_e = np.sum(ldos_e.reshape(-1, 2), axis=1)
        ldos = np.append(ldos, [ldos_e], axis=0)
        ldos_sys_p = np.append(ldos_sys_p, [ldos_sys_p_e], axis=0)
        ldos_sys_m = np.append(ldos_sys_m, [ldos_sys_m_e], axis=0)
    if show is True:
        kwant.plotter.map(sys, ldos_sys_p, num_lead_cells=3)
        kwant.plotter.map(sys, ldos_sys_m, num_lead_cells=3)
        kwant.plotter.map(sys, ldos, num_lead_cells=3)
    return ldos, ldos_sys_p, ldos_sys_m


def dos_leads(energies, band, momenta, W, show=False):
    band = np.array(band).transpose()
    dos = []
    for chempot in energies:
        __ = 0
        for i in range(len(band)):
            for j in range(len(momenta)-1):
                if (band[i][j] < chempot < band[i][j+1]) or (band[i][j+1] < chempot < band[i][j]):
                    __ = __ + 1/(abs(band[i][j] - band[i][j+1])*W*len(energies))
        dos.append(__)
    if show is True:
        plt.plot(energies, dos)
        plt.xlabel("Energy [t]")
        plt.ylabel("Dos [t^-1 a^-2]")
        plt.title("Density of States of Lead")
        plt.show()
    return dos


def s_matrix(sys, energies):
    """Compute the S-Matrix as a function of energy.
    Return a list of s of type
    <class 'kwant.solvers.common.SMatrix'>.
    """
    s = []
    for energy in energies:
        s.append(kwant.smatrix(sys, energy))
    return s


def temp_smooth(par, energies, dE, N_steps, temperature):
    if temperature == 0:
        return par
    else:
        cond_T = []
        for j in range(len(energies)):
            el = 0
            for i in range(len(energies)):
                arg = (energies[i]-energies[j])/temperature/2
                if abs(arg) > 300:
                    correction_factor = 0
                else:
                    correction_factor = 1/cosh((energies[i]-energies[j])
                                               / temperature/2)**2
                el = el + 0.25*par[i]*(dE/temperature)*correction_factor
            cond_T.append(el)
        return cond_T


def write_datafile(filename, L, W, srashba, barrier, lshift, lferro, theta,
                   phi, uwidth, ushift, T, N_steps, dE, energies, upos,
                   elapsed_time, ar_shape, lfermi_momentum, lfermi_wavelength,
                   rfermi_momentum, rfermi_wavelength, ufermi_momentum,
                   ufermi_wavelength):
    filename = filename + '.par'
    print 'Writing parameters to ' + filename
    with open(filename, 'w') as f:
        f.write(str(datetime.datetime.today()) + '\n \n')
        f.write('System (scattering region)\n')
        f.write('L = ' + str(L) + '\n')
        f.write('W = ' + str(W) + '\n')
        f.write('srashba = ' + str(srashba) + '\n \n')
        f.write('Left Lead\n')
        f.write('lshift = ' + str(lshift) + '\n')
        f.write('lferro = ' + str(lferro) + '\n')
        f.write('theta = ' + str(theta) + '\n')
        f.write('phi = ' + str(phi) + '\n \n')
        f.write('Upper Lead\n')
        f.write('uwidth = ' + str(uwidth) + '\n \n')
        f.write('Misc\n')
        f.write('T = ' + str(T) + '\n')
        f.write('N_steps = ' + str(N_steps) + '\n')
        f.write('dE = ' + str(dE) + '\n')
        f.write('energies = ' + str(energies) + '\n \n')
        f.write('upos = ' + str(upos) + '\n \n')
        f.write('Elapsed Time = ' + str(round(elapsed_time, 1))
                + ' seconds' + '\n \n')
        f.write('Array shape = ' + str(ar_shape) + '\n \n')
        f.write('Barrier = ' + str(barrier) + '\n')
        f.write('ushift = ' + str(ushift) + '\n')
        f.write('left fermi momenta = ' + str(lfermi_momentum) + '\n')
        f.write('left fermi wavelength = ' + str(lfermi_wavelength) + '\n')
        f.write('right fermi momenta = ' + str(rfermi_momentum) + '\n')
        f.write('right fermi wavelength = ' + str(rfermi_wavelength) + '\n')
        f.write('upper fermi momenta = ' + str(ufermi_momentum) + '\n')
        f.write('upper fermi wavelength = ' + str(ufermi_wavelength) + '\n')
    f.close()


def write_output(filename, upos, energies, g_ul_p, g_ul_m, g_ur_p,
                 g_ur_m, g_rl_p, g_rl_m, v_p, v_m, dv, dos_p, dos_m,
                 ldos_p, ldos_m, rdos_p, rdos_m, udos_p, udos_m):
    print 'Saving output to ' + filename + '.data'
    filename = filename + '.data'
    data = np.stack((np.array(g_ul_p), np.array(g_ul_m), np.array(g_ur_p),
                     np.array(g_ur_m), np.array(g_rl_p), np.array(g_rl_m),
                     np.array(v_p), np.array(v_m), np.array(dv),
                     np.array(dos_p), np.array(dos_m), np.array(ldos_p),
                     np.array(ldos_m), np.array(rdos_p), np.array(rdos_m),
                     np.array(udos_p), np.array(udos_m)), axis=0)
    with open(filename, 'w') as f:
        # Write shape of numpy array
        # Any line starting with '#' will be ignored by numpy.loadtxt
        f.write('# Array shape: {0}\n'.format(data.shape))
        # Write values, ulead position, and energies
        data_names = ['g_ul_p', 'g_ul_m', 'g_ur_p', 'g_ur_m', 'g_rl_p',
                      'g_rl_m', 'v_p', 'v_m', 'dv', 'dos_p','dos_m',
                      'ldos_p', 'ldos_m', 'rdos_p', 'rdos_m', 'udos_p',
                      'udos_m']
        f.write('# data = [g_ul_p, g_ul_m, g_ur_p, g_ur_m, '
                + 'g_rl_p, g_rl_m, v_p, v_m, dv, dos_p, dos_m, '
                + 'ldos_p, ldos_m, rdos_p, rdos_m, udos_p, udos_m]'
                + '\n')
        f.write('# pos = ' + str(upos) + '\n')
        f.write('# energies = ' + str(energies) + '\n')
        # Iterating through an ndimensional array produces slices along
        # the last axis. This is equivalent to data[i,:,:] in this case
        i = 0
        for data_slice in data:
            # Writing out a break to indicate different slices...
            f.write('# ' + data_names[i] + '\n')

            # The formatting string indicates that I'm writing out
            # the values in left-justified columns 7 characters in width
            # with 2 decimal places.
            np.savetxt(f, data_slice, fmt='%-2.8f')
            i = i + 1
    return data.shape


def write_dos(filename, sys, energies, ldos_sys, ldos_sys_p, ldos_sys_m):
    print 'Saving DoS to ' + filename + '.dos'
    filename = filename + '.dos'
    data = np.stack((ldos_sys, ldos_sys_p, ldos_sys_m), axis=0)

    with open(filename, 'w') as f:
        # Write shape of numpy array
        # Any line starting with '#' will be ignored by numpy.loadtxt
        f.write('# Array shape: {0}\n'.format(data.shape))
        # Write headers
        headers = ['ldos_sys', 'ldos_sys_p', 'ldos_sys_m']
        f.write('# columns = [ldos_sys, ldos_sys_p, ldos_sys_m]' + '\n')
        f.write('# energies = ' + str(energies) + '\n')
        # Iterating through an ndimensional array produces slices along
        # the last axis. This is equivalent to data[i,:,:] in this case
        i = 0
        for data_slice in data:
            # Writing out a break to indicate different slices...
            f.write('# ' + headers[i] + '\n')

            # The formatting string indicates that I'm writing out
            # the values in left-justified columns 7 characters in width
            # with 2 decimal places.
            np.savetxt(f, data_slice, fmt='%-2.8f')
            i = i + 1
    return data.shape

def write_plots(filename, index, energies, upos, g_ul_p, g_ul_m, g_ur_p,
                g_ur_m, g_rl_p, g_rl_m, v_p, v_m, dv, theta, phi):
    print 'Saving plots to ' + filename + '.pdf'
    g_ul_p_vs_upos = [g_ul_p[i][index] for i in range(len(g_ul_p))]
    g_ul_m_vs_upos = [g_ul_m[i][index] for i in range(len(g_ul_m))]
    g_ur_p_vs_upos = [g_ur_p[i][index] for i in range(len(g_ur_p))]
    g_ur_m_vs_upos = [g_ur_m[i][index] for i in range(len(g_ur_m))]
    g_rl_p_vs_upos = [g_rl_p[i][index] for i in range(len(g_rl_p))]
    g_rl_m_vs_upos = [g_rl_m[i][index] for i in range(len(g_rl_m))]
    v_p_vs_upos = [v_p[i][index] for i in range(len(v_p))]
    v_m_vs_upos = [v_m[i][index] for i in range(len(v_m))]
    dv_vs_upos = [dv[i][index] for i in range(len(dv))]

    # Get the FM orientation (initial spin polarization)
    # to use in title.
    if theta == 0 and phi == 0:
        spin = 'Sz'
    elif theta != 0 and phi == 0:
        spin = 'Sx'
    elif theta != 0 and phi != 0:
        spin = 'Sy'
    cp_str = str(energies[index])

    # Make a 3x1 subplot of conductances.
    fig1 = plt.figure()
    fig1.suptitle('Conductances (' + spin + ', $\mu$ = ' + cp_str + 't)',
                  fontsize=16, fontweight='bold')
    sub1 = fig1.add_subplot(311)
    sub1.grid(True)
    sub1.set_ylabel('$G_{UL} \ [e^2/h]$')
    sub1.set_ylim((0, 1.0))
    sub1.plot(upos, g_ul_p_vs_upos, 'b-', label='p')
    sub1.plot(upos, g_ul_m_vs_upos, 'r-', label='m')
    # plt.legend(bbox_to_anchor=(0., 1.16, 1., .102), loc=3,ncol=3,
    #            mode="expand", borderaxespad=0.)
    plt.legend(shadow=True, fancybox=True)

    sub2 = fig1.add_subplot(312)
    sub2.grid(True)
    sub2.set_ylabel('$G_{UR} \ [e^2/h]$')
    sub2.set_ylim((0, 1.5))
    sub2.plot(upos, g_ur_p_vs_upos, 'b-')
    sub2.plot(upos, g_ur_m_vs_upos, 'r-')

    sub3 = fig1.add_subplot(313)
    sub3.grid(True)
    sub3.set_ylabel('$G_{RL} \ [e^2/h]$')
    sub3.set_xlabel('Upper Lead Position [a]')
    sub3.set_ylim((0, 8))
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
    fig2.suptitle('Voltages (' + spin + ', $\mu$ = ' + cp_str + 't)',
                  fontsize=16, fontweight='bold')
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

    # Fine-tune figure; make subplots close to each other
    # and hide x ticks for all but bottom plot.
    plt.setp([a.get_xticklabels() for a in fig2.axes[:-1]], visible=False)
    plt.subplots_adjust(left=0.11, bottom=0.09, right=0.96, top=0.92,
                        wspace=None, hspace=0.15)
    plt.show()

    pp = PdfPages(filename + '.pdf')
    pp.savefig(fig1)
    pp.savefig(fig2)
    pp.close()


def make_system(a=1, t=1, L=10, W=10, srashba=0, barrier=0, lshift=0, lferro=0,
                theta=0, phi=0, uleftedge=5, uwidth=2, ushift=0, show=False):
    """Create a tight-binding system on a square lattice.

    Keyword arguments:
    a -- Lattice constant (default 1)
    t -- Hopping amplitude (default 1)
    L -- Length of scattering region (default 10)
    W -- Width of scattering region (default 10)
    srashba -- Rashba strength in scattering region (default 0)
    lshift -- Net energy shift in left lead (default 0)
    lferro -- FM strength in left lead (default 0)
    theta -- FM orientation (default 0)
    phi -- FM orientation (default 0)
    uleftedge -- Starting position of upper lead (default 5)
    uwidth -- Width of upper lead (default 2)
    """
    square = kwant.lattice.square(a)

    # Define the scattering region shape.
    def rectangle(pos):
        x, y = pos
        return (0 <= x <= L) and (0 <= y <= W)

    sys = kwant.Builder()

    # Define the scattering region.
    sys[square.shape(rectangle, (0, 0))] = 4 * t * s_0 * 0
    for j in range(0, W+1):
        sys[square(0, j)] = (4 * t * 0 + barrier) * s_0
    # Hoppings in the x-direction.
    sys[kwant.builder.HoppingKind((1, 0), square, square)] = (-t * s_0
                                                              - 1j * srashba
                                                              * s_y)
    # Hoppings in the y-direction.
    sys[kwant.builder.HoppingKind((0, 1), square, square)] = (-t * s_0
                                                              + 1j * srashba
                                                              * s_x)

    # Plot closed system.
    # kwant.plot(sys)

    # Define the left and right lead shapes.
    def lead_shape(pos):
        x, y = pos
        return (0 <= y <= W)

    # Define the left lead.
    llead = kwant.Builder(kwant.TranslationalSymmetry((-1, 0)))
    llead[square.shape(lead_shape, (0, 0))] = (4 * t * s_0 * 0
                                               + (lshift) * s_0
                                               - lferro
                                               * (sin(theta)*cos(phi)*s_x
                                                  + sin(theta)*sin(phi)*s_y
                                                  + cos(theta)*s_z))
    llead[kwant.builder.HoppingKind((1, 0), square, square)] = -t * s_0
    llead[kwant.builder.HoppingKind((0, 1), square, square)] = -t * s_0

    # Define the right lead.
    rlead = kwant.Builder(kwant.TranslationalSymmetry((1, 0)))
    rlead[square.shape(lead_shape, (0, 0))] = (4 * t * 0) * s_0
    rlead[kwant.builder.HoppingKind((1, 0), square, square)] = -t * s_0
    rlead[kwant.builder.HoppingKind((0, 1), square, square)] = -t * s_0

    # Define the upper lead.
    ulead = kwant.Builder(kwant.TranslationalSymmetry((0, 1)))
    ulead[(square(i, 0) for i in
          xrange(uleftedge, uleftedge + uwidth))] = (4 * t * 0 + ushift) * s_0
    ulead[square.neighbors()] = -t * s_0

    # Attach leads to the system
    sys.attach_lead(llead)
    sys.attach_lead(rlead)
    sys.attach_lead(ulead)
    if show is True:
        kwant.plot(sys)

    return sys, [llead, rlead, ulead]


def main():
    start_time = time.time()

    # Open a file to save parameters
    if len(sys.argv) != 2:
        filename = 'temp'
    else:
        filename = str(sys.argv[1])
    print 'Writing output to ' + filename + ".data"

    # Define parameters to make system
    L = 400
    W = 10
    srashba = 0.0
    barrier = 0.0
    lshift = 4
    lferro = 4
    theta = pi/2
    phi = 0
    uwidth = 2
    ushift = 0.0

    T = 0.0
    N_steps = 81
    dE = 0.1
    energies = [round(dE*i - 4, 1) for i in range(N_steps)]
    momenta = [-pi + 0.02*pi*i for i in range(100)]

    upos = []
    g_ul_p = []
    g_ul_m = []
    g_ur_p = []
    g_ur_m = []
    g_rl_p = []
    g_rl_m = []
    v_p = []
    v_m = []
    dv = []
    dos_p = []
    dos_m = []

    for pos in range(0, L, 1):
        upos.append(pos)
        print pos
        sys_p, leads_p = make_system(L=L, W=W, srashba=srashba, barrier=barrier,
                                     lshift=lshift, lferro=lferro,
                                     theta=theta, phi=phi,
                                     uleftedge=pos, uwidth=uwidth, ushift=ushift)
        sys_m, leads_m = make_system(L=L, W=W, srashba=srashba, barrier=barrier,
                                     lshift=lshift, lferro=-lferro,
                                     theta=theta, phi=phi,
                                     uleftedge=pos, uwidth=uwidth, ushift=ushift)
        # Finalize the system and leads to use in calculations.
        sys_p = sys_p.finalized()
        sys_m = sys_m.finalized()
        [llead_p, rlead_p, ulead_p] = [leads_p[i].finalized() for i
                                       in range(len(leads_p))]
        [llead_m, rlead_m, ulead_m] = [leads_m[i].finalized() for i
                                       in range(len(leads_m))]
        if pos == 0:
            # Calculate the band structure for each lead.
            lband_p = bandstructure(llead_p, momenta)
            lband_m = bandstructure(llead_m, momenta)
            rband_p = bandstructure(rlead_p, momenta)
            rband_m = bandstructure(rlead_m, momenta)
            uband_p = bandstructure(ulead_p, momenta)
            uband_m = bandstructure(ulead_m, momenta)
            # Calculate the minima of the transverse bands.
            lbandmin = band_minima(lband_p)
            rbandmin = band_minima(rband_p)
            ubandmin = band_minima(uband_p)
            # Don't include lbandmin, since it is
            # duplicated by rbandmin.
            # print lbandmin, rbandmin, ubandmin
            energies = energies + rbandmin + ubandmin
            energies.sort()
            print energies
            lfermi_momentum, lfermi_wavelength = fermi_scale(lband_p, momenta, energies, stop=False)
            rfermi_momentum, rfermi_wavelength = fermi_scale(rband_p, momenta, energies, stop=False)
            ufermi_momentum, ufermi_wavelength = fermi_scale(uband_p, momenta, energies, stop=False)

            # Calculate local density of states of the system
            ldos_sys, ldos_sys_p, ldos_sys_m = local_dos(sys_p, energies, L, W)

        s_p = s_matrix(sys_p, energies)
        s_m = s_matrix(sys_m, energies)
        # Zero temperature conductances and voltages.
        g_ul_p_0 = conductance(s_p, energies, 2, 0)
        g_ul_m_0 = conductance(s_m, energies, 2, 0)
        g_ur_p_0 = conductance(s_p, energies, 2, 1)
        g_ur_m_0 = conductance(s_m, energies, 2, 1)
        g_rl_p_0 = conductance(s_p, energies, 1, 0)
        g_rl_m_0 = conductance(s_m, energies, 1, 0)
        v_p_0 = floating_voltage(g_ul_p_0, g_ur_p_0)
        v_m_0 = floating_voltage(g_ul_m_0, g_ur_m_0)
        dos_p_0 = density_of_states(sys_p, energies, L, W)
        dos_m_0 = density_of_states(sys_m, energies, L, W)
        # Save the temperature-smoothed conductances and voltages.
        g_ul_p.append(temp_smooth(g_ul_p_0, energies, dE, N_steps, T))
        g_ul_m.append(temp_smooth(g_ul_m_0, energies, dE, N_steps, T))
        g_ur_p.append(temp_smooth(g_ur_p_0, energies, dE, N_steps, T))
        g_ur_m.append(temp_smooth(g_ur_m_0, energies, dE, N_steps, T))
        g_rl_p.append(temp_smooth(g_rl_p_0, energies, dE, N_steps, T))
        g_rl_m.append(temp_smooth(g_rl_m_0, energies, dE, N_steps, T))
        v_p_T = temp_smooth(v_p_0, energies, dE, N_steps, T)
        v_m_T = temp_smooth(v_m_0, energies, dE, N_steps, T)
        v_p.append(v_p_T)
        v_m.append(v_m_T)
        # Difference of temperature-smoothed voltage.
        dv.append(map(op.sub, v_p_T, v_m_T))
        dos_p.append(temp_smooth(dos_p_0,energies, dE, N_steps, T))
        dos_m.append(temp_smooth(dos_m_0,energies, dE, N_steps, T))

    # Calculate DoS in each lead.
    ldos_p = dos_leads(energies, lband_p, momenta, W)
    ldos_m = dos_leads(energies, lband_m, momenta, W)
    rdos_p = dos_leads(energies, rband_p, momenta, W)
    rdos_m = dos_leads(energies, rband_m, momenta, W)
    udos_p = dos_leads(energies, uband_p, momenta, uwidth)
    udos_m = dos_leads(energies, uband_m, momenta, uwidth)
    # Store dos in the same format as all other data.
    ldos_p = [ldos_p for i in range(len(upos))]
    ldos_m = [ldos_m for i in range(len(upos))]
    rdos_p = [rdos_p for i in range(len(upos))]
    rdos_m = [rdos_m for i in range(len(upos))]
    udos_p = [udos_p for i in range(len(upos))]
    udos_m = [udos_m for i in range(len(upos))]

    print '=============================='
    elapsed_time = time.time() - start_time
    print 'Elapsed time: ' + str(round(elapsed_time, 1)) + ' seconds'
    print '=============================='

    # Write out results to .data file
    ar_shape = write_output(filename, upos, energies, g_ul_p, g_ul_m,
                            g_ur_p, g_ur_m, g_rl_p, g_rl_m, v_p, v_m,
                            dv, dos_p, dos_m, ldos_p, ldos_m, rdos_p,
                            rdos_m, udos_p, udos_m)

    # Write out parameters to .par file
    write_datafile(filename, L, W, srashba, barrier, lshift, lferro, theta, phi,
                   uwidth, ushift, T, N_steps, dE, energies, upos,
                   elapsed_time, ar_shape, lfermi_momentum, lfermi_wavelength,
                   rfermi_momentum, rfermi_wavelength, ufermi_momentum,
                   ufermi_wavelength)

    # Write dos to .dos file.
    write_dos(filename, sys_p, energies, ldos_sys, ldos_sys_p, ldos_sys_m)

# Call the main function if the script gets executed (as
# opposed to imported).
# See <http://docs.python.org/library/__main__.html>.
if __name__ == '__main__':
    main()
