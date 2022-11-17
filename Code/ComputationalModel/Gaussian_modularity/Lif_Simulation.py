"""
Attenzione, non posso fare toggle tra input
"""
from ExternalFunctions import crea_e_o_apri, savelookuptable, savepeaktrain
import os
from datetime import date
from Lif_Initialization import *
from Lif_Plot import *

import threading

seed(val)
###########################################
#       User defined simulation parameters
###########################################
#  ---------- Simulation ----------
preptime = 10 * second  # Length of not recorded simulation to settle the network (preparation time) - minimum=1*second
simtime = 60 * second  # Length of the simulation
# zoom = 10 * second          # portion of signal you want to plot
chunk = 10  # Length of simulation chunks
should_save = True  # False=do not save output  -   True=save output
should_show = False  # False=do not show figures
DateFolder = str(date.today())  # Date in the format aaaa-mm-dd
SimulationNumber = 3  # Number of simulation in the date

comment = '\nNumber of layers:\t' + str(num_layers) + \
          '\nUpward connections:\t' + str(upward_conn) + \
          '\nDownward connections:\t' + str(downward_conn) + \
          '\nspatial config--> index changes 4 each layer' \
          '\nrefractory period = ' + str(refractory_period) + \
          '\nperc_conn_ext_exc = ' + str(perc_conn_ext_exc) + \
          '\nperc_conn_ext_inh = ' + str(perc_conn_ext_inh)

defaultclock.dt = 0.1 * ms  # Clock of the simulation

#  ---------- Input ----------
inputtype = 1  # Input Type: 0=step input; 1=noisy input
DURstim = 100 * ms  # Duration of the step input
Tstim = 500 * ms  # Initial time of the step input
input_value = '(randint(0, 5, size=n_cells) * 6.4)*mV'


###########################################


#####################################
#         Input definition
#####################################
@network_operation(dt=1 * ms)
def testinput(t):
    if inputtype == 1:
        # ------ Noisy input ------
        for pop in layers:
            pop.I_noise = (randint(0, 5, size=n_cells) * 6.4) * mV  #### Fluttuazioni Vm: modificabili


#####################################

#####################################
#         Output folder definition
#####################################
if should_save:
    crea_e_o_apri('Results_3D_LIF')
    crea_e_o_apri(DateFolder)
    crea_e_o_apri('Simulazione' + str(SimulationNumber))
    SimDir = os.getcwd()
    ImagesDir = os.path.join(SimDir, 'Images/')
    if not os.path.isdir(ImagesDir):
        os.makedirs(ImagesDir)
    LutDir = os.path.join(SimDir, 'LookupTable/')
    if not os.path.isdir(LutDir):
        os.makedirs(LutDir)
    NumConnDir = os.path.join(SimDir, 'NumConnections/')
    if not os.path.isdir(NumConnDir):
        os.makedirs(NumConnDir)
    PeakTrainsDir = os.path.join(SimDir, 'PeakTrains/')
    if not os.path.isdir(PeakTrainsDir):
        os.makedirs(PeakTrainsDir)
    WeightDir = os.path.join(SimDir, 'WeightMatrix/')
    if not os.path.isdir(WeightDir):
        os.makedirs(WeightDir)

#####################################


#####################################
#         Simulation & Online Saving
#####################################
net = Network(collect())
net.run(preptime, report='text')

# Sx
sm_layer0_sx = SpikeMonitor(l0_sx, variables='vm')
net.add(sm_layer0_sx)
if num_layers > 1:
    sm_layer1_sx = SpikeMonitor(l1_sx, variables='vm')
    net.add(sm_layer1_sx)
if num_layers > 2:
    sm_layer2_sx = SpikeMonitor(l2_sx, variables='vm')
    net.add(sm_layer2_sx)
if num_layers > 3:
    sm_layer3_sx = SpikeMonitor(l3_sx, variables='vm')
    net.add(sm_layer3_sx)
if num_layers > 4:
    sm_layer4_sx = SpikeMonitor(l4_sx, variables='vm')
    net.add(sm_layer4_sx)
if num_layers > 5:
    sm_layer5_sx = SpikeMonitor(l5_sx, variables='vm')
    net.add(sm_layer5_sx)

# Dx
sm_layer0_dx = SpikeMonitor(l0_dx, variables='vm')
net.add(sm_layer0_dx)
if num_layers > 1:
    sm_layer1_dx = SpikeMonitor(l1_dx, variables='vm')
    net.add(sm_layer1_dx)
if num_layers > 2:
    sm_layer2_dx = SpikeMonitor(l2_dx, variables='vm')
    net.add(sm_layer2_dx)
if num_layers > 3:
    sm_layer3_dx = SpikeMonitor(l3_dx, variables='vm')
    net.add(sm_layer3_dx)
if num_layers > 4:
    sm_layer4_dx = SpikeMonitor(l4_dx, variables='vm')
    net.add(sm_layer4_dx)
if num_layers > 5:
    sm_layer5_dx = SpikeMonitor(l5_dx, variables='vm')
    net.add(sm_layer5_dx)

if should_save:
    quotient = simtime / second // chunk
    remainder = simtime / second % chunk
    times = np.repeat(chunk * second, quotient).tolist()
    times.append(remainder * second)

    os.chdir(PeakTrainsDir)
    for _, tt in enumerate(times):
        # Sx
        sm_l0_sx = SpikeMonitor(l0_sx, variables='vm')
        net.add(sm_l0_sx)
        if num_layers > 1:
            sm_l1_sx = SpikeMonitor(l1_sx, variables='vm')
            net.add(sm_l1_sx)
        if num_layers > 2:
            sm_l2_sx = SpikeMonitor(l2_sx, variables='vm')
            net.add(sm_l2_sx)
        if num_layers > 3:
            sm_l3_sx = SpikeMonitor(l3_sx, variables='vm')
            net.add(sm_l3_sx)
        if num_layers > 4:
            sm_l4_sx = SpikeMonitor(l4_sx, variables='vm')
            net.add(sm_l4_sx)
        if num_layers > 5:
            sm_l5_sx = SpikeMonitor(l5_sx, variables='vm')
            net.add(sm_l5_sx)

        # Dx
        sm_l0_dx = SpikeMonitor(l0_dx, variables='vm')
        net.add(sm_l0_dx)
        if num_layers > 1:
            sm_l1_dx = SpikeMonitor(l1_dx, variables='vm')
            net.add(sm_l1_dx)
        if num_layers > 2:
            sm_l2_dx = SpikeMonitor(l2_dx, variables='vm')
            net.add(sm_l2_dx)
        if num_layers > 3:
            sm_l3_dx = SpikeMonitor(l3_dx, variables='vm')
            net.add(sm_l3_dx)
        if num_layers > 4:
            sm_l4_dx = SpikeMonitor(l4_dx, variables='vm')
            net.add(sm_l4_dx)
        if num_layers > 5:
            sm_l5_dx = SpikeMonitor(l5_dx, variables='vm')
            net.add(sm_l5_dx)

        net.run(tt, report='text')

        # Sx
        train_l0_sx = sm_l0_sx.spike_trains()
        savepeaktrain('ptl0_sx', train_l0_sx, l0_sx, n_exc, n_inh, 0)
        del sm_l0_sx
        if num_layers > 1:
            train_l1_sx = sm_l1_sx.spike_trains()
            savepeaktrain('ptl1_sx', train_l1_sx, l1_sx, n_exc, n_inh, 1)
            del sm_l1_sx
        if num_layers > 2:
            train_l2_sx = sm_l2_sx.spike_trains()
            savepeaktrain('ptl2_sx', train_l2_sx, l2_sx, n_exc, n_inh, 2)
            del sm_l2_sx
        if num_layers > 3:
            train_l3_sx = sm_l3_sx.spike_trains()
            savepeaktrain('ptl3_sx', train_l3_sx, l3_sx, n_exc, n_inh, 3)
            del sm_l3_sx
        if num_layers > 4:
            train_l4_sx = sm_l4_sx.spike_trains()
            savepeaktrain('ptl4_sx', train_l4_sx, l4_sx, n_exc, n_inh, 4)
            del sm_l4_sx
        if num_layers > 5:
            train_l5_sx = sm_l5_sx.spike_trains()
            savepeaktrain('ptl5_sx', train_l5_sx, l5_sx, n_exc, n_inh, 5)
            del sm_l5_sx

        # Dx
        train_l0_dx = sm_l0_dx.spike_trains()
        savepeaktrain('ptl0_dx', train_l0_dx, l0_dx, n_exc, n_inh, 0)
        del sm_l0_dx
        if num_layers > 1:
            train_l1_dx = sm_l1_dx.spike_trains()
            savepeaktrain('ptl1_dx', train_l1_dx, l1_dx, n_exc, n_inh, 1)
            del sm_l1_dx
        if num_layers > 2:
            train_l2_dx = sm_l2_dx.spike_trains()
            savepeaktrain('ptl2_dx', train_l2_dx, l2_dx, n_exc, n_inh, 2)
            del sm_l2_dx
        if num_layers > 3:
            train_l3_dx = sm_l3_dx.spike_trains()
            savepeaktrain('ptl3_dx', train_l3_dx, l3_dx, n_exc, n_inh, 3)
            del sm_l3_dx
        if num_layers > 4:
            train_l4_dx = sm_l4_dx.spike_trains()
            savepeaktrain('ptl4_dx', train_l4_dx, l4_dx, n_exc, n_inh, 4)
            del sm_l4_dx
        if num_layers > 5:
            train_l5_dx = sm_l5_dx.spike_trains()
            savepeaktrain('ptl5_dx', train_l5_dx, l5_dx, n_exc, n_inh, 5)
            del sm_l5_dx

    os.chdir(SimDir)
else:
    net.run(simtime, report='text')
#####################################


#####################################
#         Results display
#####################################

# ------ IMAGES ------
if should_save:
    os.chdir(ImagesDir)

# ------ Raster plot ------
# Sx
plot_raster(sm_layer0_sx, preptime, simtime, 'raster_l0_sx.tif', s=should_show, save=should_save)
if num_layers > 1:
    plot_raster(sm_layer1_sx, preptime, simtime, 'raster_l1_sx.tif', s=should_show, save=should_save)
if num_layers > 2:
    plot_raster(sm_layer2_sx, preptime, simtime, 'raster_l2_sx.tif', s=should_show, save=should_save)
if num_layers > 3:
    plot_raster(sm_layer3_sx, preptime, simtime, 'raster_l3_sx.tif', s=should_show, save=should_save)
if num_layers > 4:
    plot_raster(sm_layer4_sx, preptime, simtime, 'raster_l4_sx.tif', s=should_show, save=should_save)
if num_layers > 5:
    plot_raster(sm_layer5_sx, preptime, simtime, 'raster_l5_sx.tif', s=should_show, save=should_save)

# Dx
plot_raster(sm_layer0_dx, preptime, simtime, 'raster_l0_dx.tif', s=should_show, save=should_save)
if num_layers > 1:
    plot_raster(sm_layer1_dx, preptime, simtime, 'raster_l1_dx.tif', s=should_show, save=should_save)
if num_layers > 2:
    plot_raster(sm_layer2_dx, preptime, simtime, 'raster_l2_dx.tif', s=should_show, save=should_save)
if num_layers > 3:
    plot_raster(sm_layer3_dx, preptime, simtime, 'raster_l3_dx.tif', s=should_show, save=should_save)
if num_layers > 4:
    plot_raster(sm_layer4_dx, preptime, simtime, 'raster_l4_dx.tif', s=should_show, save=should_save)
if num_layers > 5:
    plot_raster(sm_layer5_dx, preptime, simtime, 'raster_l5_dx.tif', s=should_show, save=should_save)

# ------ Shuffled raster plot ------
# Sx
plot_raster_shuffle(l0_sx, sm_layer0_sx, preptime, simtime, 'raster_shuffle_l0_sx.tif', s=should_show, save=should_save)
if num_layers > 1:
    plot_raster_shuffle(l1_sx, sm_layer1_sx, preptime, simtime, 'raster_shuffle_l1_sx.tif', s=should_show,
                        save=should_save)
if num_layers > 2:
    plot_raster_shuffle(l2_sx, sm_layer2_sx, preptime, simtime, 'raster_shuffle_l2_sx.tif', s=should_show,
                        save=should_save)
if num_layers > 3:
    plot_raster_shuffle(l3_sx, sm_layer3_sx, preptime, simtime, 'raster_shuffle_l3_sx.tif', s=should_show,
                        save=should_save)
if num_layers > 4:
    plot_raster_shuffle(l4_sx, sm_layer4_sx, preptime, simtime, 'raster_shuffle_l4_sx.tif', s=should_show,
                        save=should_save)
if num_layers > 5:
    plot_raster_shuffle(l5_sx, sm_layer5_sx, preptime, simtime, 'raster_shuffle_l5_sx.tif', s=should_show,
                        save=should_save)

# Dx
plot_raster_shuffle(l0_dx, sm_layer0_dx, preptime, simtime, 'raster_shuffle_l0_dx.tif', s=should_show, save=should_save)
if num_layers > 1:
    plot_raster_shuffle(l1_dx, sm_layer1_dx, preptime, simtime, 'raster_shuffle_l1_dx.tif', s=should_show,
                        save=should_save)
if num_layers > 2:
    plot_raster_shuffle(l2_dx, sm_layer2_dx, preptime, simtime, 'raster_shuffle_l2_dx.tif', s=should_show,
                        save=should_save)
if num_layers > 3:
    plot_raster_shuffle(l3_dx, sm_layer3_dx, preptime, simtime, 'raster_shuffle_l3_dx.tif', s=should_show,
                        save=should_save)
if num_layers > 4:
    plot_raster_shuffle(l4_dx, sm_layer4_dx, preptime, simtime, 'raster_shuffle_l4_dx.tif', s=should_show,
                        save=should_save)
if num_layers > 5:
    plot_raster_shuffle(l5_dx, sm_layer5_dx, preptime, simtime, 'raster_shuffle_l5_dx.tif', s=should_show,
                        save=should_save)

# ------ RASTER shuffle ------
fig1 = figure(1, figsize=[6 * num_layers * cm_fig_conv, 14 * cm_fig_conv], dpi=600.0)

ax0_sx = subplot(2, num_layers, 1)
pt_l0 = sm_layer0_sx.spike_trains()
for key, __ in pt_l0.items():
    plot(pt_l0[key], repeat(l0_sx.index[key], len(pt_l0[key])), ',k')
xlim(preptime / second, (preptime + simtime) / second)

ax0_dx = subplot(2, num_layers, num_layers + 1)
pt_l0 = sm_layer0_dx.spike_trains()
for key, __ in pt_l0.items():
    plot(pt_l0[key], repeat(l0_dx.index[key], len(pt_l0[key])), ',k')
xlim(preptime / second, (preptime + simtime) / second)

if num_layers > 1:
    ax1_sx = subplot(2, num_layers, 2, sharex=ax0_sx)
    pt_l1 = sm_layer1_sx.spike_trains()
    for key, __ in pt_l1.items():
        plot(pt_l1[key], repeat(l1_sx.index[key], len(pt_l1[key])), ',k')
    xlim(preptime / second, (preptime + simtime) / second)

    ax1_dx = subplot(2, num_layers, num_layers + 2, sharex=ax0_sx)
    pt_l1 = sm_layer1_dx.spike_trains()
    for key, __ in pt_l1.items():
        plot(pt_l1[key], repeat(l1_dx.index[key], len(pt_l1[key])), ',k')
    xlim(preptime / second, (preptime + simtime) / second)
if num_layers > 2:
    ax2_sx = subplot(2, num_layers, 3, sharex=ax0_sx)
    pt_l2 = sm_layer2_sx.spike_trains()
    for key, __ in pt_l2.items():
        plot(pt_l2[key], repeat(l2_sx.index[key], len(pt_l2[key])), ',k')
    xlim(preptime / second, (preptime + simtime) / second)

    ax2_dx = subplot(2, num_layers, num_layers + 3, sharex=ax0_dx)
    pt_l2 = sm_layer2_dx.spike_trains()
    for key, __ in pt_l2.items():
        plot(pt_l2[key], repeat(l2_dx.index[key], len(pt_l2[key])), ',k')
    xlim(preptime / second, (preptime + simtime) / second)
if num_layers > 3:
    ax3_sx = subplot(2, num_layers, 4, sharex=ax0_sx)
    pt_l3 = sm_layer3_sx.spike_trains()
    for key, __ in pt_l3.items():
        plot(pt_l3[key], repeat(l3_sx.index[key], len(pt_l3[key])), ',k')
    xlim(preptime / second, (preptime + simtime) / second)

    ax3_dx = subplot(2, num_layers, num_layers + 4, sharex=ax0_dx)
    pt_l3 = sm_layer3_dx.spike_trains()
    for key, __ in pt_l3.items():
        plot(pt_l3[key], repeat(l3_dx.index[key], len(pt_l3[key])), ',k')
    xlim(preptime / second, (preptime + simtime) / second)
if num_layers > 4:
    ax4_sx = subplot(2, num_layers, 5, sharex=ax0_sx)
    pt_l4 = sm_layer4_sx.spike_trains()
    for key, __ in pt_l4.items():
        plot(pt_l4[key], repeat(l4_sx.index[key], len(pt_l4[key])), ',k')
    xlim(preptime / second, (preptime + simtime) / second)

    ax4_dx = subplot(2, num_layers, num_layers + 5, sharex=ax0_dx)
    pt_l4 = sm_layer4_dx.spike_trains()
    for key, __ in pt_l4.items():
        plot(pt_l4[key], repeat(l4_dx.index[key], len(pt_l4[key])), ',k')
    xlim(preptime / second, (preptime + simtime) / second)
if num_layers > 5:
    ax5_sx = subplot(2, num_layers, 6, sharex=ax0_sx)
    pt_l5 = sm_layer5_sx.spike_trains()
    for key, __ in pt_l5.items():
        plot(pt_l5[key], repeat(l5_sx.index[key], len(pt_l5[key])), ',k')
    xlim(preptime / second, (preptime + simtime) / second)

    ax5_dx = subplot(2, num_layers, num_layers + 6, sharex=ax0_dx)
    pt_l5 = sm_layer5_dx.spike_trains()
    for key, __ in pt_l5.items():
        plot(pt_l5[key], repeat(l5_dx.index[key], len(pt_l5[key])), ',k')
    xlim(preptime / second, (preptime + simtime) / second)

xlabel('Time (s)')
ylabel('Neuron index (position)', rotation='vertical')
gca().spines['top'].set_visible(False)
gca().spines['right'].set_visible(False)
tight_layout()

if should_show:
    show()

if should_save:
    fig1.savefig('Raster_complessivo.tif', dpi=600, transparent=True)

# ------ 3D layout ------
plot_spatial_configuration(layers_exc, layers_inh, s=should_show, save=should_save)

# os.chdir(SimDir)

# ------ VALUES ------
mfr_l0_sx = mean(sm_layer0_sx.count / simtime)
print('MFR_l0_sx: ' + str(mfr_l0_sx))
mfr_l0_dx = mean(sm_layer0_dx.count / simtime)
print('MFR_l0_dx: ' + str(mfr_l0_dx))
if num_layers > 1:
    mfr_l1_sx = mean(sm_layer1_sx.count / simtime)
    print('MFR_l1_sx: ' + str(mfr_l1_sx))
    mfr_l1_dx = mean(sm_layer1_dx.count / simtime)
    print('MFR_l1_dx: ' + str(mfr_l1_dx))
if num_layers > 2:
    mfr_l2_sx = mean(sm_layer2_sx.count / simtime)
    print('MFR_l2_sx: ' + str(mfr_l2_sx))
    mfr_l2_dx = mean(sm_layer2_dx.count / simtime)
    print('MFR_l2_dx: ' + str(mfr_l2_dx))
if num_layers > 3:
    mfr_l3_sx = mean(sm_layer3_sx.count / simtime)
    print('MFR_l3_sx: ' + str(mfr_l3_sx))
    mfr_l3_dx = mean(sm_layer3_dx.count / simtime)
    print('MFR_l3_dx: ' + str(mfr_l3_dx))
if num_layers > 4:
    mfr_l4_sx = mean(sm_layer4_sx.count / simtime)
    print('MFR_l4_sx: ' + str(mfr_l4_sx))
    mfr_l4_dx = mean(sm_layer4_dx.count / simtime)
    print('MFR_l4_dx: ' + str(mfr_l4_dx))
if num_layers > 5:
    mfr_l5_sx = mean(sm_layer5_sx.count / simtime)
    print('MFR_l5_sx: ' + str(mfr_l5_sx))
    mfr_l5_dx = mean(sm_layer5_dx.count / simtime)
    print('MFR_l5_dx: ' + str(mfr_l5_dx))
#####################################


#####################################
#         Saving Results
#####################################
if should_save:
    os.chdir(SimDir)

    # --- Saving the specifics of the simulation: duration, frequency, input equation, MFR, etc. ---
    with open("SimulationSpec", "a+") as f:
        f.write(DateFolder + '\tSimulation ' + str(SimulationNumber) + '\n')
        f.write('\nDuration: \t%s' % simtime)
        f.write('\nPreparation Time: \t%s' % preptime)
        f.write('\nSampling frequency: \t%s' % (1 / defaultclock.dt))
        f.write('\nSeed: \t%s' % val)
        f.write('\nNumber of:'
                '\n\t- Cells per layer \t%s'
                '\n\t- Exc cells \t%s'
                '\n\t- Inh cells \t%s'
                % (n_cells, n_exc, n_inh))
        f.write('\n\nInput value: ' + input_value)
        f.write('\n\nMFR_l0_sx: %s' % mfr_l0_sx)
        f.write('\nMFR_l0_dx: %s' % mfr_l0_dx)
        if num_layers > 1:
            f.write('\nMFR_l1_dx: %s' % mfr_l1_dx)
            f.write('\nMFR_l1_dx: %s' % mfr_l1_dx)
        if num_layers > 2:
            f.write('\nMFR_l2_sx: %s' % mfr_l2_sx)
            f.write('\nMFR_l2_dx: %s' % mfr_l2_dx)
        if num_layers > 3:
            f.write('\nMFR_l3_sx: %s' % mfr_l3_sx)
            f.write('\nMFR_l3_dx: %s' % mfr_l3_dx)
        if num_layers > 4:
            f.write('\nMFR_l4_dx: %s' % mfr_l4_dx)
            f.write('\nMFR_l4_dx: %s' % mfr_l4_dx)
        if num_layers > 5:
            f.write('\nMFR_l5_sx: %s' % mfr_l5_sx)
            f.write('\nMFR_l5_dx: %s' % mfr_l5_dx)
        f.write('\n\n' + comment)
        f.write('\n\n\n------------------------------\n\n\n')
        f.close()

    # --- Saving the lookup table ---
    os.chdir(LutDir)
    savelookuptable('lookuptable_l0_sx', train_l0_sx, l0_sx)
    savelookuptable('lookuptable_l0_dx', train_l0_dx, l0_dx)
    if num_layers > 1:
        savelookuptable('lookuptable_l1_sx', train_l1_sx, l1_sx)
        savelookuptable('lookuptable_l1_dx', train_l1_dx, l1_dx)
    if num_layers > 2:
        savelookuptable('lookuptable_l2_sx', train_l2_sx, l2_sx)
        savelookuptable('lookuptable_l2_dx', train_l2_dx, l2_dx)
    if num_layers > 3:
        savelookuptable('lookuptable_l3_sx', train_l3_sx, l3_sx)
        savelookuptable('lookuptable_l3_dx', train_l3_dx, l3_dx)
    if num_layers > 4:
        savelookuptable('lookuptable_l4_sx', train_l4_sx, l4_sx)
        savelookuptable('lookuptable_l4_dx', train_l4_dx, l4_dx)
    if num_layers > 5:
        savelookuptable('lookuptable_l5_sx', train_l5_sx, l5_sx)
        savelookuptable('lookuptable_l5_dx', train_l5_dx, l5_dx)
    os.chdir(SimDir)

    # --- Saving the weight matrices ---
    os.chdir(WeightDir)
    np.savetxt('WeightMatrix_s_exc_all_l0_sx.txt', WeightMatrix_s_exc_all_l0_sx)
    np.savetxt('WeightMatrix_s_inh_exc_l0_sx.txt', WeightMatrix_s_inh_exc_l0_sx)
    np.savetxt('WeightMatrix_s_exc_all_l0_dx.txt', WeightMatrix_s_exc_all_l0_dx)
    np.savetxt('WeightMatrix_s_inh_exc_l0_dx.txt', WeightMatrix_s_inh_exc_l0_dx)
    np.savetxt('WeightMatrix_s_exc_all_l0_sxdx.txt', WeightMatrix_s_exc_all_l0_sxdx)
    np.savetxt('WeightMatrix_s_inh_exc_l0_sxdx.txt', WeightMatrix_s_inh_exc_l0_sxdx)
    np.savetxt('WeightMatrix_s_exc_all_l0_dxsx.txt', WeightMatrix_s_exc_all_l0_dxsx)
    np.savetxt('WeightMatrix_s_inh_exc_l0_dxsx.txt', WeightMatrix_s_inh_exc_l0_dxsx)
    if num_layers > 1:
        np.savetxt('WeightMatrix_s_exc_all_l1_sx.txt', WeightMatrix_s_exc_all_l1_sx)
        np.savetxt('WeightMatrix_s_inh_exc_l1_sx.txt', WeightMatrix_s_inh_exc_l1_sx)
        np.savetxt('WeightMatrix_s_exc_all_l1_dx.txt', WeightMatrix_s_exc_all_l1_dx)
        np.savetxt('WeightMatrix_s_inh_exc_l1_dx.txt', WeightMatrix_s_inh_exc_l1_dx)
        if upward_conn:
            # Layer 0 to Layer 1
            np.savetxt('WeightMatrix_s_exc_all_l0l1_sx.txt', WeightMatrix_s_exc_all_l0l1_sx)
            np.savetxt('WeightMatrix_s_inh_exc_l0l1_sx.txt', WeightMatrix_s_inh_exc_l0l1_sx)
            np.savetxt('WeightMatrix_s_exc_all_l0l1_dx.txt', WeightMatrix_s_exc_all_l0l1_dx)
            np.savetxt('WeightMatrix_s_inh_exc_l0l1_dx.txt', WeightMatrix_s_inh_exc_l0l1_dx)
        if downward_conn:
            # Layer 1 to Layer 0
            np.savetxt('WeightMatrix_s_exc_all_l1l0_sx.txt', WeightMatrix_s_exc_all_l1l0_sx)
            np.savetxt('WeightMatrix_s_inh_exc_l1l0_sx.txt', WeightMatrix_s_inh_exc_l1l0_sx)
            np.savetxt('WeightMatrix_s_exc_all_l1l0_dx.txt', WeightMatrix_s_exc_all_l1l0_dx)
            np.savetxt('WeightMatrix_s_inh_exc_l1l0_dx.txt', WeightMatrix_s_inh_exc_l1l0_dx)
    if num_layers > 2:
        np.savetxt('WeightMatrix_s_exc_all_l2_sx.txt', WeightMatrix_s_exc_all_l2_sx)
        np.savetxt('WeightMatrix_s_inh_exc_l2_sx.txt', WeightMatrix_s_inh_exc_l2_sx)
        np.savetxt('WeightMatrix_s_exc_all_l2_dx.txt', WeightMatrix_s_exc_all_l2_dx)
        np.savetxt('WeightMatrix_s_inh_exc_l2_dx.txt', WeightMatrix_s_inh_exc_l2_dx)
        if upward_conn:
            # Layer 0 to Layer 2
            np.savetxt('WeightMatrix_s_exc_all_l0l2_sx.txt', WeightMatrix_s_exc_all_l0l2_sx)
            np.savetxt('WeightMatrix_s_inh_exc_l0l2_sx.txt', WeightMatrix_s_inh_exc_l0l2_sx)
            np.savetxt('WeightMatrix_s_exc_all_l0l2_dx.txt', WeightMatrix_s_exc_all_l0l2_dx)
            np.savetxt('WeightMatrix_s_inh_exc_l0l2_dx.txt', WeightMatrix_s_inh_exc_l0l2_dx)
            # Layer 1 to Layer 2
            np.savetxt('WeightMatrix_s_exc_all_l1l2_sx.txt', WeightMatrix_s_exc_all_l1l2_sx)
            np.savetxt('WeightMatrix_s_inh_exc_l1l2_sx.txt', WeightMatrix_s_inh_exc_l1l2_sx)
            np.savetxt('WeightMatrix_s_exc_all_l1l2_dx.txt', WeightMatrix_s_exc_all_l1l2_dx)
            np.savetxt('WeightMatrix_s_inh_exc_l1l2_dx.txt', WeightMatrix_s_inh_exc_l1l2_dx)
        if downward_conn:
            # Layer 2 to Layer 0
            np.savetxt('WeightMatrix_s_exc_all_l2l0_sx.txt', WeightMatrix_s_exc_all_l2l0_sx)
            np.savetxt('WeightMatrix_s_inh_exc_l2l0_sx.txt', WeightMatrix_s_inh_exc_l2l0_sx)
            np.savetxt('WeightMatrix_s_exc_all_l2l0_dx.txt', WeightMatrix_s_exc_all_l2l0_dx)
            np.savetxt('WeightMatrix_s_inh_exc_l2l0_dx.txt', WeightMatrix_s_inh_exc_l2l0_dx)
            # Layer 2 to Layer 1
            np.savetxt('WeightMatrix_s_exc_all_l2l1_sx.txt', WeightMatrix_s_exc_all_l2l1_sx)
            np.savetxt('WeightMatrix_s_inh_exc_l2l1_sx.txt', WeightMatrix_s_inh_exc_l2l1_sx)
            np.savetxt('WeightMatrix_s_exc_all_l2l1_dx.txt', WeightMatrix_s_exc_all_l2l1_dx)
            np.savetxt('WeightMatrix_s_inh_exc_l2l1_dx.txt', WeightMatrix_s_inh_exc_l2l1_dx)
    if num_layers > 3:
        np.savetxt('WeightMatrix_s_exc_all_l3_sx.txt', WeightMatrix_s_exc_all_l3_sx)
        np.savetxt('WeightMatrix_s_inh_exc_l3_sx.txt', WeightMatrix_s_inh_exc_l3_sx)
        np.savetxt('WeightMatrix_s_exc_all_l3_dx.txt', WeightMatrix_s_exc_all_l3_dx)
        np.savetxt('WeightMatrix_s_inh_exc_l3_dx.txt', WeightMatrix_s_inh_exc_l3_dx)
        if upward_conn:
            # Layer 0 to Layer 3
            np.savetxt('WeightMatrix_s_exc_all_l0l3_sx.txt', WeightMatrix_s_exc_all_l0l3_sx)
            np.savetxt('WeightMatrix_s_inh_exc_l0l3_sx.txt', WeightMatrix_s_inh_exc_l0l3_sx)
            np.savetxt('WeightMatrix_s_exc_all_l0l3_dx.txt', WeightMatrix_s_exc_all_l0l3_dx)
            np.savetxt('WeightMatrix_s_inh_exc_l0l3_dx.txt', WeightMatrix_s_inh_exc_l0l3_dx)
            # Layer 1 to Layer 3
            np.savetxt('WeightMatrix_s_exc_all_l1l3_sx.txt', WeightMatrix_s_exc_all_l1l3_sx)
            np.savetxt('WeightMatrix_s_inh_exc_l1l3_sx.txt', WeightMatrix_s_inh_exc_l1l3_sx)
            np.savetxt('WeightMatrix_s_exc_all_l1l3_dx.txt', WeightMatrix_s_exc_all_l1l3_dx)
            np.savetxt('WeightMatrix_s_inh_exc_l1l3_dx.txt', WeightMatrix_s_inh_exc_l1l3_dx)
            # Layer 2 to Layer 3
            np.savetxt('WeightMatrix_s_exc_all_l2l3_sx.txt', WeightMatrix_s_exc_all_l2l3_sx)
            np.savetxt('WeightMatrix_s_inh_exc_l2l3_sx.txt', WeightMatrix_s_inh_exc_l2l3_sx)
            np.savetxt('WeightMatrix_s_exc_all_l2l3_dx.txt', WeightMatrix_s_exc_all_l2l3_dx)
            np.savetxt('WeightMatrix_s_inh_exc_l2l3_dx.txt', WeightMatrix_s_inh_exc_l2l3_dx)
        if downward_conn:
            # Layer 3 to Layer 0
            np.savetxt('WeightMatrix_s_exc_all_l3l0_sx.txt', WeightMatrix_s_exc_all_l3l0_sx)
            np.savetxt('WeightMatrix_s_inh_exc_l3l0_sx.txt', WeightMatrix_s_inh_exc_l3l0_sx)
            np.savetxt('WeightMatrix_s_exc_all_l3l0_dx.txt', WeightMatrix_s_exc_all_l3l0_dx)
            np.savetxt('WeightMatrix_s_inh_exc_l3l0_dx.txt', WeightMatrix_s_inh_exc_l3l0_dx)
            # Layer 3 to Layer 1
            np.savetxt('WeightMatrix_s_exc_all_l3l1_sx.txt', WeightMatrix_s_exc_all_l3l1_sx)
            np.savetxt('WeightMatrix_s_inh_exc_l3l1_sx.txt', WeightMatrix_s_inh_exc_l3l1_sx)
            np.savetxt('WeightMatrix_s_exc_all_l3l1_dx.txt', WeightMatrix_s_exc_all_l3l1_dx)
            np.savetxt('WeightMatrix_s_inh_exc_l3l1_dx.txt', WeightMatrix_s_inh_exc_l3l1_dx)
            # Layer 3 to Layer 2
            np.savetxt('WeightMatrix_s_exc_all_l3l2_sx.txt', WeightMatrix_s_exc_all_l3l2_sx)
            np.savetxt('WeightMatrix_s_inh_exc_l3l2_sx.txt', WeightMatrix_s_inh_exc_l3l2_sx)
            np.savetxt('WeightMatrix_s_exc_all_l3l2_dx.txt', WeightMatrix_s_exc_all_l3l2_dx)
            np.savetxt('WeightMatrix_s_inh_exc_l3l2_dx.txt', WeightMatrix_s_inh_exc_l3l2_dx)
    if num_layers > 4:
        np.savetxt('WeightMatrix_s_exc_all_l4_sx.txt', WeightMatrix_s_exc_all_l4_sx)
        np.savetxt('WeightMatrix_s_inh_exc_l4_sx.txt', WeightMatrix_s_inh_exc_l4_sx)
        np.savetxt('WeightMatrix_s_exc_all_l4_dx.txt', WeightMatrix_s_exc_all_l4_dx)
        np.savetxt('WeightMatrix_s_inh_exc_l4_dx.txt', WeightMatrix_s_inh_exc_l4_dx)
        if upward_conn:
            # Layer 0 to Layer 4
            np.savetxt('WeightMatrix_s_exc_all_l0l4_sx.txt', WeightMatrix_s_exc_all_l0l4_sx)
            np.savetxt('WeightMatrix_s_inh_exc_l0l4_sx.txt', WeightMatrix_s_inh_exc_l0l4_sx)
            np.savetxt('WeightMatrix_s_exc_all_l0l4_dx.txt', WeightMatrix_s_exc_all_l0l4_dx)
            np.savetxt('WeightMatrix_s_inh_exc_l0l4_dx.txt', WeightMatrix_s_inh_exc_l0l4_dx)
            # Layer 1 to Layer 4
            np.savetxt('WeightMatrix_s_exc_all_l1l4_sx.txt', WeightMatrix_s_exc_all_l1l4_sx)
            np.savetxt('WeightMatrix_s_inh_exc_l1l4_sx.txt', WeightMatrix_s_inh_exc_l1l4_sx)
            # Layer 2 to Layer 4
            np.savetxt('WeightMatrix_s_exc_all_l2l4_sx.txt', WeightMatrix_s_exc_all_l2l4_sx)
            np.savetxt('WeightMatrix_s_inh_exc_l2l4_sx.txt', WeightMatrix_s_inh_exc_l2l4_sx)
            np.savetxt('WeightMatrix_s_exc_all_l2l4_dx.txt', WeightMatrix_s_exc_all_l2l4_dx)
            np.savetxt('WeightMatrix_s_inh_exc_l2l4_dx.txt', WeightMatrix_s_inh_exc_l2l4_dx)
            # Layer 3 to Layer 4
            np.savetxt('WeightMatrix_s_exc_all_l3l4_sx.txt', WeightMatrix_s_exc_all_l3l4_sx)
            np.savetxt('WeightMatrix_s_inh_exc_l3l4_sx.txt', WeightMatrix_s_inh_exc_l3l4_sx)
            np.savetxt('WeightMatrix_s_exc_all_l3l4_dx.txt', WeightMatrix_s_exc_all_l3l4_dx)
            np.savetxt('WeightMatrix_s_inh_exc_l3l4_dx.txt', WeightMatrix_s_inh_exc_l3l4_dx)
        if downward_conn:
            # Layer 4 to Layer 0
            np.savetxt('WeightMatrix_s_exc_all_l4l0_sx.txt', WeightMatrix_s_exc_all_l4l0_sx)
            np.savetxt('WeightMatrix_s_inh_exc_l4l0_sx.txt', WeightMatrix_s_inh_exc_l4l0_sx)
            np.savetxt('WeightMatrix_s_exc_all_l4l0_dx.txt', WeightMatrix_s_exc_all_l4l0_dx)
            np.savetxt('WeightMatrix_s_inh_exc_l4l0_dx.txt', WeightMatrix_s_inh_exc_l4l0_dx)
            # Layer 4 to Layer 1
            np.savetxt('WeightMatrix_s_exc_all_l4l1_sx.txt', WeightMatrix_s_exc_all_l4l1_sx)
            np.savetxt('WeightMatrix_s_inh_exc_l4l1_sx.txt', WeightMatrix_s_inh_exc_l4l1_sx)
            np.savetxt('WeightMatrix_s_exc_all_l4l1_dx.txt', WeightMatrix_s_exc_all_l4l1_dx)
            np.savetxt('WeightMatrix_s_inh_exc_l4l1_dx.txt', WeightMatrix_s_inh_exc_l4l1_dx)
            # Layer 4 to Layer 2
            np.savetxt('WeightMatrix_s_exc_all_l4l2_sx.txt', WeightMatrix_s_exc_all_l4l2_sx)
            np.savetxt('WeightMatrix_s_inh_exc_l4l2_sx.txt', WeightMatrix_s_inh_exc_l4l2_sx)
            np.savetxt('WeightMatrix_s_exc_all_l4l2_dx.txt', WeightMatrix_s_exc_all_l4l2_dx)
            np.savetxt('WeightMatrix_s_inh_exc_l4l2_dx.txt', WeightMatrix_s_inh_exc_l4l2_dx)
            # Layer 4 to Layer 3
            np.savetxt('WeightMatrix_s_exc_all_l4l3_sx.txt', WeightMatrix_s_exc_all_l4l3_sx)
            np.savetxt('WeightMatrix_s_inh_exc_l4l3_sx.txt', WeightMatrix_s_inh_exc_l4l3_sx)
            np.savetxt('WeightMatrix_s_exc_all_l4l3_dx.txt', WeightMatrix_s_exc_all_l4l3_dx)
            np.savetxt('WeightMatrix_s_inh_exc_l4l3_dx.txt', WeightMatrix_s_inh_exc_l4l3_dx)
    if num_layers > 5:
        np.savetxt('WeightMatrix_s_exc_all_l5_sx.txt', WeightMatrix_s_exc_all_l5_sx)
        np.savetxt('WeightMatrix_s_inh_exc_l5_sx.txt', WeightMatrix_s_inh_exc_l5_sx)
        np.savetxt('WeightMatrix_s_exc_all_l5_dx.txt', WeightMatrix_s_exc_all_l5_dx)
        np.savetxt('WeightMatrix_s_inh_exc_l5_dx.txt', WeightMatrix_s_inh_exc_l5_dx)
        if upward_conn:
            # Layer 0 to Layer 5
            np.savetxt('WeightMatrix_s_exc_all_l0l5_sx.txt', WeightMatrix_s_exc_all_l0l5_sx)
            np.savetxt('WeightMatrix_s_inh_exc_l0l5_sx.txt', WeightMatrix_s_inh_exc_l0l5_sx)
            np.savetxt('WeightMatrix_s_exc_all_l0l5_dx.txt', WeightMatrix_s_exc_all_l0l5_dx)
            np.savetxt('WeightMatrix_s_inh_exc_l0l5_dx.txt', WeightMatrix_s_inh_exc_l0l5_dx)
            # Layer 1 to Layer 5
            np.savetxt('WeightMatrix_s_exc_all_l1l5_sx.txt', WeightMatrix_s_exc_all_l1l5_sx)
            np.savetxt('WeightMatrix_s_inh_exc_l1l5_sx.txt', WeightMatrix_s_inh_exc_l1l5_sx)
            np.savetxt('WeightMatrix_s_exc_all_l1l5_dx.txt', WeightMatrix_s_exc_all_l1l5_dx)
            np.savetxt('WeightMatrix_s_inh_exc_l1l5_dx.txt', WeightMatrix_s_inh_exc_l1l5_dx)
            # Layer 2 to Layer 5
            np.savetxt('WeightMatrix_s_exc_all_l2l5_sx.txt', WeightMatrix_s_exc_all_l2l5_sx)
            np.savetxt('WeightMatrix_s_inh_exc_l2l5_sx.txt', WeightMatrix_s_inh_exc_l2l5_sx)
            np.savetxt('WeightMatrix_s_exc_all_l2l5_dx.txt', WeightMatrix_s_exc_all_l2l5_dx)
            np.savetxt('WeightMatrix_s_inh_exc_l2l5_dx.txt', WeightMatrix_s_inh_exc_l2l5_dx)
            # Layer 3 to Layer 5
            np.savetxt('WeightMatrix_s_exc_all_l3l5_sx.txt', WeightMatrix_s_exc_all_l3l5_sx)
            np.savetxt('WeightMatrix_s_inh_exc_l3l5_sx.txt', WeightMatrix_s_inh_exc_l3l5_sx)
            np.savetxt('WeightMatrix_s_exc_all_l3l5_dx.txt', WeightMatrix_s_exc_all_l3l5_dx)
            np.savetxt('WeightMatrix_s_inh_exc_l3l5_dx.txt', WeightMatrix_s_inh_exc_l3l5_dx)
            # Layer 4 to Layer 5
            np.savetxt('WeightMatrix_s_exc_all_l4l5_sx.txt', WeightMatrix_s_exc_all_l4l5_sx)
            np.savetxt('WeightMatrix_s_inh_exc_l4l5_sx.txt', WeightMatrix_s_inh_exc_l4l5_sx)
            np.savetxt('WeightMatrix_s_exc_all_l4l5_dx.txt', WeightMatrix_s_exc_all_l4l5_dx)
            np.savetxt('WeightMatrix_s_inh_exc_l4l5_dx.txt', WeightMatrix_s_inh_exc_l4l5_dx)
        if downward_conn:
            # Layer 5 to Layer 0
            np.savetxt('WeightMatrix_s_exc_all_l5l0_sx.txt', WeightMatrix_s_exc_all_l5l0_sx)
            np.savetxt('WeightMatrix_s_inh_exc_l5l0_sx.txt', WeightMatrix_s_inh_exc_l5l0_sx)
            np.savetxt('WeightMatrix_s_exc_all_l5l0_dx.txt', WeightMatrix_s_exc_all_l5l0_dx)
            np.savetxt('WeightMatrix_s_inh_exc_l5l0_dx.txt', WeightMatrix_s_inh_exc_l5l0_dx)
            # Layer 5 to Layer 1
            np.savetxt('WeightMatrix_s_exc_all_l5l1_sx.txt', WeightMatrix_s_exc_all_l5l1_sx)
            np.savetxt('WeightMatrix_s_inh_exc_l5l1_sx.txt', WeightMatrix_s_inh_exc_l5l1_sx)
            np.savetxt('WeightMatrix_s_exc_all_l5l1_dx.txt', WeightMatrix_s_exc_all_l5l1_dx)
            np.savetxt('WeightMatrix_s_inh_exc_l5l1_dx.txt', WeightMatrix_s_inh_exc_l5l1_dx)
            # Layer 5 to Layer 2
            np.savetxt('WeightMatrix_s_exc_all_l5l2_sx.txt', WeightMatrix_s_exc_all_l5l2_sx)
            np.savetxt('WeightMatrix_s_inh_exc_l5l2_sx.txt', WeightMatrix_s_inh_exc_l5l2_sx)
            np.savetxt('WeightMatrix_s_exc_all_l5l2_dx.txt', WeightMatrix_s_exc_all_l5l2_dx)
            np.savetxt('WeightMatrix_s_inh_exc_l5l2_dx.txt', WeightMatrix_s_inh_exc_l5l2_dx)
            # Layer 5 to Layer 3
            np.savetxt('WeightMatrix_s_exc_all_l5l3_sx.txt', WeightMatrix_s_exc_all_l5l3_sx)
            np.savetxt('WeightMatrix_s_inh_exc_l5l3_sx.txt', WeightMatrix_s_inh_exc_l5l3_sx)
            np.savetxt('WeightMatrix_s_exc_all_l5l3_dx.txt', WeightMatrix_s_exc_all_l5l3_dx)
            np.savetxt('WeightMatrix_s_inh_exc_l5l3_dx.txt', WeightMatrix_s_inh_exc_l5l3_dx)
            # Layer 5 to Layer 4
            np.savetxt('WeightMatrix_s_exc_all_l5l4_sx.txt', WeightMatrix_s_exc_all_l5l4_sx)
            np.savetxt('WeightMatrix_s_inh_exc_l5l4_sx.txt', WeightMatrix_s_inh_exc_l5l4_sx)
            np.savetxt('WeightMatrix_s_exc_all_l5l4_dx.txt', WeightMatrix_s_exc_all_l5l4_dx)
            np.savetxt('WeightMatrix_s_inh_exc_l5l4_dx.txt', WeightMatrix_s_inh_exc_l5l4_dx)
    os.chdir(SimDir)

    # --- Saving the connection numbers ---
    os.chdir(NumConnDir)
    np.savetxt('Num_conn_l0.txt', num_conn_l0_sx)
    np.savetxt('Num_conn_l0.txt', num_conn_l0_dx)
    if num_layers > 1:
        np.savetxt('Num_conn_l1.txt', num_conn_l1_sx)
        np.savetxt('Num_conn_l1.txt', num_conn_l1_dx)
    if num_layers > 2:
        np.savetxt('Num_conn_l2.txt', num_conn_l2_sx)
        np.savetxt('Num_conn_l2.txt', num_conn_l2_dx)
    if num_layers > 3:
        np.savetxt('Num_conn_l3.txt', num_conn_l3_sx)
        np.savetxt('Num_conn_l3.txt', num_conn_l3_dx)
    if num_layers > 4:
        np.savetxt('Num_conn_l4.txt', num_conn_l4_sx)
        np.savetxt('Num_conn_l4.txt', num_conn_l4_dx)
    if num_layers > 5:
        np.savetxt('Num_conn_l5.txt', num_conn_l5_sx)
        np.savetxt('Num_conn_l5.txt', num_conn_l5_dx)
    os.chdir(SimDir)
