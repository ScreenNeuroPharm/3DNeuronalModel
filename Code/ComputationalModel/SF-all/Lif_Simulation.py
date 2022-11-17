from ExternalFunctions import crea_e_o_apri, savelookuptable, savepeaktrain
import os
from Lif_Initialization import *

seed(val)
###########################################
#       User defined simulation parameters
###########################################
#  ---------- Simulation ----------
preptime = 10 * second       # Length of not recorded simulation to settle the network (preparation time) - minimum=1*second
simtime = 60 * second       # Length of the simulation
chunk = 10                  # Length of simulation chunks
savefile = True             # False=do not save output  -   True=save output
DateFolder = '2022-04-24'   # Date in the format aaaa-mm-dd
SimulationNumber = 2        # Number of simulation in the date

comment = '\nNumber of layers:\t' + str(num_layers) + \
          '\nUpward connections:\t' + str(upward_conn) + \
          '\nDownward connections:\t' + str(downward_conn) + \
          '\nspatial config--> index changes 4 each layer' \
          '\nrefractory period = ' + str(refractory_period) + \
          '\nimplementation of scale free connectivity: ' \
          '\n\t - slope: ' + str(a) + \
          '\n\t - lower degree limit: ' + str(lower_degree_limit) + \
          '\n\t - upper degree limit: ' + str(upper_degree_limit) + \
          '\n\t - number of degree values: ' + str(number_degree_values)
#          '\n\nmax conn = ' + str(max_conn_int) + \
#           '\nupward:\tL0 to L1; L0 to L2; L0 to L3; L1 to L2; L1 to L3; L2 to L3' \
#           '\ndownward:\tL1 to L0; L2 to L0; L2 to L1; L3 to L0; L3 to L1; L3 to L2' \

defaultclock.dt = 0.1 * ms      # Clock of the simulation

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
#         Simulation & Online Saving
#####################################
net = Network(collect())
net.run(preptime, report='text')

sm_layer0 = SpikeMonitor(l0, variables='vm')
net.add(sm_layer0)
if num_layers > 1:
    sm_layer1 = SpikeMonitor(l1, variables='vm')
    net.add(sm_layer1)
if num_layers > 2:
    sm_layer2 = SpikeMonitor(l2, variables='vm')
    net.add(sm_layer2)
if num_layers > 3:
    sm_layer3 = SpikeMonitor(l3, variables='vm')
    net.add(sm_layer3)
if num_layers > 4:
    sm_layer4 = SpikeMonitor(l4, variables='vm')
    net.add(sm_layer4)
if num_layers > 5:
    sm_layer5 = SpikeMonitor(l5, variables='vm')
    net.add(sm_layer5)

if savefile:
    quoziente = simtime / second // chunk
    resto = simtime / second % chunk
    times = np.repeat(chunk * second, quoziente).tolist()
    times.append(resto * second)

    crea_e_o_apri('Results_3D_LIF')
    crea_e_o_apri(DateFolder)
    crea_e_o_apri('Simulazione' + str(SimulationNumber))
    SimDir = os.getcwd()

    crea_e_o_apri('PeakTrains')
    for _, tt in enumerate(times):
        sm_l0 = SpikeMonitor(l0, variables='vm')
        net.add(sm_l0)
        if num_layers > 1:
            sm_l1 = SpikeMonitor(l1, variables='vm')
            net.add(sm_l1)
        if num_layers > 2:
            sm_l2 = SpikeMonitor(l2, variables='vm')
            net.add(sm_l2)
        if num_layers > 3:
            sm_l3 = SpikeMonitor(l3, variables='vm')
            net.add(sm_l3)
        if num_layers > 4:
            sm_l4 = SpikeMonitor(l4, variables='vm')
            net.add(sm_l4)
        if num_layers > 5:
            sm_l5 = SpikeMonitor(l5, variables='vm')
            net.add(sm_l5)

        net.run(tt, report='text')

        train_l0 = sm_l0.spike_trains()
        savepeaktrain('ptl0', train_l0, l0, n_exc, n_inh, 0)
        del sm_l0
        if num_layers > 1:
            train_l1 = sm_l1.spike_trains()
            savepeaktrain('ptl1', train_l1, l1, n_exc, n_inh, 1)
            del sm_l1
        if num_layers > 2:
            train_l2 = sm_l2.spike_trains()
            savepeaktrain('ptl2', train_l2, l2, n_exc, n_inh, 2)
            del sm_l2
        if num_layers > 3:
            train_l3 = sm_l3.spike_trains()
            savepeaktrain('ptl3', train_l3, l3, n_exc, n_inh, 3)
            del sm_l3
        if num_layers > 4:
            train_l4 = sm_l4.spike_trains()
            savepeaktrain('ptl4', train_l4, l4, n_exc, n_inh, 4)
            del sm_l4
        if num_layers > 5:
            train_l5 = sm_l5.spike_trains()
            savepeaktrain('ptl5', train_l5, l5, n_exc, n_inh, 5)
            del sm_l5
    os.chdir(SimDir)
else:
    net.run(simtime, report='text')
#####################################


#####################################
#         Results display
#####################################
# ------ RASTER ------
fig1 = figure(1, figsize=[3 * num_layers, 6.5], dpi=600.0)

ax0 = subplot(1, num_layers, 1)
plot(sm_layer0.t / ms / 1000, sm_layer0.i, ',k')
ylabel('Neuron index')
xlim(preptime / second, (preptime + simtime) / second)
xlabel('Time (s)')
if num_layers > 1:
    ax1 = subplot(1, num_layers, 2, sharex=ax0)
    plot(sm_layer1.t / ms / 1000, sm_layer1.i, ',k')
    xlim(preptime / second, (preptime + simtime) / second)
    xlabel('Time (s)')
if num_layers > 2:
    ax2 = subplot(1, num_layers, 3, sharex=ax0)
    plot(sm_layer2.t / ms / 1000, sm_layer2.i, ',k')
    xlim(preptime / second, (preptime + simtime) / second)
    xlabel('Time (s)')
if num_layers > 3:
    ax3 = subplot(1, num_layers, 4, sharex=ax1)
    plot(sm_layer3.t / ms / 1000, sm_layer3.i, ',k')
    xlim(preptime / second, (preptime + simtime) / second)
    xlabel('Time (s)')
if num_layers > 4:
    ax4 = subplot(1, num_layers, 5, sharex=ax1)
    plot(sm_layer4.t / ms / 1000, sm_layer4.i, ',k')
    xlim(preptime / second, (preptime + simtime) / second)
    xlabel('Time (s)')
if num_layers > 5:
    ax5 = subplot(1, num_layers, 6, sharex=ax1)
    plot(sm_layer5.t / ms / 1000, sm_layer5.i, ',k')
    xlim(preptime / second, (preptime + simtime) / second)
    xlabel('Time (s)')

tight_layout()
# show()

# ------ RASTER shuffle ------
fig2 = figure(2, figsize=[10, 6.5], dpi=600.0)

ax0 = subplot(1, num_layers, 1)
pt_l0 = sm_layer0.spike_trains()
for key, __ in pt_l0.items():
    plot(pt_l0[key], repeat(l0.index[key], len(pt_l0[key])), ',b')
ylabel('Neuron index (position)')
xlim(preptime / second, (preptime + simtime) / second)
xlabel('Time (s)')
if num_layers > 1:
    ax1 = subplot(1, num_layers, 2, sharex=ax1)
    pt_l1 = sm_layer1.spike_trains()
    for key, __ in pt_l1.items():
        plot(pt_l1[key], repeat(l1.index[key], len(pt_l1[key])), ',b')
    xlim(preptime / second, (preptime + simtime) / second)
    xlabel('Time (s)')
if num_layers > 2:
    ax2 = subplot(1, num_layers, 3, sharex=ax0)
    pt_l2 = sm_layer2.spike_trains()
    for key, __ in pt_l2.items():
        plot(pt_l2[key], repeat(l2.index[key], len(pt_l2[key])), ',b')
    xlim(preptime / second, (preptime + simtime) / second)
    xlabel('Time (s)')
if num_layers > 3:
    ax3 = subplot(1, num_layers, 4, sharex=ax0)
    pt_l3 = sm_layer3.spike_trains()
    for key, __ in pt_l3.items():
        plot(pt_l3[key], repeat(l3.index[key], len(pt_l3[key])), ',b')
    xlim(preptime / second, (preptime + simtime) / second)
    xlabel('Time (s)')
if num_layers > 4:
    ax4 = subplot(1, num_layers, 5, sharex=ax0)
    pt_l4 = sm_layer4.spike_trains()
    for key, __ in pt_l4.items():
        plot(pt_l4[key], repeat(l4.index[key], len(pt_l4[key])), ',b')
    xlim(preptime / second, (preptime + simtime) / second)
    xlabel('Time (s)')
if num_layers > 5:
    ax5 = subplot(1, num_layers, 6, sharex=ax0)
    pt_l5 = sm_layer5.spike_trains()
    for key, __ in pt_l5.items():
        plot(pt_l5[key], repeat(l5.index[key], len(pt_l5[key])), ',b')
    xlim(preptime / second, (preptime + simtime) / second)
    xlabel('Time (s)')

tight_layout()
# show()

# ------ 3D layout ------
fig3 = figure(3, figsize=[15, 10], dpi=600.0)
ax = axes(projection='3d')
for layer in layers_exc:
    ax.scatter(layer.x / umeter, layer.y / umeter, layer.z / umeter, c='r', marker='.')
for layer in layers_inh:
    ax.scatter(layer.x / umeter, layer.y / umeter, layer.z / umeter, c='b', marker='.')
ax.set_xlabel('X (µm)')
ax.set_ylabel('Y (µm)')
ax.set_zlabel('Z (µm)')
axis('auto')

tight_layout()
# show()

# ------ VALUES ------
mfr_l0 = mean(sm_layer0.count / simtime)
print('MFR_l0: ' + str(mfr_l0))
if num_layers > 1:
    mfr_l1 = mean(sm_layer1.count / simtime)
    print('MFR_l1: ' + str(mfr_l1))
if num_layers > 2:
    mfr_l2 = mean(sm_layer2.count / simtime)
    print('MFR_l2: ' + str(mfr_l2))
if num_layers > 3:
    mfr_l3 = mean(sm_layer3.count / simtime)
    print('MFR_l3: ' + str(mfr_l3))
if num_layers > 4:
    mfr_l4 = mean(sm_layer4.count / simtime)
    print('MFR_l4: ' + str(mfr_l4))
if num_layers > 5:
    mfr_l5 = mean(sm_layer5.count / simtime)
    print('MFR_l5: ' + str(mfr_l5))
#####################################


#####################################
#         Saving Results
#####################################
if savefile:
    os.chdir(SimDir)

    # --- Saving the specifics of the simulation: duration, frequency, input equation, MFR, etc. ---
    with open("SimulationSpec", "a+") as f:
        f.write(DateFolder + '\tSimulation ' + str(SimulationNumber) + '\n')
        f.write('\nDuration: \t%s' % simtime)
        f.write('\nSampling frequency: \t%s' % (1 / defaultclock.dt))
        f.write('\nSeed: \t%s' % val)
        f.write('\nNumber of:'
                '\n\t- Cells per layer \t%s'
                '\n\t- Exc cells \t%s'
                '\n\t- Inh cells \t%s'
                % (n_cells, n_exc, n_inh))
        f.write('\n\nInput value: ' + input_value)
        f.write('\n\nMFR_l0: %s' % mfr_l0)
        if num_layers > 1:
            f.write('\nMFR_l1: %s' % mfr_l1)
        if num_layers > 2:
            f.write('\nMFR_l2: %s' % mfr_l2)
        if num_layers > 3:
            f.write('\nMFR_l3: %s' % mfr_l3)
        if num_layers > 4:
            f.write('\nMFR_l4: %s' % mfr_l4)
        if num_layers > 5:
            f.write('\nMFR_l5: %s' % mfr_l5)
        f.write('\n\n' + comment)
        f.write('\n\n\n------------------------------\n\n\n')
        f.close()

    # --- Saving the lookup table ---
    crea_e_o_apri('LookupTable')
    savelookuptable('lookuptable_l0', train_l0, l0)
    if num_layers > 1:
        savelookuptable('lookuptable_l1', train_l1, l1)
    if num_layers > 2:
        savelookuptable('lookuptable_l2', train_l2, l2)
    if num_layers > 3:
        savelookuptable('lookuptable_l3', train_l3, l3)
    if num_layers > 4:
        savelookuptable('lookuptable_l4', train_l4, l4)
    if num_layers > 5:
        savelookuptable('lookuptable_l5', train_l5, l5)
    os.chdir(SimDir)

    # --- Saving the weight matrices ---
    crea_e_o_apri('WeightMatrix')
    # Layer 0
    np.savetxt('WeightMatrix_s_exc_all_l0.txt', WeightMatrix_s_exc_all_l0)
    np.savetxt('WeightMatrix_s_inh_exc_l0.txt', WeightMatrix_s_inh_exc_l0)
    # Layer 1
    if num_layers > 1:
        np.savetxt('WeightMatrix_s_exc_all_l1.txt', WeightMatrix_s_exc_all_l1)
        np.savetxt('WeightMatrix_s_inh_exc_l1.txt', WeightMatrix_s_inh_exc_l1)
        if upward_conn:
            # Layer 0 to Layer 1
            np.savetxt('WeightMatrix_s_exc_all_l0l1.txt', WeightMatrix_s_exc_all_l0l1)
            np.savetxt('WeightMatrix_s_inh_exc_l0l1.txt', WeightMatrix_s_inh_exc_l0l1)
        if downward_conn:
            # Layer 1 to Layer 0
            np.savetxt('WeightMatrix_s_exc_all_l1l0.txt', WeightMatrix_s_exc_all_l1l0)
            np.savetxt('WeightMatrix_s_inh_exc_l1l0.txt', WeightMatrix_s_inh_exc_l1l0)
    # Layer 2
    if num_layers > 2:
        np.savetxt('WeightMatrix_s_exc_all_l2.txt', WeightMatrix_s_exc_all_l2)
        np.savetxt('WeightMatrix_s_inh_exc_l2.txt', WeightMatrix_s_inh_exc_l2)
        if upward_conn:
            # Layer 0 to Layer 2
            np.savetxt('WeightMatrix_s_exc_all_l0l2.txt', WeightMatrix_s_exc_all_l0l2)
            np.savetxt('WeightMatrix_s_inh_exc_l0l2.txt', WeightMatrix_s_inh_exc_l0l2)
            # Layer 1 to Layer 2
            np.savetxt('WeightMatrix_s_exc_all_l1l2.txt', WeightMatrix_s_exc_all_l1l2)
            np.savetxt('WeightMatrix_s_inh_exc_l1l2.txt', WeightMatrix_s_inh_exc_l1l2)
        if downward_conn:
            # Layer 2 to Layer 0
            np.savetxt('WeightMatrix_s_exc_all_l2l0.txt', WeightMatrix_s_exc_all_l2l0)
            np.savetxt('WeightMatrix_s_inh_exc_l2l0.txt', WeightMatrix_s_inh_exc_l2l0)
            # Layer 2 to Layer 1
            np.savetxt('WeightMatrix_s_exc_all_l2l1.txt', WeightMatrix_s_exc_all_l2l1)
            np.savetxt('WeightMatrix_s_inh_exc_l2l1.txt', WeightMatrix_s_inh_exc_l2l1)
    # Layer 3
    if num_layers > 3:
        np.savetxt('WeightMatrix_s_exc_all_l3.txt', WeightMatrix_s_exc_all_l3)
        np.savetxt('WeightMatrix_s_inh_exc_l3.txt', WeightMatrix_s_inh_exc_l3)
        if upward_conn:
            # Layer 0 to Layer 3
            np.savetxt('WeightMatrix_s_exc_all_l0l3.txt', WeightMatrix_s_exc_all_l0l3)
            np.savetxt('WeightMatrix_s_inh_exc_l0l3.txt', WeightMatrix_s_inh_exc_l0l3)
            # Layer 1 to Layer 3
            np.savetxt('WeightMatrix_s_exc_all_l1l3.txt', WeightMatrix_s_exc_all_l1l3)
            np.savetxt('WeightMatrix_s_inh_exc_l1l3.txt', WeightMatrix_s_inh_exc_l1l3)
            # Layer 2 to Layer 3
            np.savetxt('WeightMatrix_s_exc_all_l2l3.txt', WeightMatrix_s_exc_all_l2l3)
            np.savetxt('WeightMatrix_s_inh_exc_l2l3.txt', WeightMatrix_s_inh_exc_l2l3)
        if downward_conn:
            # Layer 3 to Layer 0
            np.savetxt('WeightMatrix_s_exc_all_l3l0.txt', WeightMatrix_s_exc_all_l3l0)
            np.savetxt('WeightMatrix_s_inh_exc_l3l0.txt', WeightMatrix_s_inh_exc_l3l0)
            # Layer 3 to Layer 1
            np.savetxt('WeightMatrix_s_exc_all_l3l1.txt', WeightMatrix_s_exc_all_l3l1)
            np.savetxt('WeightMatrix_s_inh_exc_l3l1.txt', WeightMatrix_s_inh_exc_l3l1)
            # Layer 3 to Layer 2
            np.savetxt('WeightMatrix_s_exc_all_l3l2.txt', WeightMatrix_s_exc_all_l3l2)
            np.savetxt('WeightMatrix_s_inh_exc_l3l2.txt', WeightMatrix_s_inh_exc_l3l2)
    # Layer 4
    if num_layers > 4:
        np.savetxt('WeightMatrix_s_exc_all_l4.txt', WeightMatrix_s_exc_all_l4)
        np.savetxt('WeightMatrix_s_inh_exc_l4.txt', WeightMatrix_s_inh_exc_l4)
        if upward_conn:
            # Layer 0 to Layer 4
            np.savetxt('WeightMatrix_s_exc_all_l0l4.txt', WeightMatrix_s_exc_all_l0l4)
            np.savetxt('WeightMatrix_s_inh_exc_l0l4.txt', WeightMatrix_s_inh_exc_l0l4)
            # Layer 1 to Layer 4
            np.savetxt('WeightMatrix_s_exc_all_l1l4.txt', WeightMatrix_s_exc_all_l1l4)
            np.savetxt('WeightMatrix_s_inh_exc_l1l4.txt', WeightMatrix_s_inh_exc_l1l4)
            # Layer 2 to Layer 4
            np.savetxt('WeightMatrix_s_exc_all_l2l4.txt', WeightMatrix_s_exc_all_l2l4)
            np.savetxt('WeightMatrix_s_inh_exc_l2l4.txt', WeightMatrix_s_inh_exc_l2l4)
            # Layer 3 to Layer 4
            np.savetxt('WeightMatrix_s_exc_all_l3l4.txt', WeightMatrix_s_exc_all_l3l4)
            np.savetxt('WeightMatrix_s_inh_exc_l3l4.txt', WeightMatrix_s_inh_exc_l3l4)
        if downward_conn:
            # Layer 4 to Layer 0
            np.savetxt('WeightMatrix_s_exc_all_l4l0.txt', WeightMatrix_s_exc_all_l4l0)
            np.savetxt('WeightMatrix_s_inh_exc_l4l0.txt', WeightMatrix_s_inh_exc_l4l0)
            # Layer 4 to Layer 1
            np.savetxt('WeightMatrix_s_exc_all_l4l1.txt', WeightMatrix_s_exc_all_l4l1)
            np.savetxt('WeightMatrix_s_inh_exc_l4l1.txt', WeightMatrix_s_inh_exc_l4l1)
            # Layer 4 to Layer 2
            np.savetxt('WeightMatrix_s_exc_all_l4l2.txt', WeightMatrix_s_exc_all_l4l2)
            np.savetxt('WeightMatrix_s_inh_exc_l4l2.txt', WeightMatrix_s_inh_exc_l4l2)
            # Layer 4 to Layer 3
            np.savetxt('WeightMatrix_s_exc_all_l4l3.txt', WeightMatrix_s_exc_all_l4l3)
            np.savetxt('WeightMatrix_s_inh_exc_l4l3.txt', WeightMatrix_s_inh_exc_l4l3)
    # Layer 5
    if num_layers > 5:
        np.savetxt('WeightMatrix_s_exc_all_l5.txt', WeightMatrix_s_exc_all_l5)
        np.savetxt('WeightMatrix_s_inh_exc_l5.txt', WeightMatrix_s_inh_exc_l5)
        if upward_conn:
            # Layer 0 to Layer 5
            np.savetxt('WeightMatrix_s_exc_all_l0l5.txt', WeightMatrix_s_exc_all_l0l5)
            np.savetxt('WeightMatrix_s_inh_exc_l0l5.txt', WeightMatrix_s_inh_exc_l0l5)
            # Layer 1 to Layer 5
            np.savetxt('WeightMatrix_s_exc_all_l1l5.txt', WeightMatrix_s_exc_all_l1l5)
            np.savetxt('WeightMatrix_s_inh_exc_l1l5.txt', WeightMatrix_s_inh_exc_l1l5)
            # Layer 2 to Layer 5
            np.savetxt('WeightMatrix_s_exc_all_l2l5.txt', WeightMatrix_s_exc_all_l2l5)
            np.savetxt('WeightMatrix_s_inh_exc_l2l5.txt', WeightMatrix_s_inh_exc_l2l5)
            # Layer 3 to Layer 5
            np.savetxt('WeightMatrix_s_exc_all_l3l5.txt', WeightMatrix_s_exc_all_l3l5)
            np.savetxt('WeightMatrix_s_inh_exc_l3l5.txt', WeightMatrix_s_inh_exc_l3l5)
            # Layer 4 to Layer 5
            np.savetxt('WeightMatrix_s_exc_all_l4l5.txt', WeightMatrix_s_exc_all_l4l5)
            np.savetxt('WeightMatrix_s_inh_exc_l4l5.txt', WeightMatrix_s_inh_exc_l4l5)
        if downward_conn:
            # Layer 5 to Layer 0
            np.savetxt('WeightMatrix_s_exc_all_l5l0.txt', WeightMatrix_s_exc_all_l5l0)
            np.savetxt('WeightMatrix_s_inh_exc_l5l0.txt', WeightMatrix_s_inh_exc_l5l0)
            # Layer 5 to Layer 1
            np.savetxt('WeightMatrix_s_exc_all_l5l1.txt', WeightMatrix_s_exc_all_l5l1)
            np.savetxt('WeightMatrix_s_inh_exc_l5l1.txt', WeightMatrix_s_inh_exc_l5l1)
            # Layer 5 to Layer 2
            np.savetxt('WeightMatrix_s_exc_all_l5l2.txt', WeightMatrix_s_exc_all_l5l2)
            np.savetxt('WeightMatrix_s_inh_exc_l5l2.txt', WeightMatrix_s_inh_exc_l5l2)
            # Layer 5 to Layer 3
            np.savetxt('WeightMatrix_s_exc_all_l5l3.txt', WeightMatrix_s_exc_all_l5l3)
            np.savetxt('WeightMatrix_s_inh_exc_l5l3.txt', WeightMatrix_s_inh_exc_l5l3)
            # Layer 5 to Layer 4
            np.savetxt('WeightMatrix_s_exc_all_l5l4.txt', WeightMatrix_s_exc_all_l5l4)
            np.savetxt('WeightMatrix_s_inh_exc_l5l4.txt', WeightMatrix_s_inh_exc_l5l4)
    os.chdir(SimDir)

    # --- Saving the connection numbers ---
    crea_e_o_apri('NumConnections')
    np.savetxt('Num_conn_l0.txt', num_conn_l0)
    if num_layers > 1:
        np.savetxt('Num_conn_l1.txt', num_conn_l1)
    if num_layers > 2:
        np.savetxt('Num_conn_l2.txt', num_conn_l2)
    if num_layers > 3:
        np.savetxt('Num_conn_l3.txt', num_conn_l3)
    if num_layers > 4:
        np.savetxt('Num_conn_l4.txt', num_conn_l4)
    if num_layers > 5:
        np.savetxt('Num_conn_l5.txt', num_conn_l5)
    os.chdir(SimDir)

    # --- Saving the images ---
    if inputtype == 1:
        fname1 = "1HP_noisySTIM_" + str(int(simtime / second)) + "s_SIM" + str(SimulationNumber) + ".tif"
        fname2 = "1HP_noisySTIM_" + str(int(simtime / second)) + "s_SIM" + str(SimulationNumber) + "_shuffle.tif"
    elif inputtype == 0:
        fname1 = "1HP_stepSTIM_SIM" + str(SimulationNumber) + ".tif"
        fname2 = "1HP_stepSTIM_SIM" + str(SimulationNumber) + "_shuffle.tif"
    fig1.savefig(fname1, dpi=600, transparent=True)
    fig2.savefig(fname2, dpi=600, transparent=True)
    fig3.savefig('3D layout', dpi=600, transparent=True)
