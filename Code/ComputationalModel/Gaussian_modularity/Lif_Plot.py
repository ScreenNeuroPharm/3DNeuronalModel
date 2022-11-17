from Lif_Parameters import *


def plot_raster(spike_monitor, preptime, simtime, filename, s=False, save=False):
    raster = figure(figsize=[5 * cm_fig_conv, 6 * cm_fig_conv], dpi=600.0)
    plot(spike_monitor.t / ms / 1000, spike_monitor.i, ",b")
    xlabel('Time (s)')
    xlim(preptime / second, (preptime + simtime) / second)
    ylabel('Neuron index (position)', rotation='vertical')
    gca().spines['top'].set_visible(False)
    gca().spines['right'].set_visible(False)
    tight_layout()
    if s:
        show()
    if save:
        raster.savefig(filename, dpi=600, transparent=True)


# fig1 = figure(1, figsize=[5 * num_layers * cm_fig_conv, 6 * cm_fig_conv], dpi=600.0)
# fig1 = figure(1, figsize=[3.3 * cm_fig_conv, 6.6 * cm_fig_conv], dpi=600.0)
# ax0 = subplot(1, num_layers, 1)
# plot(sm_layer0.t / ms / 1000, sm_layer0.i, ',k')
# ylabel('Neuron index')
# xlim(preptime / second, (preptime + simtime) / second)
# xlabel('Time (s)')
# if num_layers > 1:
#     ax1 = subplot(1, num_layers, 2, sharex=ax0)
#     plot(sm_layer1.t / ms / 1000, sm_layer1.i, ',k')
#     xlim(preptime / second, (preptime + simtime) / second)
#     xlabel('Time (s)')
# if num_layers > 2:
#     ax2 = subplot(1, num_layers, 3, sharex=ax0)
#     plot(sm_layer2.t / ms / 1000, sm_layer2.i, ',k')
#     xlim(preptime / second, (preptime + simtime) / second)
#     xlabel('Time (s)')
# if num_layers > 3:
#     ax3 = subplot(1, num_layers, 4, sharex=ax1)
#     plot(sm_layer3.t / ms / 1000, sm_layer3.i, ',k')
#     xlim(preptime / second, (preptime + simtime) / second)
#     xlabel('Time (s)')
# if num_layers > 4:
#     ax4 = subplot(1, num_layers, 5, sharex=ax1)
#     plot(sm_layer4.t / ms / 1000, sm_layer4.i, ',k')
#     xlim(preptime / second, (preptime + simtime) / second)
#     xlabel('Time (s)')
# if num_layers > 5:
#     ax5 = subplot(1, num_layers, 6, sharex=ax1)
#     plot(sm_layer5.t / ms / 1000, sm_layer5.i, ',k')
#     xlim(preptime / second, (preptime + simtime) / second)
#     xlabel('Time (s)')
# tight_layout()
# # show()


def plot_raster_shuffle(layer, spike_monitor, preptime, simtime, filename, s=False, save=False):
    raster_shuffle = figure(figsize=[5 * cm_fig_conv, 6 * cm_fig_conv], dpi=600.0)
    pt = spike_monitor.spike_trains()
    for key, __ in pt.items():
        plot(pt[key], repeat(layer.index[key], len(pt[key])), ',k')
    xlabel('Time (s)')
    xlim(preptime / second, (preptime + simtime) / second)
    ylabel('Neuron index (position)', rotation='vertical')
    gca().spines['top'].set_visible(False)
    gca().spines['right'].set_visible(False)
    tight_layout()
    if s:
        show()
    if save:
        raster_shuffle.savefig(filename, dpi=600, transparent=True)


def plot_spatial_configuration(layers_exc, layers_inh, s=False, save=False):
    spatial_config = figure(figsize=[10 * cm_fig_conv, 8 * cm_fig_conv], dpi=600.0)
    ax = axes(projection='3d')
    for layer in layers_exc:
        ax.scatter(layer.x / umeter, layer.y / umeter, layer.z / umeter, c='r', marker='.', s=3)
    for layer in layers_inh:
        ax.scatter(layer.x / umeter, layer.y / umeter, layer.z / umeter, c='b', marker='.', s=3)
    ax.set_xlabel('X (µm)')
    ax.set_ylabel('Y (µm)')
    ax.set_zlabel('Z (µm)')
    axis('auto')
    tight_layout()
    if s:
        show()
    if save:
        spatial_config.savefig('3D layout.tif', dpi=600, transparent=True)

# ------ 3D layout ------
# fig3 = figure(3, figsize=[10 * cm_fig_conv, 8 * cm_fig_conv], dpi=600.0)
# ax = axes(projection='3d')
# for layer in layers_exc:
#     ax.scatter(layer.x / umeter, layer.y / umeter, layer.z / umeter, c='r', marker='.', s=3)
# for layer in layers_inh:
#     ax.scatter(layer.x / umeter, layer.y / umeter, layer.z / umeter, c='b', marker='.', s=3)
# ax.set_xlabel('X (µm)')
# ax.set_ylabel('Y (µm)')
# ax.set_zlabel('Z (µm)')
# axis('auto')
# tight_layout()
# # show()

# if should_save:
    # --- Saving the images ---
    # if inputtype == 1:
    #     fname1 = "1HP_noisySTIM_" + str(int(simtime / second)) + "s_SIM" + str(SimulationNumber) + ".tif"
    #     fname2 = "1HP_noisySTIM_" + str(int(simtime / second)) + "s_SIM" + str(SimulationNumber) + "_shuffle.tif"
    # elif inputtype == 0:
    #     fname1 = "1HP_stepSTIM_SIM" + str(SimulationNumber) + ".tif"
    #     fname2 = "1HP_stepSTIM_SIM" + str(SimulationNumber) + "_shuffle.tif"
    # fig1.savefig(fname1, dpi=600, transparent=True)
    # fig2.savefig(fname2, dpi=600, transparent=True)
    # fig3.savefig('3D layout.tif', dpi=600, transparent=True)
