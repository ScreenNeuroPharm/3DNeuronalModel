from Lif_Parameters import *
from Lif_Equations import eqs_LIF
seed(val)
#####################################
# Initialize the neuronal populations - SX
#####################################
# ------- Layer 0 -------
# l0 = NeuronGroup(2*n_cells, eqs_LIF, method='exponential_euler', threshold='vm>v_th', reset='vm=v_r', refractory=refractory_period)
l0_sx = NeuronGroup(n_cells, eqs_LIF, method='exponential_euler', threshold='vm>v_th', reset='vm=v_r', refractory=refractory_period)
l0_dx = NeuronGroup(n_cells, eqs_LIF, method='exponential_euler', threshold='vm>v_th', reset='vm=v_r', refractory=refractory_period)
# l0_sx = l0[:n_cells]
l0_sx_exc = l0_sx[:n_exc]
l0_sx_inh = l0_sx[n_exc:]
# l0_dx = l0[n_cells:]
l0_dx_exc = l0_dx[:n_exc]
l0_dx_inh = l0_dx[n_exc:]

# layers = [l0]
layers = [l0_sx, l0_dx]       # per ciclo automatico in connection con zip
layers_sx = [l0_sx]
layers_dx = [l0_dx]
layers_exc = [l0_sx_exc, l0_dx_exc]
layers_inh = [l0_sx_inh, l0_dx_inh]

# ------- Layer 1 -------
if num_layers > 1:
    l1_sx = NeuronGroup(n_cells, eqs_LIF, method='exponential_euler', threshold='vm>v_th', reset='vm=v_r', refractory=refractory_period)
    l1_dx = NeuronGroup(n_cells, eqs_LIF, method='exponential_euler', threshold='vm>v_th', reset='vm=v_r', refractory=refractory_period)

    # l1_sx = l1[:n_cells]
    l1_sx_exc = l1_sx[:n_exc]
    l1_sx_inh = l1_sx[n_exc:]
    # l1_dx = l1[n_cells:]
    l1_dx_exc = l1_dx[:n_exc]
    l1_dx_inh = l1_dx[n_exc:]

    # layers.append(l1)
    layers.extend([l1_sx, l1_dx])
    layers_sx.append(l1_sx)
    layers_dx.append(l1_dx)
    layers_exc.extend([l1_sx_exc, l1_dx_exc])
    layers_inh.extend([l1_sx_inh, l1_dx_inh])
    
# ------- Layer 2 -------
if num_layers > 2:
    l2_sx = NeuronGroup(n_cells, eqs_LIF, method='exponential_euler', threshold='vm>v_th', reset='vm=v_r', refractory=refractory_period)
    l2_dx = NeuronGroup(n_cells, eqs_LIF, method='exponential_euler', threshold='vm>v_th', reset='vm=v_r', refractory=refractory_period)
    
    # l2_sx = l2[:n_cells]
    l2_sx_exc = l2_sx[:n_exc]
    l2_sx_inh = l2_sx[n_exc:]
    # l2_dx = l2[n_cells:]
    l2_dx_exc = l2_dx[:n_exc]
    l2_dx_inh = l2_dx[n_exc:]

    # layers.append(l2)
    layers.extend([l2_sx, l2_dx])
    layers_sx.append(l2_sx)
    layers_dx.append(l2_dx)
    layers_exc.extend([l2_sx_exc, l2_dx_exc])
    layers_inh.extend([l2_sx_inh, l2_dx_inh])

# ------- Layer 3 -------
if num_layers > 3:
    # l3 = NeuronGroup(2*n_cells, eqs_LIF, method='exponential_euler', threshold='vm>v_th', reset='vm=v_r', refractory=refractory_period)
    l3_sx = NeuronGroup(n_cells, eqs_LIF, method='exponential_euler', threshold='vm>v_th', reset='vm=v_r', refractory=refractory_period)
    l3_dx = NeuronGroup(n_cells, eqs_LIF, method='exponential_euler', threshold='vm>v_th', reset='vm=v_r', refractory=refractory_period)
    
    # l3_sx = l3[:n_cells]
    l3_sx_exc = l3_sx[:n_exc]
    l3_sx_inh = l3_sx[n_exc:]
    # l3_dx = l3[n_cells:]
    l3_dx_exc = l3_dx[:n_exc]
    l3_dx_inh = l3_dx[n_exc:]

    # layers.append(l3)
    layers.extend([l3_sx, l3_dx])
    layers_sx.append(l3_sx)
    layers_dx.append(l3_dx)
    layers_exc.extend([l3_sx_exc, l3_dx_exc])
    layers_inh.extend([l3_sx_inh, l3_dx_inh])
    
# ------- Layer 4 -------
if num_layers > 4:
    # l4 = NeuronGroup(2*n_cells, eqs_LIF, method='exponential_euler', threshold='vm>v_th', reset='vm=v_r', refractory=refractory_period)
    l4_sx = NeuronGroup(n_cells, eqs_LIF, method='exponential_euler', threshold='vm>v_th', reset='vm=v_r', refractory=refractory_period)
    l4_dx = NeuronGroup(n_cells, eqs_LIF, method='exponential_euler', threshold='vm>v_th', reset='vm=v_r', refractory=refractory_period)

    # l4_sx = l4[:n_cells]
    l4_sx_exc = l4_sx[:n_exc]
    l4_sx_inh = l4_sx[n_exc:]
    # l4_dx = l4[n_cells:]
    l4_dx_exc = l4_dx[:n_exc]
    l4_dx_inh = l4_dx[n_exc:]

    # layers.append(l4)
    layers.extend([l4_sx, l4_dx])
    layers_sx.append(l4_sx)
    layers_dx.append(l4_dx)
    layers_exc.extend([l4_sx_exc, l4_dx_exc])
    layers_inh.extend([l4_sx_inh, l4_dx_inh])
    
# ------- Layer 5 -------
if num_layers > 5:
    # l5 = NeuronGroup(2*n_cells, eqs_LIF, method='exponential_euler', threshold='vm>v_th', reset='vm=v_r', refractory=refractory_period)
    l5_sx = NeuronGroup(n_cells, eqs_LIF, method='exponential_euler', threshold='vm>v_th', reset='vm=v_r', refractory=refractory_period)
    l5_dx = NeuronGroup(n_cells, eqs_LIF, method='exponential_euler', threshold='vm>v_th', reset='vm=v_r', refractory=refractory_period)
    
    # l5_sx = l5[:n_cells]
    l5_sx_exc = l5_sx[:n_exc]
    l5_sx_inh = l5_sx[n_exc:]
    # l5_dx = l5[n_cells:]
    l5_dx_exc = l5_dx[:n_exc]
    l5_dx_inh = l5_dx[n_exc:]

    # layers.append(l5)
    layers.extend([l5_sx, l5_dx])
    layers_sx.append(l5_sx)
    layers_dx.append(l5_dx)
    layers_exc.extend([l5_sx_exc, l5_dx_exc])
    layers_inh.extend([l5_sx_inh, l5_dx_inh])
#####################################


#####################################
# Initialize the grid positions
#####################################
for pp, pop_ll in enumerate(layers):
    pos_index = arange(n_cells)
    np.random.shuffle(pos_index)

    x_pos = np.zeros_like(pos_index)
    y_pos = np.zeros_like(pos_index)

    for index in np.arange(n_cells):
        x_pos[index] = (pos_index[index] // rows) * (grid_dist / umeter) - rows / 2.0 * (grid_dist / umeter)
        y_pos[index] = (pos_index[index] % rows) * (grid_dist / umeter) - cols / 2.0 * (grid_dist / umeter)

    grid_dist_noise = np.random.randint(-(grid_dist / 7) / umeter, (grid_dist / 7) / umeter, size=2*n_cells) * umeter
    # pop_ll[:n_cells].index = pos_index
    # pop_ll[:n_cells].x = x_pos * umeter + grid_dist_noise[pop_ll[:n_cells].i] - x_trasl
    # pop_ll[:n_cells].y = y_pos * umeter + grid_dist_noise[pop_ll[:n_cells].i]
    pop_ll.index = pos_index
    pop_ll.x = x_pos * umeter + grid_dist_noise[pop_ll.i] - x_trasl * (1 - (pp % 2)) + x_trasl * (pp % 2)
    pop_ll.y = y_pos * umeter + grid_dist_noise[pop_ll.i]
    # pop_ll[n_cells:].index = pos_index
    # pop_ll[n_cells:].x = x_pos * umeter + grid_dist_noise[pop_ll[n_cells:].i] + x_trasl
    # pop_ll[n_cells:].y = y_pos * umeter + grid_dist_noise[pop_ll[n_cells:].i]
    pop_ll.z = (pp//2)*dist_pop
    pop_ll.index_up = True

# for pop in layers_sx:
#     pop.x = pop.x - x_trasl
#
# for pop in layers_dx:
#     pop.x = pop.x + x_trasl
#####################################


####################################
####################################
if __name__ == '__main__':
    fig1 = figure(figsize=[10*cm_fig_conv, 8*cm_fig_conv], dpi=600.0)
    ax = axes(projection='3d')
    for layer in layers_exc:
        ax.scatter(layer.x / umeter, layer.y / umeter, layer.z / umeter, c='r', marker='.', s=3)
    for layer in layers_inh:
        ax.scatter(layer.x / umeter, layer.y / umeter, layer.z / umeter, c='b', marker='.', s=3)
    ax.set_xlabel('x (µm)')
    ax.set_ylabel('y (µm)')
    ax.set_zlabel('z (µm)')
    axis('auto')
    tight_layout()
    show()

    fig2 = figure(2, figsize=[6*cm_fig_conv, 6*cm_fig_conv], dpi=600.0)
    # plot(l0_sx_exc.x / umeter, l0_sx_exc.y / umeter, 'r.', l0_sx_inh.x / umeter, l0_sx_inh.y / umeter, 'b.', markersize='1')
    plot(l0_sx_exc.x / umeter, l0_sx_exc.y / umeter, 'r.', markersize='1')
    # plot(l0_dx_exc.x / umeter, l0_dx_exc.y / umeter, 'r.', l0_dx_inh.x / umeter, l0_dx_inh.y / umeter, 'b.', markersize='1')
    plot(l0_dx_inh.x / umeter, l0_dx_inh.y / umeter, 'b.', markersize='1')
    xlabel('x (µm)')
    ylabel('y (µm)', rotation='vertical')
    gca().spines['top'].set_visible(False)
    gca().spines['right'].set_visible(False)
    axis('equal')
    tight_layout()
    show()

    # fig1.savefig('3D layout_3l', dpi=600, transparent=True)
    # fig2.savefig('2D layout_nnn', dpi=600, transparent=True)
    print("Lif_SpatialConfiguration - END")
else:
    print("Finished Lif_SpatialConfiguration ")
####################################
####################################
