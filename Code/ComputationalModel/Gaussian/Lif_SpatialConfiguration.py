from Lif_Parameters import *
from Lif_Equations import eqs_LIF
seed(val)
#####################################
# Initialize the neuronal populations
#####################################
# ------- Layer 0 -------
l0 = NeuronGroup(n_cells, eqs_LIF, method='exponential_euler', threshold='vm>v_th', reset='vm=v_r', refractory=refractory_period)
l0_exc = l0[:n_exc]
l0_inh = l0[n_exc:]
layers = [l0]
layers_exc = [l0_exc]
layers_inh = [l0_inh]

# ------- Layer 1 -------
if num_layers > 1:
    l1 = NeuronGroup(n_cells, eqs_LIF, method='exponential_euler', threshold='vm>v_th', reset='vm=v_r', refractory=refractory_period)
    l1_exc = l1[:n_exc]
    l1_inh = l1[n_exc:]
    layers.append(l1)
    layers_exc.append(l1_exc)
    layers_inh.append(l1_inh)
    
# ------- Layer 2 -------
if num_layers > 2:
    l2 = NeuronGroup(n_cells, eqs_LIF, method='exponential_euler', threshold='vm>v_th', reset='vm=v_r', refractory=refractory_period)
    l2_exc = l2[:n_exc]
    l2_inh = l2[n_exc:]
    layers.append(l2)
    layers_exc.append(l2_exc)
    layers_inh.append(l2_inh)

# ------- Layer 3 -------
if num_layers > 3:
    l3 = NeuronGroup(n_cells, eqs_LIF, method='exponential_euler', threshold='vm>v_th', reset='vm=v_r', refractory=refractory_period)
    l3_exc = l3[:n_exc]
    l3_inh = l3[n_exc:]
    layers.append(l3)
    layers_exc.append(l3_exc)
    layers_inh.append(l3_inh)
    
# ------- Layer 4 -------
if num_layers > 4:
    l4 = NeuronGroup(n_cells, eqs_LIF, method='exponential_euler', threshold='vm>v_th', reset='vm=v_r', refractory=refractory_period)
    l4_exc = l4[:n_exc]
    l4_inh = l4[n_exc:]
    layers.append(l4)
    layers_exc.append(l4_exc)
    layers_inh.append(l4_inh)
    
# ------- Layer 5 -------
if num_layers > 5:
    l5 = NeuronGroup(n_cells, eqs_LIF, method='exponential_euler', threshold='vm>v_th', reset='vm=v_r', refractory=refractory_period)
    l5_exc = l5[:n_exc]
    l5_inh = l5[n_exc:]
    layers.append(l5)
    layers_exc.append(l5_exc)
    layers_inh.append(l5_inh)
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

    grid_dist_noise = np.random.randint(-(grid_dist / 7) / umeter, (grid_dist / 7) / umeter, size=n_cells) * umeter
    pop_ll.x = x_pos * umeter + grid_dist_noise[pop_ll.i]
    pop_ll.y = y_pos * umeter + grid_dist_noise[pop_ll.i]
    pop_ll.z = pp*dist_pop
    pop_ll.index = pos_index
    pop_ll.index_up = True
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
    ax.set_xlabel('x (??m)')
    ax.set_ylabel('y (??m)')
    ax.set_zlabel('z (??m)')
    axis('auto')
    tight_layout()
    show()

    fig2 = figure(2, figsize=[6*cm_fig_conv, 6*cm_fig_conv], dpi=600.0)
    plot(l0_exc.x / umeter, l0_exc.y / umeter, 'r.', l0_inh.x / umeter, l0_inh.y / umeter, 'b.', markersize='1')
    xlabel('x (??m)')
    ylabel('y (??m)', rotation='vertical')
    gca().spines['top'].set_visible(False)
    gca().spines['right'].set_visible(False)
    axis('equal')
    tight_layout()
    show()

    fig1.savefig('3D layout_3l', dpi=600, transparent=True)
    # fig2.savefig('2D layout_nnn', dpi=600, transparent=True)
    print("Lif_SpatialConfiguration - END")
else:
    print("Finished Lif_SpatialConfiguration ")
####################################
####################################
