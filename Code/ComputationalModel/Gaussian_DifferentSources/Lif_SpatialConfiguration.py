from Lif_Parameters import *
from Lif_Equations import eqs_LIF
from math import ceil
import random

seed(val)

#####################################
# Initialize the neuronal populations
#####################################
# ------- Layer 0 -------
l0 = NeuronGroup(n_cells, eqs_LIF, method='exponential_euler', threshold='vm>v_th', reset='vm=v_r',
                 refractory=refractory_period)
l0_exc = l0[:n_exc]
l0_inh = l0[n_exc:]
layers = [l0]
layers_exc = [l0_exc]
layers_inh = [l0_inh]

# ------- Layer 1 -------
if num_layers > 1:
    l1 = NeuronGroup(n_cells, eqs_LIF, method='exponential_euler', threshold='vm>v_th', reset='vm=v_r',
                     refractory=refractory_period)
    l1_exc = l1[:n_exc]
    l1_inh = l1[n_exc:]
    layers.append(l1)
    layers_exc.append(l1_exc)
    layers_inh.append(l1_inh)

# ------- Layer 2 -------
if num_layers > 2:
    l2 = NeuronGroup(n_cells, eqs_LIF, method='exponential_euler', threshold='vm>v_th', reset='vm=v_r',
                     refractory=refractory_period)
    l2_exc = l2[:n_exc]
    l2_inh = l2[n_exc:]
    layers.append(l2)
    layers_exc.append(l2_exc)
    layers_inh.append(l2_inh)

# ------- Layer 3 -------
if num_layers > 3:
    l3 = NeuronGroup(n_cells, eqs_LIF, method='exponential_euler', threshold='vm>v_th', reset='vm=v_r',
                     refractory=refractory_period)
    l3_exc = l3[:n_exc]
    l3_inh = l3[n_exc:]
    layers.append(l3)
    layers_exc.append(l3_exc)
    layers_inh.append(l3_inh)

# ------- Layer 4 -------
if num_layers > 4:
    l4 = NeuronGroup(n_cells, eqs_LIF, method='exponential_euler', threshold='vm>v_th', reset='vm=v_r',
                     refractory=refractory_period)
    l4_exc = l4[:n_exc]
    l4_inh = l4[n_exc:]
    layers.append(l4)
    layers_exc.append(l4_exc)
    layers_inh.append(l4_inh)

# ------- Layer 5 -------
if num_layers > 5:
    l5 = NeuronGroup(n_cells, eqs_LIF, method='exponential_euler', threshold='vm>v_th', reset='vm=v_r',
                     refractory=refractory_period)
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
    pop_ll.z = pp * dist_pop
    pop_ll.index = pos_index
    pop_ll.index_up = False
#####################################

#####################################
# Initialize 3D connections
#####################################
x_coord = np.arange(-429, 429 + 1, (grid_dist / umeter))
# y_coord = np.arange(-442, 390+1, (grid_dist / umeter))

indices_up_layers = []
for ll in layers:  # può essere trasformato in funzione --> così da poter scegliere una casistica diversa per ogni layer
    if type_connectivity == 'random':
        n_up = 359
        r_centre = 0 * umeter
        n_centre = int(n_up)
        indices_up = random.choices(arange(n_cells), k=n_centre)
    elif type_connectivity == 'central':
        n_up = 0  # parametrizzato su raggio
        r_centre = 278 * umeter
        n_centre = 1
        indices_up = []
        x_c = x_coord[rows // 2]
        index_y_c = ceil(rows / 2)
        index_centre = ((x_c / (grid_dist / umeter) + rows / 2) * rows) + index_y_c
        i_centre = np.where(ll.index == index_centre)
        x_c = ll[i_centre].x / umeter
        y_c = ll[i_centre].y / umeter
        for kk in arange(n_cells):
            x_kk = ll[kk].x / umeter
            y_kk = ll[kk].y / umeter
            if (x_kk - x_c) ** 2 + (y_kk - y_c) ** 2 <= (r_centre / umeter) ** 2:
                ll[kk].index_up = True
                indices_up.append(kk)
    elif type_connectivity == 'multiple_centres':
        n_up = 0  # paramtrizzato su raggio
        r_centre = 129 * umeter
        n_centre = 5
        indices_up = []
        x_c = np.asarray(list(x_coord[ii] for ii in
                              [ceil(1 * rows / 5), ceil(1 * rows / 5), ceil(rows / 2), ceil(4 * rows / 5),
                               ceil(4 * rows / 5)]))
        index_y_c = [ceil(1 * rows / 5), ceil(4 * rows / 5), ceil(rows / 2), ceil(1 * rows / 5), ceil(4 * rows / 5)]
        index_centre = ((x_c / (grid_dist / umeter) + rows / 2) * rows) + index_y_c
        i_centre = []
        x_c = []
        y_c = []
        for jj in index_centre:
            i_c = np.where(ll.index == jj)
            i_centre.append(i_c)
            x_c.append(ll[i_c].x / umeter)
            y_c.append(ll[i_c].y / umeter)
            for kk in arange(n_cells):
                x_kk = ll[kk].x / umeter
                y_kk = ll[kk].y / umeter
                if (x_kk - x_c[-1]) ** 2 + (y_kk - y_c[-1]) ** 2 <= (r_centre / umeter) ** 2:
                    ll[kk].index_up = True
                    indices_up.append(kk)
    else:
        print('ERROR - %s configuration not yet implemented' % type_connectivity)
    indices_up_layers.append(indices_up)

####################################
####################################
if __name__ == '__main__':
    # fig1 = figure(figsize=[15, 10], dpi=600.0)
    # ax = axes(projection='3d')
    # for layer in layers_exc:
    #     ax.scatter(layer.x / umeter, layer.y / umeter, layer.z / umeter, c='r', marker='.')
    # for layer in layers_inh:
    #     ax.scatter(layer.x / umeter, layer.y / umeter, layer.z / umeter, c='b', marker='.')
    #
    # ax.set_xlabel('X Label')
    # ax.set_ylabel('Y Label')
    # ax.set_zlabel('Z Label')
    # axis('auto')
    # tight_layout()
    # show()
    n_up = len(indices_up_layers[0])
    print(n_up)
    # print(indices_up_layers[0])
    # print(indices_up_layers[1])
    # for cc in arange(n_cells):
    #     print(l0[cc].index_up)
    #     if cc in indices_up:
    #         print('check True')

    fig2 = figure(2, figsize=[5*cm_fig_conv, 5*cm_fig_conv], dpi=600.0)
    plot(l0_exc.x / umeter, l0_exc.y / umeter, 'r.', l0_inh.x / umeter, l0_inh.y / umeter, 'b.', markersize=0.6)
    for cc in arange(n_cells):
        if cc in indices_up_layers[0]:
            l0[cc].index_up = True
            plot(l0[cc].x / umeter, l0[cc].y / umeter, 'go', markersize=1)
            # print('value: ' + str(l0[kk].x / umeter) + ' ' + str(l0[kk].y / umeter))
            # count += 1
    # if type_connectivity == 'central' or type_connectivity == 'multiple_centres':
    #     for vv in i_centre:
    #         plot(l0[vv].x / umeter, l0[vv].y / umeter, 'k*', markersize=10)

    xlabel('x (µm)')
    ylabel('y (µm)', rotation='vertical')
    axis('equal')
    gca().spines['top'].set_visible(False)
    gca().spines['right'].set_visible(False)
    tight_layout()
    show()

    # fig1.savefig('3D layout', dpi=600, transparent=True)
    fig2.savefig('3D-connectivity_random_newnewnewnew', dpi=600, transparent=True)
    print("Lif_SpatialConfiguration - END")
else:
    print("Finished Lif_SpatialConfiguration ")
####################################
####################################
