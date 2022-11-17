from Lif_SpatialConfiguration import *
from ExternalFunctions import connecting_int, createWeightMatrix, connecting_3D, connecting_ext
import random
import itertools

seed(val)


def conneting_population(pop: NeuronGroup, exc: NeuronGroup, inh: NeuronGroup, p_max: float):
    """
      Connecting the single layer. Eliminates excessive connections
      :param pop: the whole population
      :param exc: excitatory neurons
      :param inh: inhibitory neurons
      :param p_max: maximus probability of connection between 2 neurons
      :return: s_exc_all, s_inh_exc
      """

    s_exc_all = connecting_int(exc, pop, 'I_synEX+= alpha_ex * w', sigma, V_ax_P, pMAX, 1000 * mV)  # modificata
    # s_exc_all = connecting_int(exc, pop, 'I_synEX+= alpha_ex * w', sigma, V_ax_P, pMAX, 1000 * mV, val)
    s_inh_exc = connecting_int(inh, exc, 'I_synIN-= alpha_in * w', sigma, V_ax_I, pMAX, 1000 * mV)  # modificata
    # s_inh_exc = connecting_int(inh, exc, 'I_synIN-= alpha_in * w', sigma, V_ax_I, pMAX, 1000 * mV, val)
    for ss in [s_exc_all, s_inh_exc]:
        ll = len(ss.N_outgoing_pre)
        diff = 0
        count = 0
        k = 0

        for idx in range(ll):
            count = count + ss.N_outgoing_pre[idx]
            if ss.N_outgoing_pre[idx] > max_conn_int[idx]:
                diff = ss.N_outgoing_pre[idx] - max_conn_int[idx]  # connessioni da eliminare
                k = count - diff
                r = []
                r = random.sample(range(ss.N_outgoing_pre[idx]),
                                  diff)  # prendo tanti indici quanto le connessioni da mettere a peso zero
                for index1 in range(len(r)):
                    ss.w[r[index1] + k - max_conn_int[idx]] = 0
    return s_exc_all, s_inh_exc


def conneting_layers(pop_source: NeuronGroup, exc_source: NeuronGroup, inh_source: NeuronGroup,
                     pop_target: NeuronGroup, exc_target: NeuronGroup, inh_target: NeuronGroup,
                     p_max: float):
    """
      Connecting the single layer. Eliminates excessive connections
      :param pop_source: source whole population
      :param exc_source: source excitatory neurons
      :param inh_source: source inhibitory neurons
      :param pop_target: target whole population
      :param exc_target: target excitatory neurons
      :param inh_target: target inhibitory neurons
      :param p_max: maximus probability of connection between 2 neurons
      :return: s_exc_all, s_inh_exc
      """
    s_exc_all = connecting_3D(exc_source, pop_target, 'I_synEX+=alpha_ex * w', sigma, sigma_3D, V_ax_P, p_max,
                              1000 * mV)
    s_inh_exc = connecting_3D(inh_source, exc_target, 'I_synIN-=alpha_in * w', sigma, sigma_3D, V_ax_I, p_max,
                              1000 * mV)
    # s_exc_all = connecting_ext(exc_source, pop_target, 'I_synEX+=alpha_ex * w', sigma, sigma_3D, V_ax_P, p_max, 1000 * mV, val)
    # s_inh_exc = connecting_ext(inh_source, exc_target, 'I_synIN-=alpha_in * w', sigma, sigma_3D, V_ax_I, p_max, 1000 * mV, val)
    # CAMBIATA
    # for ss in [s_exc_all, s_inh_exc]:
    #     ll = len(ss.N_outgoing_pre)
    #     diff = 0
    #     count = 0
    #     k = 0
    #     for idx in range(ll):
    #         count = count + ss.N_outgoing_pre[idx]
    #         if ss.N_outgoing_pre[idx] > max_conn_ext[idx]:
    #             diff = ss.N_outgoing_pre[idx] - max_conn_ext[idx]  # connessioni da eliminare
    #             k = count - diff
    #             r = []
    #             r = random.sample(range(ss.N_outgoing_pre[idx]), diff)  # prendo tanti indici quanto le connessioni da mettere a peso zero
    #             for index1 in range(len(r)):
    #                 ss.w[r[index1] + k - max_conn_ext[idx]] = 0
    ll = len(s_exc_all.N_outgoing_pre)
    diff = 0
    count = 0
    k = 0
    for idx in range(ll):
        count = count + s_exc_all.N_outgoing_pre[idx]
        if s_exc_all.N_outgoing_pre[idx] > max_conn_ext_exc[idx]:
            diff = s_exc_all.N_outgoing_pre[idx] - max_conn_ext_exc[idx]  # connessioni da eliminare
            k = count - diff
            r = []
            r = random.sample(range(s_exc_all.N_outgoing_pre[idx]),
                              diff)  # prendo tanti indici quanto le connessioni da mettere a peso zero
            for index1 in range(len(r)):
                s_exc_all.w[r[index1] + k - max_conn_ext_exc[idx]] = 0

    ll = len(s_inh_exc.N_outgoing_pre)
    diff = 0
    count = 0
    k = 0
    for idx in range(ll):
        count = count + s_inh_exc.N_outgoing_pre[idx]
        if s_inh_exc.N_outgoing_pre[idx] > max_conn_ext_inh[idx]:
            diff = s_inh_exc.N_outgoing_pre[idx] - max_conn_ext_inh[idx]  # connessioni da eliminare
            k = count - diff
            r = []
            r = random.sample(range(s_inh_exc.N_outgoing_pre[idx]),
                              diff)  # prendo tanti indici quanto le connessioni da mettere a peso zero
            for index1 in range(len(r)):
                s_inh_exc.w[r[index1] + k - max_conn_ext_inh[idx]] = 0
    return s_exc_all, s_inh_exc


def connecting_modules(pop_source: NeuronGroup, exc_source: NeuronGroup, inh_source: NeuronGroup,
                            pop_target: NeuronGroup, exc_target: NeuronGroup, inh_target: NeuronGroup,
                            p_max: float):
    """
      Connecting the single layer. Eliminates excessive connections
      :param pop_source: source whole population
      :param exc_source: source excitatory neurons
      :param inh_source: source inhibitory neurons
      :param pop_target: target whole population
      :param exc_target: target excitatory neurons
      :param inh_target: target inhibitory neurons
      :param p_max: maximus probability of connection between 2 neurons
      :return: s_exc_all, s_inh_exc
      """
    s_exc_all = connecting_ext(exc_source, pop_target, 'I_synEX+=alpha_ex * w', sigma_2D_x, sigma_2D_y, V_ax_P, p_max,
                              1000 * mV)
    s_inh_exc = connecting_ext(inh_source, exc_target, 'I_synIN-=alpha_in * w', sigma_2D_x, sigma_2D_y, V_ax_I, p_max,
                              1000 * mV)

    ll = len(s_exc_all.N_outgoing_pre)
    diff = 0
    count = 0
    k = 0
    for idx in range(ll):
        count = count + s_exc_all.N_outgoing_pre[idx]
        if s_exc_all.N_outgoing_pre[idx] > max_conn_ext_exc[idx]:
            diff = s_exc_all.N_outgoing_pre[idx] - max_conn_ext_exc[idx]  # connessioni da eliminare
            k = count - diff
            r = []
            r = random.sample(range(s_exc_all.N_outgoing_pre[idx]),
                              diff)  # prendo tanti indici quanto le connessioni da mettere a peso zero
            for index1 in range(len(r)):
                s_exc_all.w[r[index1] + k - max_conn_ext_exc[idx]] = 0

    ll = len(s_inh_exc.N_outgoing_pre)
    diff = 0
    count = 0
    k = 0
    for idx in range(ll):
        count = count + s_inh_exc.N_outgoing_pre[idx]
        if s_inh_exc.N_outgoing_pre[idx] > max_conn_ext_inh[idx]:
            diff = s_inh_exc.N_outgoing_pre[idx] - max_conn_ext_inh[idx]  # connessioni da eliminare
            k = count - diff
            r = []
            r = random.sample(range(s_inh_exc.N_outgoing_pre[idx]),
                              diff)  # prendo tanti indici quanto le connessioni da mettere a peso zero
            for index1 in range(len(r)):
                s_inh_exc.w[r[index1] + k - max_conn_ext_inh[idx]] = 0
    return s_exc_all, s_inh_exc


###########################################
#       CONNECTIVITY Implementation
###########################################
num_conn_l0_sx, num_conn_l0_dx = np.zeros([n_cells, num_layers]), np.zeros([n_cells, num_layers])
num_conn_l1_sx, num_conn_l1_dx = np.zeros([n_cells, num_layers]), np.zeros([n_cells, num_layers])
num_conn_l2_sx, num_conn_l2_dx = np.zeros([n_cells, num_layers]), np.zeros([n_cells, num_layers])
num_conn_l3_sx, num_conn_l3_dx = np.zeros([n_cells, num_layers]), np.zeros([n_cells, num_layers])
num_conn_l4_sx, num_conn_l4_dx = np.zeros([n_cells, num_layers]), np.zeros([n_cells, num_layers])
num_conn_l5_sx, num_conn_l5_dx = np.zeros([n_cells, num_layers]), np.zeros([n_cells, num_layers])

# s_within = []
# for (ll, lexc, linh) in zip(layers, layers_exc, layers_inh):
#     s_ea, s_ie = conneting_population(ll, lexc, linh, pMAX)
#     s_within.append(s_ea, s_ie)

# Layer 0
# SX
s_exc_all_l0_sx, s_inh_exc_l0_sx = conneting_population(l0_sx, l0_sx_exc, l0_sx_inh, pMAX)
WeightMatrix_s_exc_all_l0_sx, num_conn_e = createWeightMatrix(l0_sx_exc, l0_sx, s_exc_all_l0_sx)
WeightMatrix_s_inh_exc_l0_sx, num_conn_i = createWeightMatrix(l0_sx_inh, l0_sx_exc, s_inh_exc_l0_sx)
num_conn_l0_sx[:, 0] = num_conn_e + num_conn_i
# DX
s_exc_all_l0_dx, s_inh_exc_l0_dx = conneting_population(l0_dx, l0_dx_exc, l0_dx_inh, pMAX)
WeightMatrix_s_exc_all_l0_dx, num_conn_e = createWeightMatrix(l0_dx_exc, l0_dx, s_exc_all_l0_dx)
WeightMatrix_s_inh_exc_l0_dx, num_conn_i = createWeightMatrix(l0_dx_inh, l0_dx_exc, s_inh_exc_l0_dx)
num_conn_l0_dx[:, 0] = num_conn_e + num_conn_i

# Connecting modules (through l0: sx->dx)
s_exc_all_l0_sxdx, s_inh_exc_l0_sxdx = conneting_layers(l0_sx, l0_sx_exc, l0_sx_inh,
                                                        l0_dx, l0_dx_exc, l0_dx_inh, pMAXext)
WeightMatrix_s_exc_all_l0_sxdx, num_conn_e = createWeightMatrix(l0_sx_exc, l0_dx, s_exc_all_l0_sxdx)
WeightMatrix_s_inh_exc_l0_sxdx, num_conn_i = createWeightMatrix(l0_sx_inh, l0_dx_exc, s_inh_exc_l0_sxdx)
# num_conn_l0_dx[:, 0] = num_conn_e + num_conn_i

# Connecting modules (through l0: dx->sx)
s_exc_all_l0_dxsx, s_inh_exc_l0_dxsx = conneting_layers(l0_dx, l0_dx_exc, l0_dx_inh,
                                                        l0_sx, l0_sx_exc, l0_sx_inh, pMAXext)
WeightMatrix_s_exc_all_l0_dxsx, num_conn_e = createWeightMatrix(l0_dx_exc, l0_sx, s_exc_all_l0_dxsx)
WeightMatrix_s_inh_exc_l0_dxsx, num_conn_i = createWeightMatrix(l0_dx_inh, l0_sx_exc, s_inh_exc_l0_dxsx)
# num_conn_l0_dx[:, 0] = num_conn_e + num_conn_i


# Layer 1
if num_layers > 1:
    # SX
    s_exc_all_l1_sx, s_inh_exc_l1_sx = conneting_population(l1_sx, l1_sx_exc, l1_sx_inh, pMAX)
    WeightMatrix_s_exc_all_l1_sx, num_conn_e = createWeightMatrix(l1_sx_exc, l1_sx, s_exc_all_l1_sx)
    WeightMatrix_s_inh_exc_l1_sx, num_conn_i = createWeightMatrix(l1_sx_inh, l1_sx_exc, s_inh_exc_l1_sx)
    num_conn_l1_sx[:, 1] = num_conn_e + num_conn_i
    if upward_conn:
        # Layer 0 to Layer 1
        s_exc_all_l0l1_sx, s_inh_exc_l0l1_sx = conneting_layers(l0_sx, l0_sx_exc, l0_sx_inh,
                                                                l1_sx, l1_sx_exc, l1_sx_inh, pMAXext)
        WeightMatrix_s_exc_all_l0l1_sx, num_conn_e = createWeightMatrix(l0_sx_exc, l1_sx, s_exc_all_l0l1_sx)
        WeightMatrix_s_inh_exc_l0l1_sx, num_conn_i = createWeightMatrix(l0_sx_inh, l1_sx_exc, s_inh_exc_l0l1_sx)
        num_conn_l0_sx[:, 1] = num_conn_e + num_conn_i
    if downward_conn:
        # Layer 1 to Layer 0
        s_exc_all_l1l0_sx, s_inh_exc_l1l0_sx = conneting_layers(l1_sx, l1_sx_exc, l1_sx_inh,
                                                                l0_sx, l0_sx_exc, l0_sx_inh, pMAXext)
        WeightMatrix_s_exc_all_l1l0_sx, num_conn_e = createWeightMatrix(l1_sx_exc, l0_sx, s_exc_all_l1l0_sx)
        WeightMatrix_s_inh_exc_l1l0_sx, num_conn_i = createWeightMatrix(l1_sx_inh, l0_sx_exc, s_inh_exc_l1l0_sx)
        num_conn_l1_sx[:, 0] = num_conn_e + num_conn_i
    # DX
    s_exc_all_l1_dx, s_inh_exc_l1_dx = conneting_population(l1_dx, l1_dx_exc, l1_dx_inh, pMAX)
    WeightMatrix_s_exc_all_l1_dx, num_conn_e = createWeightMatrix(l1_dx_exc, l1_dx, s_exc_all_l1_dx)
    WeightMatrix_s_inh_exc_l1_dx, num_conn_i = createWeightMatrix(l1_dx_inh, l1_dx_exc, s_inh_exc_l1_dx)
    num_conn_l1_dx[:, 1] = num_conn_e + num_conn_i
    if upward_conn:
        # Layer 0 to Layer 1
        s_exc_all_l0l1_dx, s_inh_exc_l0l1_dx = conneting_layers(l0_dx, l0_dx_exc, l0_dx_inh,
                                                                l1_dx, l1_dx_exc, l1_dx_inh, pMAXext)
        WeightMatrix_s_exc_all_l0l1_dx, num_conn_e = createWeightMatrix(l0_dx_exc, l1_dx, s_exc_all_l0l1_dx)
        WeightMatrix_s_inh_exc_l0l1_dx, num_conn_i = createWeightMatrix(l0_dx_inh, l1_dx_exc, s_inh_exc_l0l1_dx)
        num_conn_l0_dx[:, 1] = num_conn_e + num_conn_i
    if downward_conn:
        # Layer 1 to Layer 0
        s_exc_all_l1l0_dx, s_inh_exc_l1l0_dx = conneting_layers(l1_dx, l1_dx_exc, l1_dx_inh,
                                                                l0_dx, l0_dx_exc, l0_dx_inh, pMAXext)
        WeightMatrix_s_exc_all_l1l0_dx, num_conn_e = createWeightMatrix(l1_dx_exc, l0_dx, s_exc_all_l1l0_dx)
        WeightMatrix_s_inh_exc_l1l0_dx, num_conn_i = createWeightMatrix(l1_dx_inh, l0_dx_exc, s_inh_exc_l1l0_dx)
        num_conn_l1_dx[:, 0] = num_conn_e + num_conn_i

# Layer 2
if num_layers > 2:
    # SX
    s_exc_all_l2_sx, s_inh_exc_l2_sx = conneting_population(l2_sx, l2_sx_exc, l2_sx_inh, pMAX)
    WeightMatrix_s_exc_all_l2_sx, num_conn_e = createWeightMatrix(l2_sx_exc, l2_sx, s_exc_all_l2_sx)
    WeightMatrix_s_inh_exc_l2_sx, num_conn_i = createWeightMatrix(l2_sx_inh, l2_sx_exc, s_inh_exc_l2_sx)
    num_conn_l2_sx[:, 2] = num_conn_e + num_conn_i
    if upward_conn:
        # Layer 0 to Layer 2
        s_exc_all_l0l2_sx, s_inh_exc_l0l2_sx = conneting_layers(l0_sx, l0_sx_exc, l0_sx_inh,
                                                          l2_sx, l2_sx_exc, l2_sx_inh, pMAXext)
        WeightMatrix_s_exc_all_l0l2_sx, num_conn_e = createWeightMatrix(l0_sx_exc, l2_sx, s_exc_all_l0l2_sx)
        WeightMatrix_s_inh_exc_l0l2_sx, num_conn_i = createWeightMatrix(l0_sx_inh, l2_sx_exc, s_inh_exc_l0l2_sx)
        num_conn_l0_sx[:, 2] = num_conn_e + num_conn_i

        # Layer 1 to Layer 2
        s_exc_all_l1l2_sx, s_inh_exc_l1l2_sx = conneting_layers(l1_sx, l1_sx_exc, l1_sx_inh,
                                                          l2_sx, l2_sx_exc, l2_sx_inh, pMAXext)
        WeightMatrix_s_exc_all_l1l2_sx, num_conn_e = createWeightMatrix(l1_sx_exc, l2_sx, s_exc_all_l1l2_sx)
        WeightMatrix_s_inh_exc_l1l2_sx, num_conn_i = createWeightMatrix(l1_sx_inh, l2_sx_exc, s_inh_exc_l1l2_sx)
        num_conn_l1_sx[:, 2] = num_conn_e + num_conn_i
    if downward_conn:
        # Layer 2 to Layer 0
        s_exc_all_l2l0_sx, s_inh_exc_l2l0_sx = conneting_layers(l2_sx, l2_sx_exc, l2_sx_inh,
                                                          l0_sx, l0_sx_exc, l0_sx_inh, pMAXext)
        WeightMatrix_s_exc_all_l2l0_sx, num_conn_e = createWeightMatrix(l2_sx_exc, l0_sx, s_exc_all_l2l0_sx)
        WeightMatrix_s_inh_exc_l2l0_sx, num_conn_i = createWeightMatrix(l2_sx_inh, l0_sx_exc, s_inh_exc_l2l0_sx)
        num_conn_l2_sx[:, 0] = num_conn_e + num_conn_i
        # Layer 2 to Layer 1
        s_exc_all_l2l1_sx, s_inh_exc_l2l1_sx = conneting_layers(l2_sx, l2_sx_exc, l2_sx_inh,
                                                          l1_sx, l1_sx_exc, l1_sx_inh, pMAXext)
        WeightMatrix_s_exc_all_l2l1_sx, num_conn_e = createWeightMatrix(l2_sx_exc, l1_sx, s_exc_all_l2l1_sx)
        WeightMatrix_s_inh_exc_l2l1_sx, num_conn_i = createWeightMatrix(l2_sx_inh, l1_sx_exc, s_inh_exc_l2l1_sx)
        num_conn_l2_sx[:, 1] = num_conn_e + num_conn_i
    # DX
    s_exc_all_l2_dx, s_inh_exc_l2_dx = conneting_population(l2_dx, l2_dx_exc, l2_dx_inh, pMAX)
    WeightMatrix_s_exc_all_l2_dx, num_conn_e = createWeightMatrix(l2_dx_exc, l2_dx, s_exc_all_l2_dx)
    WeightMatrix_s_inh_exc_l2_dx, num_conn_i = createWeightMatrix(l2_dx_inh, l2_dx_exc, s_inh_exc_l2_dx)
    num_conn_l2_dx[:, 2] = num_conn_e + num_conn_i
    if upward_conn:
        # Layer 0 to Layer 2
        s_exc_all_l0l2_dx, s_inh_exc_l0l2_dx = conneting_layers(l0_dx, l0_dx_exc, l0_dx_inh,
                                                                l2_dx, l2_dx_exc, l2_dx_inh, pMAXext)
        WeightMatrix_s_exc_all_l0l2_dx, num_conn_e = createWeightMatrix(l0_dx_exc, l2_dx, s_exc_all_l0l2_dx)
        WeightMatrix_s_inh_exc_l0l2_dx, num_conn_i = createWeightMatrix(l0_dx_inh, l2_dx_exc, s_inh_exc_l0l2_dx)
        num_conn_l0_dx[:, 2] = num_conn_e + num_conn_i

        # Layer 1 to Layer 2
        s_exc_all_l1l2_dx, s_inh_exc_l1l2_dx = conneting_layers(l1_dx, l1_dx_exc, l1_dx_inh,
                                                                l2_dx, l2_dx_exc, l2_dx_inh, pMAXext)
        WeightMatrix_s_exc_all_l1l2_dx, num_conn_e = createWeightMatrix(l1_dx_exc, l2_dx, s_exc_all_l1l2_dx)
        WeightMatrix_s_inh_exc_l1l2_dx, num_conn_i = createWeightMatrix(l1_dx_inh, l2_dx_exc, s_inh_exc_l1l2_dx)
        num_conn_l1_dx[:, 2] = num_conn_e + num_conn_i
    if downward_conn:
        # Layer 2 to Layer 0
        s_exc_all_l2l0_dx, s_inh_exc_l2l0_dx = conneting_layers(l2_dx, l2_dx_exc, l2_dx_inh,
                                                                l0_dx, l0_dx_exc, l0_dx_inh, pMAXext)
        WeightMatrix_s_exc_all_l2l0_dx, num_conn_e = createWeightMatrix(l2_dx_exc, l0_dx, s_exc_all_l2l0_dx)
        WeightMatrix_s_inh_exc_l2l0_dx, num_conn_i = createWeightMatrix(l2_dx_inh, l0_dx_exc, s_inh_exc_l2l0_dx)
        num_conn_l2_dx[:, 0] = num_conn_e + num_conn_i
        # Layer 2 to Layer 1
        s_exc_all_l2l1_dx, s_inh_exc_l2l1_dx = conneting_layers(l2_dx, l2_dx_exc, l2_dx_inh,
                                                                l1_dx, l1_dx_exc, l1_dx_inh, pMAXext)
        WeightMatrix_s_exc_all_l2l1_dx, num_conn_e = createWeightMatrix(l2_dx_exc, l1_dx, s_exc_all_l2l1_dx)
        WeightMatrix_s_inh_exc_l2l1_dx, num_conn_i = createWeightMatrix(l2_dx_inh, l1_dx_exc, s_inh_exc_l2l1_dx)
        num_conn_l2_dx[:, 1] = num_conn_e + num_conn_i

# Layer 3
if num_layers > 3:
    # SX
    s_exc_all_l3_sx, s_inh_exc_l3_sx = conneting_population(l3_sx, l3_sx_exc, l3_sx_inh, pMAX)
    WeightMatrix_s_exc_all_l3_sx, num_conn_e = createWeightMatrix(l3_sx_exc, l3_sx, s_exc_all_l3_sx)
    WeightMatrix_s_inh_exc_l3_sx, num_conn_i = createWeightMatrix(l3_sx_inh, l3_sx_exc, s_inh_exc_l3_sx)
    num_conn_l3_sx[:, 3] = num_conn_e + num_conn_i
    if upward_conn:
        # Layer 0 to Layer 3
        s_exc_all_l0l3_sx, s_inh_exc_l0l3_sx = conneting_layers(l0_sx, l0_sx_exc, l0_sx_inh,
                                                          l3_sx, l3_sx_exc, l3_sx_inh, pMAXext)
        WeightMatrix_s_exc_all_l0l3_sx, num_conn_e = createWeightMatrix(l0_sx_exc, l3_sx, s_exc_all_l0l3_sx)
        WeightMatrix_s_inh_exc_l0l3_sx, num_conn_i = createWeightMatrix(l0_sx_inh, l3_sx_exc, s_inh_exc_l0l3_sx)
        num_conn_l0_sx[:, 3] = num_conn_e + num_conn_i
        # Layer 1 to Layer 3
        s_exc_all_l1l3_sx, s_inh_exc_l1l3_sx = conneting_layers(l1_sx, l1_sx_exc, l1_sx_inh,
                                                          l3_sx, l3_sx_exc, l3_sx_inh, pMAXext)
        WeightMatrix_s_exc_all_l1l3_sx, num_conn_e = createWeightMatrix(l1_sx_exc, l3_sx, s_exc_all_l1l3_sx)
        WeightMatrix_s_inh_exc_l1l3_sx, num_conn_i = createWeightMatrix(l1_sx_inh, l3_sx_exc, s_inh_exc_l1l3_sx)
        num_conn_l1_sx[:, 3] = num_conn_e + num_conn_i
        # Layer 2 to Layer 3
        s_exc_all_l2l3_sx, s_inh_exc_l2l3_sx = conneting_layers(l2_sx, l2_sx_exc, l2_sx_inh,
                                                          l3_sx, l3_sx_exc, l3_sx_inh, pMAXext)
        WeightMatrix_s_exc_all_l2l3_sx, num_conn_e = createWeightMatrix(l2_sx_exc, l3_sx, s_exc_all_l2l3_sx)
        WeightMatrix_s_inh_exc_l2l3_sx, num_conn_i = createWeightMatrix(l2_sx_inh, l3_sx_exc, s_inh_exc_l2l3_sx)
        num_conn_l2_sx[:, 3] = num_conn_e + num_conn_i
    if downward_conn:
        # Layer 3 to Layer 0
        s_exc_all_l3l0_sx, s_inh_exc_l3l0_sx = conneting_layers(l3_sx, l3_sx_exc, l3_sx_inh,
                                                          l0_sx, l0_sx_exc, l0_sx_inh, pMAXext)
        WeightMatrix_s_exc_all_l3l0_sx, num_conn_e = createWeightMatrix(l3_sx_exc, l0_sx, s_exc_all_l3l0_sx)
        WeightMatrix_s_inh_exc_l3l0_sx, num_conn_i = createWeightMatrix(l3_sx_inh, l0_sx_exc, s_inh_exc_l3l0_sx)
        num_conn_l3_sx[:, 0] = num_conn_e + num_conn_i
        # Layer 3 to Layer 1
        s_exc_all_l3l1_sx, s_inh_exc_l3l1_sx = conneting_layers(l3_sx, l3_sx_exc, l3_sx_inh,
                                                          l1_sx, l1_sx_exc, l1_sx_inh, pMAXext)
        WeightMatrix_s_exc_all_l3l1_sx, num_conn_e = createWeightMatrix(l3_sx_exc, l1_sx, s_exc_all_l3l1_sx)
        WeightMatrix_s_inh_exc_l3l1_sx, num_conn_i = createWeightMatrix(l3_sx_inh, l1_sx_exc, s_inh_exc_l3l1_sx)
        num_conn_l3_sx[:, 1] = num_conn_e + num_conn_i
        # Layer 3 to Layer 2
        s_exc_all_l3l2_sx, s_inh_exc_l3l2_sx = conneting_layers(l3_sx, l3_sx_exc, l3_sx_inh,
                                                          l2_sx, l2_sx_exc, l2_sx_inh, pMAXext)
        WeightMatrix_s_exc_all_l3l2_sx, num_conn_e = createWeightMatrix(l3_sx_exc, l2_sx, s_exc_all_l3l2_sx)
        WeightMatrix_s_inh_exc_l3l2_sx, num_conn_i = createWeightMatrix(l3_sx_inh, l2_sx_exc, s_inh_exc_l3l2_sx)
        num_conn_l3_sx[:, 2] = num_conn_e + num_conn_i
    # DX
    s_exc_all_l3_dx, s_inh_exc_l3_dx = conneting_population(l3_dx, l3_dx_exc, l3_dx_inh, pMAX)
    WeightMatrix_s_exc_all_l3_dx, num_conn_e = createWeightMatrix(l3_dx_exc, l3_dx, s_exc_all_l3_dx)
    WeightMatrix_s_inh_exc_l3_dx, num_conn_i = createWeightMatrix(l3_dx_inh, l3_dx_exc, s_inh_exc_l3_dx)
    num_conn_l3_dx[:, 3] = num_conn_e + num_conn_i
    if upward_conn:
        # Layer 0 to Layer 3
        s_exc_all_l0l3_dx, s_inh_exc_l0l3_dx = conneting_layers(l0_dx, l0_dx_exc, l0_dx_inh,
                                                                l3_dx, l3_dx_exc, l3_dx_inh, pMAXext)
        WeightMatrix_s_exc_all_l0l3_dx, num_conn_e = createWeightMatrix(l0_dx_exc, l3_dx, s_exc_all_l0l3_dx)
        WeightMatrix_s_inh_exc_l0l3_dx, num_conn_i = createWeightMatrix(l0_dx_inh, l3_dx_exc, s_inh_exc_l0l3_dx)
        num_conn_l0_dx[:, 3] = num_conn_e + num_conn_i
        # Layer 1 to Layer 3
        s_exc_all_l1l3_dx, s_inh_exc_l1l3_dx = conneting_layers(l1_dx, l1_dx_exc, l1_dx_inh,
                                                                l3_dx, l3_dx_exc, l3_dx_inh, pMAXext)
        WeightMatrix_s_exc_all_l1l3_dx, num_conn_e = createWeightMatrix(l1_dx_exc, l3_dx, s_exc_all_l1l3_dx)
        WeightMatrix_s_inh_exc_l1l3_dx, num_conn_i = createWeightMatrix(l1_dx_inh, l3_dx_exc, s_inh_exc_l1l3_dx)
        num_conn_l1_dx[:, 3] = num_conn_e + num_conn_i
        # Layer 2 to Layer 3
        s_exc_all_l2l3_dx, s_inh_exc_l2l3_dx = conneting_layers(l2_dx, l2_dx_exc, l2_dx_inh,
                                                                l3_dx, l3_dx_exc, l3_dx_inh, pMAXext)
        WeightMatrix_s_exc_all_l2l3_dx, num_conn_e = createWeightMatrix(l2_dx_exc, l3_dx, s_exc_all_l2l3_dx)
        WeightMatrix_s_inh_exc_l2l3_dx, num_conn_i = createWeightMatrix(l2_dx_inh, l3_dx_exc, s_inh_exc_l2l3_dx)
        num_conn_l2_dx[:, 3] = num_conn_e + num_conn_i
    if downward_conn:
        # Layer 3 to Layer 0
        s_exc_all_l3l0_dx, s_inh_exc_l3l0_dx = conneting_layers(l3_dx, l3_dx_exc, l3_dx_inh,
                                                                l0_dx, l0_dx_exc, l0_dx_inh, pMAXext)
        WeightMatrix_s_exc_all_l3l0_dx, num_conn_e = createWeightMatrix(l3_dx_exc, l0_dx, s_exc_all_l3l0_dx)
        WeightMatrix_s_inh_exc_l3l0_dx, num_conn_i = createWeightMatrix(l3_dx_inh, l0_dx_exc, s_inh_exc_l3l0_dx)
        num_conn_l3_dx[:, 0] = num_conn_e + num_conn_i
        # Layer 3 to Layer 1
        s_exc_all_l3l1_dx, s_inh_exc_l3l1_dx = conneting_layers(l3_dx, l3_dx_exc, l3_dx_inh,
                                                                l1_dx, l1_dx_exc, l1_dx_inh, pMAXext)
        WeightMatrix_s_exc_all_l3l1_dx, num_conn_e = createWeightMatrix(l3_dx_exc, l1_dx, s_exc_all_l3l1_dx)
        WeightMatrix_s_inh_exc_l3l1_dx, num_conn_i = createWeightMatrix(l3_dx_inh, l1_dx_exc, s_inh_exc_l3l1_dx)
        num_conn_l3_dx[:, 1] = num_conn_e + num_conn_i
        # Layer 3 to Layer 2
        s_exc_all_l3l2_dx, s_inh_exc_l3l2_dx = conneting_layers(l3_dx, l3_dx_exc, l3_dx_inh,
                                                                l2_dx, l2_dx_exc, l2_dx_inh, pMAXext)
        WeightMatrix_s_exc_all_l3l2_dx, num_conn_e = createWeightMatrix(l3_dx_exc, l2_dx, s_exc_all_l3l2_dx)
        WeightMatrix_s_inh_exc_l3l2_dx, num_conn_i = createWeightMatrix(l3_dx_inh, l2_dx_exc, s_inh_exc_l3l2_dx)
        num_conn_l3_dx[:, 2] = num_conn_e + num_conn_i

# Layer 4
if num_layers > 4:
    # SX
    s_exc_all_l4_sx, s_inh_exc_l4_sx = conneting_population(l4_sx, l4_sx_exc, l4_sx_inh, pMAX)
    WeightMatrix_s_exc_all_l4_sx, num_conn_e = createWeightMatrix(l4_sx_exc, l4_sx, s_exc_all_l4_sx)
    WeightMatrix_s_inh_exc_l4_sx, num_conn_i = createWeightMatrix(l4_sx_inh, l4_sx_exc, s_inh_exc_l4_sx)
    num_conn_l4_sx[:, 4] = num_conn_e + num_conn_i
    if upward_conn:
        # Layer 0 to Layer 4
        s_exc_all_l0l4_sx, s_inh_exc_l0l4_sx = conneting_layers(l0_sx, l0_sx_exc, l0_sx_inh,
                                                          l4_sx, l4_sx_exc, l4_sx_inh, pMAXext)
        WeightMatrix_s_exc_all_l0l4_sx, num_conn_e = createWeightMatrix(l0_sx_exc, l4_sx, s_exc_all_l0l4_sx)
        WeightMatrix_s_inh_exc_l0l4_sx, num_conn_i = createWeightMatrix(l0_sx_inh, l4_sx_exc, s_inh_exc_l0l4_sx)
        num_conn_l0_sx[:, 4] = num_conn_e + num_conn_i
        # Layer 1 to Layer 4
        s_exc_all_l1l4_sx, s_inh_exc_l1l4_sx = conneting_layers(l1_sx, l1_sx_exc, l1_sx_inh,
                                                          l4_sx, l4_sx_exc, l4_sx_inh, pMAXext)
        WeightMatrix_s_exc_all_l1l4_sx, num_conn_e = createWeightMatrix(l1_sx_exc, l4_sx, s_exc_all_l1l4_sx)
        WeightMatrix_s_inh_exc_l1l4_sx, num_conn_i = createWeightMatrix(l1_sx_inh, l4_sx_exc, s_inh_exc_l1l4_sx)
        num_conn_l1_sx[:, 4] = num_conn_e + num_conn_i
        # Layer 2 to Layer 4
        s_exc_all_l2l4_sx, s_inh_exc_l2l4_sx = conneting_layers(l2_sx, l2_sx_exc, l2_sx_inh,
                                                          l4_sx, l4_sx_exc, l4_sx_inh, pMAXext)
        WeightMatrix_s_exc_all_l2l4_sx, num_conn_e = createWeightMatrix(l2_sx_exc, l4_sx, s_exc_all_l2l4_sx)
        WeightMatrix_s_inh_exc_l2l4_sx, num_conn_i = createWeightMatrix(l2_sx_inh, l4_sx_exc, s_inh_exc_l2l4_sx)
        num_conn_l2_sx[:, 4] = num_conn_e + num_conn_i
        # Layer 3 to Layer 4
        s_exc_all_l3l4_sx, s_inh_exc_l3l4_sx = conneting_layers(l3_sx, l3_sx_exc, l3_sx_inh,
                                                          l4_sx, l4_sx_exc, l4_sx_inh, pMAXext)
        WeightMatrix_s_exc_all_l3l4_sx, num_conn_e = createWeightMatrix(l3_sx_exc, l4_sx, s_exc_all_l3l4_sx)
        WeightMatrix_s_inh_exc_l3l4_sx, num_conn_i = createWeightMatrix(l3_sx_inh, l4_sx_exc, s_inh_exc_l3l4_sx)
        num_conn_l3_sx[:, 4] = num_conn_e + num_conn_i
    if downward_conn:
        # Layer 4 to Layer 0
        s_exc_all_l4l0_sx, s_inh_exc_l4l0_sx = conneting_layers(l4_sx, l4_sx_exc, l4_sx_inh,
                                                          l0_sx, l0_sx_exc, l0_sx_inh, pMAXext)
        WeightMatrix_s_exc_all_l4l0_sx, num_conn_e = createWeightMatrix(l4_sx_exc, l0_sx, s_exc_all_l4l0_sx)
        WeightMatrix_s_inh_exc_l4l0_sx, num_conn_i = createWeightMatrix(l4_sx_inh, l0_sx_exc, s_inh_exc_l4l0_sx)
        num_conn_l4_sx[:, 0] = num_conn_e + num_conn_i
        # Layer 4 to Layer 1
        s_exc_all_l4l1_sx, s_inh_exc_l4l1_sx = conneting_layers(l4_sx, l4_sx_exc, l4_sx_inh,
                                                          l1_sx, l1_sx_exc, l1_sx_inh, pMAXext)
        WeightMatrix_s_exc_all_l4l1_sx, num_conn_e = createWeightMatrix(l4_sx_exc, l1_sx, s_exc_all_l4l1_sx)
        WeightMatrix_s_inh_exc_l4l1_sx, num_conn_i = createWeightMatrix(l4_sx_inh, l1_sx_exc, s_inh_exc_l4l1_sx)
        num_conn_l4_sx[:, 1] = num_conn_e + num_conn_i
        # Layer 4 to Layer 2
        s_exc_all_l4l2_sx, s_inh_exc_l4l2_sx = conneting_layers(l4_sx, l4_sx_exc, l4_sx_inh,
                                                          l2_sx, l2_sx_exc, l2_sx_inh, pMAXext)
        WeightMatrix_s_exc_all_l4l2_sx, num_conn_e = createWeightMatrix(l4_sx_exc, l2_sx, s_exc_all_l4l2_sx)
        WeightMatrix_s_inh_exc_l4l2_sx, num_conn_i = createWeightMatrix(l4_sx_inh, l2_sx_exc, s_inh_exc_l4l2_sx)
        num_conn_l4_sx[:, 2] = num_conn_e + num_conn_i
        # Layer 4 to Layer 3
        s_exc_all_l4l3_sx, s_inh_exc_l4l3_sx = conneting_layers(l4_sx, l4_sx_exc, l4_sx_inh,
                                                          l3_sx, l3_sx_exc, l3_sx_inh, pMAXext)
        WeightMatrix_s_exc_all_l4l3_sx, num_conn_e = createWeightMatrix(l4_sx_exc, l3_sx, s_exc_all_l4l3_sx)
        WeightMatrix_s_inh_exc_l4l3_sx, num_conn_i = createWeightMatrix(l4_sx_inh, l3_sx_exc, s_inh_exc_l4l3_sx)
        num_conn_l4_sx[:, 3] = num_conn_e + num_conn_i
    # DX
    s_exc_all_l4_dx, s_inh_exc_l4_dx = conneting_population(l4_dx, l4_dx_exc, l4_dx_inh, pMAX)
    WeightMatrix_s_exc_all_l4_dx, num_conn_e = createWeightMatrix(l4_dx_exc, l4_dx, s_exc_all_l4_dx)
    WeightMatrix_s_inh_exc_l4_dx, num_conn_i = createWeightMatrix(l4_dx_inh, l4_dx_exc, s_inh_exc_l4_dx)
    num_conn_l4_dx[:, 4] = num_conn_e + num_conn_i
    if upward_conn:
        # Layer 0 to Layer 4
        s_exc_all_l0l4_dx, s_inh_exc_l0l4_dx = conneting_layers(l0_dx, l0_dx_exc, l0_dx_inh,
                                                                l4_dx, l4_dx_exc, l4_dx_inh, pMAXext)
        WeightMatrix_s_exc_all_l0l4_dx, num_conn_e = createWeightMatrix(l0_dx_exc, l4_dx, s_exc_all_l0l4_dx)
        WeightMatrix_s_inh_exc_l0l4_dx, num_conn_i = createWeightMatrix(l0_dx_inh, l4_dx_exc, s_inh_exc_l0l4_dx)
        num_conn_l0_dx[:, 4] = num_conn_e + num_conn_i
        # Layer 1 to Layer 4
        s_exc_all_l1l4_dx, s_inh_exc_l1l4_dx = conneting_layers(l1_dx, l1_dx_exc, l1_dx_inh,
                                                                l4_dx, l4_dx_exc, l4_dx_inh, pMAXext)
        WeightMatrix_s_exc_all_l1l4_dx, num_conn_e = createWeightMatrix(l1_dx_exc, l4_dx, s_exc_all_l1l4_dx)
        WeightMatrix_s_inh_exc_l1l4_dx, num_conn_i = createWeightMatrix(l1_dx_inh, l4_dx_exc, s_inh_exc_l1l4_dx)
        num_conn_l1_dx[:, 4] = num_conn_e + num_conn_i
        # Layer 2 to Layer 4
        s_exc_all_l2l4_dx, s_inh_exc_l2l4_dx = conneting_layers(l2_dx, l2_dx_exc, l2_dx_inh,
                                                                l4_dx, l4_dx_exc, l4_dx_inh, pMAXext)
        WeightMatrix_s_exc_all_l2l4_dx, num_conn_e = createWeightMatrix(l2_dx_exc, l4_dx, s_exc_all_l2l4_dx)
        WeightMatrix_s_inh_exc_l2l4_dx, num_conn_i = createWeightMatrix(l2_dx_inh, l4_dx_exc, s_inh_exc_l2l4_dx)
        num_conn_l2_dx[:, 4] = num_conn_e + num_conn_i
        # Layer 3 to Layer 4
        s_exc_all_l3l4_dx, s_inh_exc_l3l4_dx = conneting_layers(l3_dx, l3_dx_exc, l3_dx_inh,
                                                                l4_dx, l4_dx_exc, l4_dx_inh, pMAXext)
        WeightMatrix_s_exc_all_l3l4_dx, num_conn_e = createWeightMatrix(l3_dx_exc, l4_dx, s_exc_all_l3l4_dx)
        WeightMatrix_s_inh_exc_l3l4_dx, num_conn_i = createWeightMatrix(l3_dx_inh, l4_dx_exc, s_inh_exc_l3l4_dx)
        num_conn_l3_dx[:, 4] = num_conn_e + num_conn_i
    if downward_conn:
        # Layer 4 to Layer 0
        s_exc_all_l4l0_dx, s_inh_exc_l4l0_dx = conneting_layers(l4_dx, l4_dx_exc, l4_dx_inh,
                                                                l0_dx, l0_dx_exc, l0_dx_inh, pMAXext)
        WeightMatrix_s_exc_all_l4l0_dx, num_conn_e = createWeightMatrix(l4_dx_exc, l0_dx, s_exc_all_l4l0_dx)
        WeightMatrix_s_inh_exc_l4l0_dx, num_conn_i = createWeightMatrix(l4_dx_inh, l0_dx_exc, s_inh_exc_l4l0_dx)
        num_conn_l4_dx[:, 0] = num_conn_e + num_conn_i
        # Layer 4 to Layer 1
        s_exc_all_l4l1_dx, s_inh_exc_l4l1_dx = conneting_layers(l4_dx, l4_dx_exc, l4_dx_inh,
                                                                l1_dx, l1_dx_exc, l1_dx_inh, pMAXext)
        WeightMatrix_s_exc_all_l4l1_dx, num_conn_e = createWeightMatrix(l4_dx_exc, l1_dx, s_exc_all_l4l1_dx)
        WeightMatrix_s_inh_exc_l4l1_dx, num_conn_i = createWeightMatrix(l4_dx_inh, l1_dx_exc, s_inh_exc_l4l1_dx)
        num_conn_l4_dx[:, 1] = num_conn_e + num_conn_i
        # Layer 4 to Layer 2
        s_exc_all_l4l2_dx, s_inh_exc_l4l2_dx = conneting_layers(l4_dx, l4_dx_exc, l4_dx_inh,
                                                                l2_dx, l2_dx_exc, l2_dx_inh, pMAXext)
        WeightMatrix_s_exc_all_l4l2_dx, num_conn_e = createWeightMatrix(l4_dx_exc, l2_dx, s_exc_all_l4l2_dx)
        WeightMatrix_s_inh_exc_l4l2_dx, num_conn_i = createWeightMatrix(l4_dx_inh, l2_dx_exc, s_inh_exc_l4l2_dx)
        num_conn_l4_dx[:, 2] = num_conn_e + num_conn_i
        # Layer 4 to Layer 3
        s_exc_all_l4l3_dx, s_inh_exc_l4l3_dx = conneting_layers(l4_dx, l4_dx_exc, l4_dx_inh,
                                                                l3_dx, l3_dx_exc, l3_dx_inh, pMAXext)
        WeightMatrix_s_exc_all_l4l3_dx, num_conn_e = createWeightMatrix(l4_dx_exc, l3_dx, s_exc_all_l4l3_dx)
        WeightMatrix_s_inh_exc_l4l3_dx, num_conn_i = createWeightMatrix(l4_dx_inh, l3_dx_exc, s_inh_exc_l4l3_dx)
        num_conn_l4_dx[:, 3] = num_conn_e + num_conn_i

# Layer 5
if num_layers > 5:
    # SX
    s_exc_all_l5_sx, s_inh_exc_l5_sx = conneting_population(l5_sx, l5_sx_exc, l5_sx_inh, pMAX)
    WeightMatrix_s_exc_all_l5_sx, num_conn_e = createWeightMatrix(l5_sx_exc, l5_sx, s_exc_all_l5_sx)
    WeightMatrix_s_inh_exc_l5_sx, num_conn_i = createWeightMatrix(l5_sx_inh, l5_sx_exc, s_inh_exc_l5_sx)
    num_conn_l5_sx[:, 5] = num_conn_e + num_conn_i
    if upward_conn:
        # Layer 0 to Layer 5
        s_exc_all_l0l5_sx, s_inh_exc_l0l5_sx = conneting_layers(l0_sx, l0_sx_exc, l0_sx_inh,
                                                          l5_sx, l5_sx_exc, l5_sx_inh, pMAXext)
        WeightMatrix_s_exc_all_l0l5_sx, num_conn_e = createWeightMatrix(l0_sx_exc, l5_sx, s_exc_all_l0l5_sx)
        WeightMatrix_s_inh_exc_l0l5_sx, num_conn_i = createWeightMatrix(l0_sx_inh, l5_sx_exc, s_inh_exc_l0l5_sx)
        num_conn_l0_sx[:, 5] = num_conn_e + num_conn_i
        # Layer 1 to Layer 5
        s_exc_all_l1l5_sx, s_inh_exc_l1l5_sx = conneting_layers(l1_sx, l1_sx_exc, l1_sx_inh,
                                                          l5_sx, l5_sx_exc, l5_sx_inh, pMAXext)
        WeightMatrix_s_exc_all_l1l5_sx, num_conn_e = createWeightMatrix(l1_sx_exc, l5_sx, s_exc_all_l1l5_sx)
        WeightMatrix_s_inh_exc_l1l5_sx, num_conn_i = createWeightMatrix(l1_sx_inh, l5_sx_exc, s_inh_exc_l1l5_sx)
        num_conn_l1_sx[:, 5] = num_conn_e + num_conn_i
        # Layer 2 to Layer 5
        s_exc_all_l2l5_sx, s_inh_exc_l2l5_sx = conneting_layers(l2_sx, l2_sx_exc, l2_sx_inh,
                                                          l5_sx, l5_sx_exc, l5_sx_inh, pMAXext)
        WeightMatrix_s_exc_all_l2l5_sx, num_conn_e = createWeightMatrix(l2_sx_exc, l5_sx, s_exc_all_l2l5_sx)
        WeightMatrix_s_inh_exc_l2l5_sx, num_conn_i = createWeightMatrix(l2_sx_inh, l5_sx_exc, s_inh_exc_l2l5_sx)
        num_conn_l2_sx[:, 5] = num_conn_e + num_conn_i
        # Layer 3 to Layer 5
        s_exc_all_l3l5_sx, s_inh_exc_l3l5_sx = conneting_layers(l3_sx, l3_sx_exc, l3_sx_inh,
                                                          l5_sx, l5_sx_exc, l5_sx_inh, pMAXext)
        WeightMatrix_s_exc_all_l3l5_sx, num_conn_e = createWeightMatrix(l3_sx_exc, l5_sx, s_exc_all_l3l5_sx)
        WeightMatrix_s_inh_exc_l3l5_sx, num_conn_i = createWeightMatrix(l3_sx_inh, l5_sx_exc, s_inh_exc_l3l5_sx)
        num_conn_l3_sx[:, 5] = num_conn_e + num_conn_i
        # Layer 4 to Layer 5
        s_exc_all_l4l5_sx, s_inh_exc_l4l5_sx = conneting_layers(l4, l4_sx_exc, l4_sx_inh,
                                                          l5_sx, l5_sx_exc, l5_sx_inh, pMAXext)
        WeightMatrix_s_exc_all_l4l5_sx, num_conn_e = createWeightMatrix(l4_sx_exc, l5_sx, s_exc_all_l4l5_sx)
        WeightMatrix_s_inh_exc_l4l5_sx, num_conn_i = createWeightMatrix(l4_sx_inh, l5_sx_exc, s_inh_exc_l4l5_sx)
        num_conn_l4_sx[:, 5] = num_conn_e + num_conn_i
    if downward_conn:
        # Layer 5 to Layer 0
        s_exc_all_l5l0_sx, s_inh_exc_l5l0_sx = conneting_layers(l5_sx, l5_sx_exc, l5_sx_inh,
                                                          l0_sx, l0_sx_exc, l0_sx_inh, pMAXext)
        WeightMatrix_s_exc_all_l5l0_sx, num_conn_e = createWeightMatrix(l5_sx_exc, l0_sx, s_exc_all_l5l0_sx)
        WeightMatrix_s_inh_exc_l5l0_sx, num_conn_i = createWeightMatrix(l5_sx_inh, l0_sx_exc, s_inh_exc_l5l0_sx)
        num_conn_l5_sx[:, 0] = num_conn_e + num_conn_i
        # Layer 5 to Layer 1
        s_exc_all_l5l1_sx, s_inh_exc_l5l1_sx = conneting_layers(l5_sx, l5_sx_exc, l5_sx_inh,
                                                          l1_sx, l1_sx_exc, l1_sx_inh, pMAXext)
        WeightMatrix_s_exc_all_l5l1_sx, num_conn_e = createWeightMatrix(l5_sx_exc, l1_sx, s_exc_all_l5l1_sx)
        WeightMatrix_s_inh_exc_l5l1_sx, num_conn_i = createWeightMatrix(l5_sx_inh, l1_sx_exc, s_inh_exc_l5l1_sx)
        num_conn_l5_sx[:, 1] = num_conn_e + num_conn_i
        # Layer 5 to Layer 2
        s_exc_all_l5l2_sx, s_inh_exc_l5l2_sx = conneting_layers(l5_sx, l5_sx_exc, l5_sx_inh,
                                                          l2_sx, l2_sx_exc, l2_sx_inh, pMAXext)
        WeightMatrix_s_exc_all_l5l2_sx, num_conn_e = createWeightMatrix(l5_sx_exc, l2_sx, s_exc_all_l5l2_sx)
        WeightMatrix_s_inh_exc_l5l2_sx, num_conn_i = createWeightMatrix(l5_sx_inh, l2_sx_exc, s_inh_exc_l5l2_sx)
        num_conn_l5_sx[:, 2] = num_conn_e + num_conn_i
        # Layer 5 to Layer 3
        s_exc_all_l5l3_sx, s_inh_exc_l5l3_sx = conneting_layers(l5_sx, l5_sx_exc, l5_sx_inh,
                                                          l3_sx, l3_sx_exc, l3_sx_inh, pMAXext)
        WeightMatrix_s_exc_all_l5l3_sx, num_conn_e = createWeightMatrix(l5_sx_exc, l3_sx, s_exc_all_l5l3_sx)
        WeightMatrix_s_inh_exc_l5l3_sx, num_conn_i = createWeightMatrix(l5_sx_inh, l3_sx_exc, s_inh_exc_l5l3_sx)
        num_conn_l5_sx[:, 3] = num_conn_e + num_conn_i
        # Layer 5 to Layer 4
        s_exc_all_l5l4_sx, s_inh_exc_l5l4_sx = conneting_layers(l5_sx, l5_sx_exc, l5_sx_inh,
                                                          l4_sx, l4_sx_exc, l4_sx_inh, pMAXext)
        WeightMatrix_s_exc_all_l5l4_sx, num_conn_e = createWeightMatrix(l5_sx_exc, l4_sx, s_exc_all_l5l4_sx)
        WeightMatrix_s_inh_exc_l5l4_sx, num_conn_i = createWeightMatrix(l5_sx_inh, l4_sx_exc, s_inh_exc_l5l4_sx)
        num_conn_l5_sx[:, 4] = num_conn_e + num_conn_i
    # DX
    s_exc_all_l5_dx, s_inh_exc_l5_dx = conneting_population(l5_dx, l5_dx_exc, l5_dx_inh, pMAX)
    WeightMatrix_s_exc_all_l5_dx, num_conn_e = createWeightMatrix(l5_dx_exc, l5_dx, s_exc_all_l5_dx)
    WeightMatrix_s_inh_exc_l5_dx, num_conn_i = createWeightMatrix(l5_dx_inh, l5_dx_exc, s_inh_exc_l5_dx)
    num_conn_l5_dx[:, 5] = num_conn_e + num_conn_i
    if upward_conn:
        # Layer 0 to Layer 5
        s_exc_all_l0l5_dx, s_inh_exc_l0l5_dx = conneting_layers(l0_dx, l0_dx_exc, l0_dx_inh,
                                                                l5_dx, l5_dx_exc, l5_dx_inh, pMAXext)
        WeightMatrix_s_exc_all_l0l5_dx, num_conn_e = createWeightMatrix(l0_dx_exc, l5_dx, s_exc_all_l0l5_dx)
        WeightMatrix_s_inh_exc_l0l5_dx, num_conn_i = createWeightMatrix(l0_dx_inh, l5_dx_exc, s_inh_exc_l0l5_dx)
        num_conn_l0_dx[:, 5] = num_conn_e + num_conn_i
        # Layer 1 to Layer 5
        s_exc_all_l1l5_dx, s_inh_exc_l1l5_dx = conneting_layers(l1_dx, l1_dx_exc, l1_dx_inh,
                                                                l5_dx, l5_dx_exc, l5_dx_inh, pMAXext)
        WeightMatrix_s_exc_all_l1l5_dx, num_conn_e = createWeightMatrix(l1_dx_exc, l5_dx, s_exc_all_l1l5_dx)
        WeightMatrix_s_inh_exc_l1l5_dx, num_conn_i = createWeightMatrix(l1_dx_inh, l5_dx_exc, s_inh_exc_l1l5_dx)
        num_conn_l1_dx[:, 5] = num_conn_e + num_conn_i
        # Layer 2 to Layer 5
        s_exc_all_l2l5_dx, s_inh_exc_l2l5_dx = conneting_layers(l2_dx, l2_dx_exc, l2_dx_inh,
                                                                l5_dx, l5_dx_exc, l5_dx_inh, pMAXext)
        WeightMatrix_s_exc_all_l2l5_dx, num_conn_e = createWeightMatrix(l2_dx_exc, l5_dx, s_exc_all_l2l5_dx)
        WeightMatrix_s_inh_exc_l2l5_dx, num_conn_i = createWeightMatrix(l2_dx_inh, l5_dx_exc, s_inh_exc_l2l5_dx)
        num_conn_l2_dx[:, 5] = num_conn_e + num_conn_i
        # Layer 3 to Layer 5
        s_exc_all_l3l5_dx, s_inh_exc_l3l5_dx = conneting_layers(l3_dx, l3_dx_exc, l3_dx_inh,
                                                                l5_dx, l5_dx_exc, l5_dx_inh, pMAXext)
        WeightMatrix_s_exc_all_l3l5_dx, num_conn_e = createWeightMatrix(l3_dx_exc, l5_dx, s_exc_all_l3l5_dx)
        WeightMatrix_s_inh_exc_l3l5_dx, num_conn_i = createWeightMatrix(l3_dx_inh, l5_dx_exc, s_inh_exc_l3l5_dx)
        num_conn_l3_dx[:, 5] = num_conn_e + num_conn_i
        # Layer 4 to Layer 5
        s_exc_all_l4l5_dx, s_inh_exc_l4l5_dx = conneting_layers(l4, l4_dx_exc, l4_dx_inh,
                                                                l5_dx, l5_dx_exc, l5_dx_inh, pMAXext)
        WeightMatrix_s_exc_all_l4l5_dx, num_conn_e = createWeightMatrix(l4_dx_exc, l5_dx, s_exc_all_l4l5_dx)
        WeightMatrix_s_inh_exc_l4l5_dx, num_conn_i = createWeightMatrix(l4_dx_inh, l5_dx_exc, s_inh_exc_l4l5_dx)
        num_conn_l4_dx[:, 5] = num_conn_e + num_conn_i
    if downward_conn:
        # Layer 5 to Layer 0
        s_exc_all_l5l0_dx, s_inh_exc_l5l0_dx = conneting_layers(l5_dx, l5_dx_exc, l5_dx_inh,
                                                                l0_dx, l0_dx_exc, l0_dx_inh, pMAXext)
        WeightMatrix_s_exc_all_l5l0_dx, num_conn_e = createWeightMatrix(l5_dx_exc, l0_dx, s_exc_all_l5l0_dx)
        WeightMatrix_s_inh_exc_l5l0_dx, num_conn_i = createWeightMatrix(l5_dx_inh, l0_dx_exc, s_inh_exc_l5l0_dx)
        num_conn_l5_dx[:, 0] = num_conn_e + num_conn_i
        # Layer 5 to Layer 1
        s_exc_all_l5l1_dx, s_inh_exc_l5l1_dx = conneting_layers(l5_dx, l5_dx_exc, l5_dx_inh,
                                                                l1_dx, l1_dx_exc, l1_dx_inh, pMAXext)
        WeightMatrix_s_exc_all_l5l1_dx, num_conn_e = createWeightMatrix(l5_dx_exc, l1_dx, s_exc_all_l5l1_dx)
        WeightMatrix_s_inh_exc_l5l1_dx, num_conn_i = createWeightMatrix(l5_dx_inh, l1_dx_exc, s_inh_exc_l5l1_dx)
        num_conn_l5_dx[:, 1] = num_conn_e + num_conn_i
        # Layer 5 to Layer 2
        s_exc_all_l5l2_dx, s_inh_exc_l5l2_dx = conneting_layers(l5_dx, l5_dx_exc, l5_dx_inh,
                                                                l2_dx, l2_dx_exc, l2_dx_inh, pMAXext)
        WeightMatrix_s_exc_all_l5l2_dx, num_conn_e = createWeightMatrix(l5_dx_exc, l2_dx, s_exc_all_l5l2_dx)
        WeightMatrix_s_inh_exc_l5l2_dx, num_conn_i = createWeightMatrix(l5_dx_inh, l2_dx_exc, s_inh_exc_l5l2_dx)
        num_conn_l5_dx[:, 2] = num_conn_e + num_conn_i
        # Layer 5 to Layer 3
        s_exc_all_l5l3_dx, s_inh_exc_l5l3_dx = conneting_layers(l5_dx, l5_dx_exc, l5_dx_inh,
                                                                l3_dx, l3_dx_exc, l3_dx_inh, pMAXext)
        WeightMatrix_s_exc_all_l5l3_dx, num_conn_e = createWeightMatrix(l5_dx_exc, l3_dx, s_exc_all_l5l3_dx)
        WeightMatrix_s_inh_exc_l5l3_dx, num_conn_i = createWeightMatrix(l5_dx_inh, l3_dx_exc, s_inh_exc_l5l3_dx)
        num_conn_l5_dx[:, 3] = num_conn_e + num_conn_i
        # Layer 5 to Layer 4
        s_exc_all_l5l4_dx, s_inh_exc_l5l4_dx = conneting_layers(l5_dx, l5_dx_exc, l5_dx_inh,
                                                                l4_dx, l4_dx_exc, l4_dx_inh, pMAXext)
        WeightMatrix_s_exc_all_l5l4_dx, num_conn_e = createWeightMatrix(l5_dx_exc, l4_dx, s_exc_all_l5l4_dx)
        WeightMatrix_s_inh_exc_l5l4_dx, num_conn_i = createWeightMatrix(l5_dx_inh, l4_dx_exc, s_inh_exc_l5l4_dx)
        num_conn_l5_dx[:, 4] = num_conn_e + num_conn_i
#####################################


#####################################
if __name__ == '__main__':
    print("Lif Connections - END")
    # print(len(s_exc_all_l0.N_outgoing_pre), len(s_inh_exc_l0.N_outgoing_pre[n_exc:]))
    # riga_l0 = WeightMatrix_s_exc_all_l0[0, :]
    # print(np.count_nonzero(riga_l0)-np.count_nonzero(np.isnan(riga_l0)))
    #
    # for ii in range(WeightMatrix_s_exc_all_l0.shape[0]):
    #     print(np.count_nonzero(WeightMatrix_s_exc_all_l0[ii, :])

    fig1 = figure(figsize=[5 * cm_fig_conv, 5 * cm_fig_conv], dpi=600.0)
    neuron_idx_exc = np.random.randint(0, n_exc)
    plot(l0_sx_exc.x[neuron_idx_exc] / umeter, l0_sx_exc.y[neuron_idx_exc] / umeter, 'o', mec='g', mfc='g')
    scatter(l0_sx.x[s_exc_all_l0_sx.j[neuron_idx_exc, :]] / umeter,
            l0_sx.y[s_exc_all_l0_sx.j[neuron_idx_exc, :]] / umeter,
            s=((s_exc_all_l0_sx.w[neuron_idx_exc, :] / mV) * 5) ** 6, marker='.', c='r')
    # scatter(l0_dx.x[s_exc_all_l0_sxdx.j[neuron_idx_exc, :]] / umeter,
    #         l0_dx.y[s_exc_all_l0_sxdx.j[neuron_idx_exc, :]] / umeter,
    #         s=((s_exc_all_l0_sxdx.w[neuron_idx_exc, :] / mV) * 5) ** 6, marker='.', c='r')
    title('exc')
    xlabel('x (µm)')
    ylabel('y (µm)', rotation='vertical')
    gca().spines['top'].set_visible(False)
    gca().spines['right'].set_visible(False)
    axis('equal')
    tight_layout()
    show()

    fig2 = figure(2, figsize=[5 * cm_fig_conv, 5 * cm_fig_conv], dpi=600.0)
    neuron_idx_inh = np.random.randint(0, n_inh)
    plot(l0_dx_inh.x[neuron_idx_inh] / umeter, l0_dx_inh.y[neuron_idx_inh] / umeter, 'o', mec='g', mfc='g')
    # scatter(l0_dx_exc.x[s_inh_exc_l0_dx.j[neuron_idx_inh, :]] / umeter,
    #         l0_dx_exc.y[s_inh_exc_l0_dx.j[neuron_idx_inh, :]] / umeter,
    #         s=((s_inh_exc_l0_dx.w[neuron_idx_inh, :] / mV) * 5) ** 6, marker='.', c='b')
    scatter(l0_sx_exc.x[s_inh_exc_l0_dxsx.j[neuron_idx_inh, :]] / umeter,
            l0_sx_exc.y[s_inh_exc_l0_dxsx.j[neuron_idx_inh, :]] / umeter,
            s=((s_inh_exc_l0_dxsx.w[neuron_idx_inh, :] / mV) * 100) ** 100, marker='.', c='b')


    title('inh')
    xlabel('x (µm)')
    ylabel('y (µm)', rotation='vertical')
    gca().spines['top'].set_visible(False)
    gca().spines['right'].set_visible(False)
    axis('equal')
    tight_layout()
    show()

    fig3 = figure(3, figsize=[10*cm_fig_conv, 8*cm_fig_conv], dpi=600.0)
    ax = axes(projection='3d')
    neuron_idx_exc = np.random.randint(0, n_exc)
    ax.plot(l0_sx_exc.x[neuron_idx_exc] / umeter,
            l0_sx_exc.y[neuron_idx_exc] / umeter,
            l0_sx_exc.z[neuron_idx_exc] / umeter,
            'o', mec='g', mfc='g')
    ax.scatter(l0_sx.x[s_exc_all_l0_sx.j[neuron_idx_exc, :]] / umeter,
               l0_sx.y[s_exc_all_l0_sx.j[neuron_idx_exc, :]] / umeter,
               l0_sx.z[s_exc_all_l0_sx.j[neuron_idx_exc, :]] / umeter,
               s=((s_exc_all_l0_sx.w[neuron_idx_exc, :] / mV) * 5) ** 6, marker='.', c='r')
    ax.scatter(l1_sx.x[s_exc_all_l0l1_sx.j[neuron_idx_exc, :]] / umeter,
               l1_sx.y[s_exc_all_l0l1_sx.j[neuron_idx_exc, :]] / umeter,
               l1_sx.z[s_exc_all_l0l1_sx.j[neuron_idx_exc, :]] / umeter)
               # s=((s_exc_all_l0l1_sx.w[neuron_idx_exc, :] / mV) * 5) ** 6, marker='.', c='r')
    # scatter(l0_dx.x[s_exc_all_l0_sxdx.j[neuron_idx_exc, :]] / umeter,
    #         l0_dx.y[s_exc_all_l0_sxdx.j[neuron_idx_exc, :]] / umeter,
    #         s=((s_exc_all_l0_sxdx.w[neuron_idx_exc, :] / mV) * 5) ** 6, marker='.', c='r')
    title('exc')
    ax.set_xlabel('x (µm)')
    ax.set_ylabel('y (µm)')
    ax.set_zlabel('z (µm)')
    axis('auto')
    tight_layout()
    show()

    # fig4 = figure(4, figsize=[5 * cm_fig_conv, 5 * cm_fig_conv], dpi=600.0)
    # neuron_idx_inh = np.random.randint(0, n_inh)
    # scatter(l0_dx_inh.x[neuron_idx_inh] / umeter, l0_dx_inh.y[neuron_idx_inh] / umeter, l0_dx_inh.z[neuron_idx_inh] / umeter,  'o', mec='g', mfc='g')
    # scatter(l0_dx_exc.x[s_inh_exc_l0_dx.j[neuron_idx_inh, :]] / umeter,
    #         l0_dx_exc.y[s_inh_exc_l0_dx.j[neuron_idx_inh, :]] / umeter,
    #         s=((s_inh_exc_l0_dx.w[neuron_idx_inh, :] / mV) * 5) ** 6, marker='.', c='b')
    #
    # title('inh')
    # xlabel('x (µm)')
    # ylabel('y (µm)', rotation='vertical')
    # gca().spines['top'].set_visible(False)
    # gca().spines['right'].set_visible(False)
    # axis('equal')
    # tight_layout()
    # show()

    # fig1.savefig('exc_8perc_connessioni_6_p', dpi=600, transparent=True)
    # fig2.savefig('inh_8perc_connessioni_6_p', dpi=600, transparent=True)
else:
    print("Finished Lif Connections")
#####################################
#####################################
