from Lif_SpatialConfiguration import *
from ExternalFunctions import connecting_int, createWeightMatrix, connecting_3D
import random

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

    s_exc_all = connecting_int(exc, pop, 'I_synEX+= alpha_ex * w', sigma, V_ax_P, pMAX, 1000 * mV, val)
    s_inh_exc = connecting_int(inh, exc, 'I_synIN-= alpha_in * w', sigma, V_ax_I, pMAX, 1000 * mV, val)
    
    for ss in [s_exc_all, s_inh_exc]:
        ll = len(ss.N_outgoing_pre)
        diff = 0
        count = 0
        k = 0

        for idx in range(ll):
            count = count + ss.N_outgoing_pre[idx]
            if ss.N_outgoing_pre[idx] > max_conn_int[idx]:
                diff = ss.N_outgoing_pre[idx] - max_conn_int[idx]   # connessioni da eliminare
                k = count - diff
                r = []
                r = random.sample(range(ss.N_outgoing_pre[idx]), diff)      # prendo tanti indici quanto le connessioni da mettere a peso zero
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

    s_exc_all = connecting_3D(exc_source, pop_target, 'I_synEX+=alpha_ex * w', sigma, sigma_3D, V_ax_P, p_max, 1000 * mV, val)
    s_inh_exc = connecting_3D(inh_source, exc_target, 'I_synIN-=alpha_in * w', sigma, sigma_3D, V_ax_I, p_max, 1000 * mV, val)

    for ss in [s_exc_all, s_inh_exc]:
        ll = len(ss.N_outgoing_pre)
        diff = 0
        count = 0
        k = 0

        for idx in range(ll):
            count = count + ss.N_outgoing_pre[idx]
            if ss.N_outgoing_pre[idx] > max_conn_ext[idx]:
                diff = ss.N_outgoing_pre[idx] - max_conn_ext[idx]  # connessioni da eliminare
                k = count - diff
                r = []
                r = random.sample(range(ss.N_outgoing_pre[idx]), diff)  # prendo tanti indici quanto le connessioni da mettere a peso zero
                for index1 in range(len(r)):
                    ss.w[r[index1] + k - max_conn_ext[idx]] = 0
    return s_exc_all, s_inh_exc


###########################################
#       CONNECTIVITY Implementation
###########################################
num_conn_l0 = np.zeros([n_cells, num_layers])
num_conn_l1 = np.zeros([n_cells, num_layers])
num_conn_l2 = np.zeros([n_cells, num_layers])
num_conn_l3 = np.zeros([n_cells, num_layers])
num_conn_l4 = np.zeros([n_cells, num_layers])
num_conn_l5 = np.zeros([n_cells, num_layers])

# Layer 0
s_exc_all_l0, s_inh_exc_l0 = conneting_population(l0, l0_exc, l0_inh, pMAX)
WeightMatrix_s_exc_all_l0, num_conn_e = createWeightMatrix(l0_exc, l0, s_exc_all_l0)
WeightMatrix_s_inh_exc_l0, num_conn_i = createWeightMatrix(l0_inh, l0_exc, s_inh_exc_l0)
num_conn_l0[:, 0] = num_conn_e + num_conn_i

# Layer 1
if num_layers > 1:
    s_exc_all_l1, s_inh_exc_l1 = conneting_population(l1, l1_exc, l1_inh, pMAX)
    WeightMatrix_s_exc_all_l1, num_conn_e = createWeightMatrix(l1_exc, l1, s_exc_all_l1)
    WeightMatrix_s_inh_exc_l1, num_conn_i = createWeightMatrix(l1_inh, l1_exc, s_inh_exc_l1)
    num_conn_l1[:, 1] = num_conn_e + num_conn_i
    if upward_conn:
        # Layer 0 to Layer 1
        s_exc_all_l0l1, s_inh_exc_l0l1 = conneting_layers(l0, l0_exc, l0_inh,
                                                          l1, l1_exc, l1_inh, pMAXext)
        WeightMatrix_s_exc_all_l0l1, num_conn_e = createWeightMatrix(l0_exc, l1, s_exc_all_l0l1)
        WeightMatrix_s_inh_exc_l0l1, num_conn_i = createWeightMatrix(l0_inh, l1_exc, s_inh_exc_l0l1)
        num_conn_l0[:, 1] = num_conn_e + num_conn_i
    if downward_conn:
        # Layer 1 to Layer 0
        s_exc_all_l1l0, s_inh_exc_l1l0 = conneting_layers(l1, l1_exc, l1_inh,
                                                          l0, l0_exc, l0_inh, pMAXext)
        WeightMatrix_s_exc_all_l1l0, num_conn_e = createWeightMatrix(l1_exc, l0, s_exc_all_l1l0)
        WeightMatrix_s_inh_exc_l1l0, num_conn_i = createWeightMatrix(l1_inh, l0_exc, s_inh_exc_l1l0)
        num_conn_l1[:, 0] = num_conn_e + num_conn_i

# Layer 2
if num_layers > 2:
    s_exc_all_l2, s_inh_exc_l2 = conneting_population(l2, l2_exc, l2_inh, pMAX)
    WeightMatrix_s_exc_all_l2, num_conn_e = createWeightMatrix(l2_exc, l2, s_exc_all_l2)
    WeightMatrix_s_inh_exc_l2, num_conn_i = createWeightMatrix(l2_inh, l2_exc, s_inh_exc_l2)
    num_conn_l2[:, 2] = num_conn_e + num_conn_i
    if upward_conn:
        # Layer 0 to Layer 2
        s_exc_all_l0l2, s_inh_exc_l0l2 = conneting_layers(l0, l0_exc, l0_inh,
                                                          l2, l2_exc, l2_inh, pMAXext)
        WeightMatrix_s_exc_all_l0l2, num_conn_e = createWeightMatrix(l0_exc, l2, s_exc_all_l0l2)
        WeightMatrix_s_inh_exc_l0l2, num_conn_i = createWeightMatrix(l0_inh, l2_exc, s_inh_exc_l0l2)
        num_conn_l0[:, 2] = num_conn_e + num_conn_i

        # Layer 1 to Layer 2
        s_exc_all_l1l2, s_inh_exc_l1l2 = conneting_layers(l1, l1_exc, l1_inh,
                                                          l2, l2_exc, l2_inh, pMAXext)
        WeightMatrix_s_exc_all_l1l2, num_conn_e = createWeightMatrix(l1_exc, l2, s_exc_all_l1l2)
        WeightMatrix_s_inh_exc_l1l2, num_conn_i = createWeightMatrix(l1_inh, l2_exc, s_inh_exc_l1l2)
        num_conn_l1[:, 2] = num_conn_e + num_conn_i
    if downward_conn:
        # Layer 2 to Layer 0
        s_exc_all_l2l0, s_inh_exc_l2l0 = conneting_layers(l2, l2_exc, l2_inh,
                                                          l0, l0_exc, l0_inh, pMAXext)
        WeightMatrix_s_exc_all_l2l0, num_conn_e = createWeightMatrix(l2_exc, l0, s_exc_all_l2l0)
        WeightMatrix_s_inh_exc_l2l0, num_conn_i = createWeightMatrix(l2_inh, l0_exc, s_inh_exc_l2l0)
        num_conn_l2[:, 0] = num_conn_e + num_conn_i
        # Layer 2 to Layer 1
        s_exc_all_l2l1, s_inh_exc_l2l1 = conneting_layers(l2, l2_exc, l2_inh,
                                                          l1, l1_exc, l1_inh, pMAXext)
        WeightMatrix_s_exc_all_l2l1, num_conn_e = createWeightMatrix(l2_exc, l1, s_exc_all_l2l1)
        WeightMatrix_s_inh_exc_l2l1, num_conn_i = createWeightMatrix(l2_inh, l1_exc, s_inh_exc_l2l1)
        num_conn_l2[:, 1] = num_conn_e + num_conn_i

# Layer 3
if num_layers > 3:
    s_exc_all_l3, s_inh_exc_l3 = conneting_population(l3, l3_exc, l3_inh, pMAX)
    WeightMatrix_s_exc_all_l3, num_conn_e = createWeightMatrix(l3_exc, l3, s_exc_all_l3)
    WeightMatrix_s_inh_exc_l3, num_conn_i = createWeightMatrix(l3_inh, l3_exc, s_inh_exc_l3)
    num_conn_l3[:, 3] = num_conn_e + num_conn_i
    if upward_conn:
        # Layer 0 to Layer 3
        s_exc_all_l0l3, s_inh_exc_l0l3 = conneting_layers(l0, l0_exc, l0_inh,
                                                          l3, l3_exc, l3_inh, pMAXext)
        WeightMatrix_s_exc_all_l0l3, num_conn_e = createWeightMatrix(l0_exc, l3, s_exc_all_l0l3)
        WeightMatrix_s_inh_exc_l0l3, num_conn_i = createWeightMatrix(l0_inh, l3_exc, s_inh_exc_l0l3)
        num_conn_l0[:, 3] = num_conn_e + num_conn_i
        # Layer 1 to Layer 3
        s_exc_all_l1l3, s_inh_exc_l1l3 = conneting_layers(l1, l1_exc, l1_inh,
                                                          l3, l3_exc, l3_inh, pMAXext)
        WeightMatrix_s_exc_all_l1l3, num_conn_e = createWeightMatrix(l1_exc, l3, s_exc_all_l1l3)
        WeightMatrix_s_inh_exc_l1l3, num_conn_i = createWeightMatrix(l1_inh, l3_exc, s_inh_exc_l1l3)
        num_conn_l1[:, 3] = num_conn_e + num_conn_i
        # Layer 2 to Layer 3
        s_exc_all_l2l3, s_inh_exc_l2l3 = conneting_layers(l2, l2_exc, l2_inh,
                                                          l3, l3_exc, l3_inh, pMAXext)
        WeightMatrix_s_exc_all_l2l3, num_conn_e = createWeightMatrix(l2_exc, l3, s_exc_all_l2l3)
        WeightMatrix_s_inh_exc_l2l3, num_conn_i = createWeightMatrix(l2_inh, l3_exc, s_inh_exc_l2l3)
        num_conn_l2[:, 3] = num_conn_e + num_conn_i
    if downward_conn:
        # Layer 3 to Layer 0
        s_exc_all_l3l0, s_inh_exc_l3l0 = conneting_layers(l3, l3_exc, l3_inh,
                                                          l0, l0_exc, l0_inh, pMAXext)
        WeightMatrix_s_exc_all_l3l0, num_conn_e = createWeightMatrix(l3_exc, l0, s_exc_all_l3l0)
        WeightMatrix_s_inh_exc_l3l0, num_conn_i = createWeightMatrix(l3_inh, l0_exc, s_inh_exc_l3l0)
        num_conn_l3[:, 0] = num_conn_e + num_conn_i
        # Layer 3 to Layer 1
        s_exc_all_l3l1, s_inh_exc_l3l1 = conneting_layers(l3, l3_exc, l3_inh,
                                                          l1, l1_exc, l1_inh, pMAXext)
        WeightMatrix_s_exc_all_l3l1, num_conn_e = createWeightMatrix(l3_exc, l1, s_exc_all_l3l1)
        WeightMatrix_s_inh_exc_l3l1, num_conn_i = createWeightMatrix(l3_inh, l1_exc, s_inh_exc_l3l1)
        num_conn_l3[:, 1] = num_conn_e + num_conn_i
        # Layer 3 to Layer 2
        s_exc_all_l3l2, s_inh_exc_l3l2 = conneting_layers(l3, l3_exc, l3_inh,
                                                          l2, l2_exc, l2_inh, pMAXext)
        WeightMatrix_s_exc_all_l3l2, num_conn_e = createWeightMatrix(l3_exc, l2, s_exc_all_l3l2)
        WeightMatrix_s_inh_exc_l3l2, num_conn_i = createWeightMatrix(l3_inh, l2_exc, s_inh_exc_l3l2)
        num_conn_l3[:, 2] = num_conn_e + num_conn_i

# Layer 4
if num_layers > 4:
    s_exc_all_l4, s_inh_exc_l4 = conneting_population(l4, l4_exc, l4_inh, pMAX)
    WeightMatrix_s_exc_all_l4, num_conn_e = createWeightMatrix(l4_exc, l4, s_exc_all_l4)
    WeightMatrix_s_inh_exc_l4, num_conn_i = createWeightMatrix(l4_inh, l4_exc, s_inh_exc_l4)
    num_conn_l4[:, 4] = num_conn_e + num_conn_i
    if upward_conn:
        # Layer 0 to Layer 4
        s_exc_all_l0l4, s_inh_exc_l0l4 = conneting_layers(l0, l0_exc, l0_inh,
                                                          l4, l4_exc, l4_inh, pMAXext)
        WeightMatrix_s_exc_all_l0l4, num_conn_e = createWeightMatrix(l0_exc, l4, s_exc_all_l0l4)
        WeightMatrix_s_inh_exc_l0l4, num_conn_i = createWeightMatrix(l0_inh, l4_exc, s_inh_exc_l0l4)
        num_conn_l0[:, 4] = num_conn_e + num_conn_i
        # Layer 1 to Layer 4
        s_exc_all_l1l4, s_inh_exc_l1l4 = conneting_layers(l1, l1_exc, l1_inh,
                                                          l4, l4_exc, l4_inh, pMAXext)
        WeightMatrix_s_exc_all_l1l4, num_conn_e = createWeightMatrix(l1_exc, l4, s_exc_all_l1l4)
        WeightMatrix_s_inh_exc_l1l4, num_conn_i = createWeightMatrix(l1_inh, l4_exc, s_inh_exc_l1l4)
        num_conn_l1[:, 4] = num_conn_e + num_conn_i
        # Layer 2 to Layer 4
        s_exc_all_l2l4, s_inh_exc_l2l4 = conneting_layers(l2, l2_exc, l2_inh,
                                                          l4, l4_exc, l4_inh, pMAXext)
        WeightMatrix_s_exc_all_l2l4, num_conn_e = createWeightMatrix(l2_exc, l4, s_exc_all_l2l4)
        WeightMatrix_s_inh_exc_l2l4, num_conn_i = createWeightMatrix(l2_inh, l4_exc, s_inh_exc_l2l4)
        num_conn_l2[:, 4] = num_conn_e + num_conn_i
        # Layer 3 to Layer 4
        s_exc_all_l3l4, s_inh_exc_l3l4 = conneting_layers(l3, l3_exc, l3_inh,
                                                          l4, l4_exc, l4_inh, pMAXext)
        WeightMatrix_s_exc_all_l3l4, num_conn_e = createWeightMatrix(l3_exc, l4, s_exc_all_l3l4)
        WeightMatrix_s_inh_exc_l3l4, num_conn_i = createWeightMatrix(l3_inh, l4_exc, s_inh_exc_l3l4)
        num_conn_l3[:, 4] = num_conn_e + num_conn_i
    if downward_conn:
        # Layer 4 to Layer 0
        s_exc_all_l4l0, s_inh_exc_l4l0 = conneting_layers(l4, l4_exc, l4_inh,
                                                          l0, l0_exc, l0_inh, pMAXext)
        WeightMatrix_s_exc_all_l4l0, num_conn_e = createWeightMatrix(l4_exc, l0, s_exc_all_l4l0)
        WeightMatrix_s_inh_exc_l4l0, num_conn_i = createWeightMatrix(l4_inh, l0_exc, s_inh_exc_l4l0)
        num_conn_l4[:, 0] = num_conn_e + num_conn_i
        # Layer 4 to Layer 1
        s_exc_all_l4l1, s_inh_exc_l4l1 = conneting_layers(l4, l4_exc, l4_inh,
                                                          l1, l1_exc, l1_inh, pMAXext)
        WeightMatrix_s_exc_all_l4l1, num_conn_e = createWeightMatrix(l4_exc, l1, s_exc_all_l4l1)
        WeightMatrix_s_inh_exc_l4l1, num_conn_i = createWeightMatrix(l4_inh, l1_exc, s_inh_exc_l4l1)
        num_conn_l4[:, 1] = num_conn_e + num_conn_i
        # Layer 4 to Layer 2
        s_exc_all_l4l2, s_inh_exc_l4l2 = conneting_layers(l4, l4_exc, l4_inh,
                                                          l2, l2_exc, l2_inh, pMAXext)
        WeightMatrix_s_exc_all_l4l2, num_conn_e = createWeightMatrix(l4_exc, l2, s_exc_all_l4l2)
        WeightMatrix_s_inh_exc_l4l2, num_conn_i = createWeightMatrix(l4_inh, l2_exc, s_inh_exc_l4l2)
        num_conn_l4[:, 2] = num_conn_e + num_conn_i
        # Layer 4 to Layer 3
        s_exc_all_l4l3, s_inh_exc_l4l3 = conneting_layers(l4, l4_exc, l4_inh,
                                                          l3, l3_exc, l3_inh, pMAXext)
        WeightMatrix_s_exc_all_l4l3, num_conn_e = createWeightMatrix(l4_exc, l3, s_exc_all_l4l3)
        WeightMatrix_s_inh_exc_l4l3, num_conn_i = createWeightMatrix(l4_inh, l3_exc, s_inh_exc_l4l3)
        num_conn_l4[:, 3] = num_conn_e + num_conn_i

# Layer 5
if num_layers > 5:
    s_exc_all_l5, s_inh_exc_l5 = conneting_population(l5, l5_exc, l5_inh, pMAX)
    WeightMatrix_s_exc_all_l5, num_conn_e = createWeightMatrix(l5_exc, l5, s_exc_all_l5)
    WeightMatrix_s_inh_exc_l5, num_conn_i = createWeightMatrix(l5_inh, l5_exc, s_inh_exc_l5)
    num_conn_l5[:, 5] = num_conn_e + num_conn_i
    if upward_conn:
        # Layer 0 to Layer 5
        s_exc_all_l0l5, s_inh_exc_l0l5 = conneting_layers(l0, l0_exc, l0_inh,
                                                          l5, l5_exc, l5_inh, pMAXext)
        WeightMatrix_s_exc_all_l0l5, num_conn_e = createWeightMatrix(l0_exc, l5, s_exc_all_l0l5)
        WeightMatrix_s_inh_exc_l0l5, num_conn_i = createWeightMatrix(l0_inh, l5_exc, s_inh_exc_l0l5)
        num_conn_l0[:, 5] = num_conn_e + num_conn_i
        # Layer 1 to Layer 5
        s_exc_all_l1l5, s_inh_exc_l1l5 = conneting_layers(l1, l1_exc, l1_inh,
                                                          l5, l5_exc, l5_inh, pMAXext)
        WeightMatrix_s_exc_all_l1l5, num_conn_e = createWeightMatrix(l1_exc, l5, s_exc_all_l1l5)
        WeightMatrix_s_inh_exc_l1l5, num_conn_i = createWeightMatrix(l1_inh, l5_exc, s_inh_exc_l1l5)
        num_conn_l1[:, 5] = num_conn_e + num_conn_i
        # Layer 2 to Layer 5
        s_exc_all_l2l5, s_inh_exc_l2l5 = conneting_layers(l2, l2_exc, l2_inh,
                                                          l5, l5_exc, l5_inh, pMAXext)
        WeightMatrix_s_exc_all_l2l5, num_conn_e = createWeightMatrix(l2_exc, l5, s_exc_all_l2l5)
        WeightMatrix_s_inh_exc_l2l5, num_conn_i = createWeightMatrix(l2_inh, l5_exc, s_inh_exc_l2l5)
        num_conn_l2[:, 5] = num_conn_e + num_conn_i
        # Layer 3 to Layer 5
        s_exc_all_l3l5, s_inh_exc_l3l5 = conneting_layers(l3, l3_exc, l3_inh,
                                                          l5, l5_exc, l5_inh, pMAXext)
        WeightMatrix_s_exc_all_l3l5, num_conn_e = createWeightMatrix(l3_exc, l5, s_exc_all_l3l5)
        WeightMatrix_s_inh_exc_l3l5, num_conn_i = createWeightMatrix(l3_inh, l5_exc, s_inh_exc_l3l5)
        num_conn_l3[:, 5] = num_conn_e + num_conn_i
        # Layer 4 to Layer 5
        s_exc_all_l4l5, s_inh_exc_l4l5 = conneting_layers(l4, l4_exc, l4_inh,
                                                          l5, l5_exc, l5_inh, pMAXext)
        WeightMatrix_s_exc_all_l4l5, num_conn_e = createWeightMatrix(l4_exc, l5, s_exc_all_l4l5)
        WeightMatrix_s_inh_exc_l4l5, num_conn_i = createWeightMatrix(l4_inh, l5_exc, s_inh_exc_l4l5)
        num_conn_l4[:, 5] = num_conn_e + num_conn_i
    if downward_conn:
        # Layer 5 to Layer 0
        s_exc_all_l5l0, s_inh_exc_l5l0 = conneting_layers(l5, l5_exc, l5_inh,
                                                          l0, l0_exc, l0_inh, pMAXext)
        WeightMatrix_s_exc_all_l5l0, num_conn_e = createWeightMatrix(l5_exc, l0, s_exc_all_l5l0)
        WeightMatrix_s_inh_exc_l5l0, num_conn_i = createWeightMatrix(l5_inh, l0_exc, s_inh_exc_l5l0)
        num_conn_l5[:, 0] = num_conn_e + num_conn_i
        # Layer 5 to Layer 1
        s_exc_all_l5l1, s_inh_exc_l5l1 = conneting_layers(l5, l5_exc, l5_inh,
                                                          l1, l1_exc, l1_inh, pMAXext)
        WeightMatrix_s_exc_all_l5l1, num_conn_e = createWeightMatrix(l5_exc, l1, s_exc_all_l5l1)
        WeightMatrix_s_inh_exc_l5l1, num_conn_i = createWeightMatrix(l5_inh, l1_exc, s_inh_exc_l5l1)
        num_conn_l5[:, 1] = num_conn_e + num_conn_i
        # Layer 5 to Layer 2
        s_exc_all_l5l2, s_inh_exc_l5l2 = conneting_layers(l5, l5_exc, l5_inh,
                                                          l2, l2_exc, l2_inh, pMAXext)
        WeightMatrix_s_exc_all_l5l2, num_conn_e = createWeightMatrix(l5_exc, l2, s_exc_all_l5l2)
        WeightMatrix_s_inh_exc_l5l2, num_conn_i = createWeightMatrix(l5_inh, l2_exc, s_inh_exc_l5l2)
        num_conn_l5[:, 2] = num_conn_e + num_conn_i
        # Layer 5 to Layer 3
        s_exc_all_l5l3, s_inh_exc_l5l3 = conneting_layers(l5, l5_exc, l5_inh,
                                                          l3, l3_exc, l3_inh, pMAXext)
        WeightMatrix_s_exc_all_l5l3, num_conn_e = createWeightMatrix(l5_exc, l3, s_exc_all_l5l3)
        WeightMatrix_s_inh_exc_l5l3, num_conn_i = createWeightMatrix(l5_inh, l3_exc, s_inh_exc_l5l3)
        num_conn_l5[:, 3] = num_conn_e + num_conn_i
        # Layer 5 to Layer 4
        s_exc_all_l5l4, s_inh_exc_l5l4 = conneting_layers(l5, l5_exc, l5_inh,
                                                          l4, l4_exc, l4_inh, pMAXext)
        WeightMatrix_s_exc_all_l5l4, num_conn_e = createWeightMatrix(l5_exc, l4, s_exc_all_l5l4)
        WeightMatrix_s_inh_exc_l5l4, num_conn_i = createWeightMatrix(l5_inh, l4_exc, s_inh_exc_l5l4)
        num_conn_l5[:, 4] = num_conn_e + num_conn_i
#####################################


#####################################
if __name__ == '__main__':
    figure(1)
    neuron_idx_exc = np.random.randint(0, n_exc) # da modificare !!!!!!!!!!!!!!!!!!!!!!
    plot(l0_exc.x[neuron_idx_exc] / umeter, l0_exc.y[neuron_idx_exc] / umeter, 'o', mec='g', mfc='none')
    scatter(l0.x[s_exc_all_l0.j[neuron_idx_exc, :]] / umeter, l0.y[s_exc_all_l0.j[neuron_idx_exc, :]] / umeter,
            s=((s_exc_all_l0.w[neuron_idx_exc, :] / mV) * 5)**7.5, marker='.', c='r')
    title('exc')
    xlabel('x')
    ylabel('y', rotation='horizontal')
    axis('equal')
    tight_layout()
    show()

    figure(2)
    neuron_idx_inh = np.random.randint(0, n_inh) # da modificare!!!!!!!!!!!!!!!!
    plot(l0_inh.x[neuron_idx_inh] / umeter, l0_inh.y[neuron_idx_inh] / umeter, 'o', mec='y', mfc='none')
    scatter(l0_exc.x[s_inh_exc_l0.j[neuron_idx_inh, :]] / umeter, l0_exc.y[s_inh_exc_l0.j[neuron_idx_inh, :]] / umeter,
            s=((s_inh_exc_l0.w[neuron_idx_inh, :] / mV) * 5)**7.5, marker='.', c='b')
    title('inh')
    xlabel('x')
    ylabel('y', rotation='horizontal')
    axis('equal')
    tight_layout()
    show()

    print("Lif Connections - END")
else:
    print("Finished Lif Connections")
#####################################
#####################################
