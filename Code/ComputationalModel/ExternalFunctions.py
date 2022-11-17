import os
import numpy as np
import matplotlib.pyplot as plt
import math
from brian2 import Synapses, umeter, ms, Unit, NeuronGroup


def crea_e_o_apri(x: str):
    actual_folder = os.getcwd()
    try:
        os.makedirs(actual_folder + "\\" + x)
    except FileExistsError:
        # print('The folder already exists\n')
        pass
    except FileNotFoundError:
        # print('The folder already exists\n')
        pass
    os.chdir(x)


def createWeightMatrix(pop1: NeuronGroup, pop2: NeuronGroup, conn: Synapses):
    """Creates the matrix with the synaptic weights. Rows=sources - Columns=targets"""
    # weight_matrix = np.zeros((len(pop1), len(pop2)))
    weight_matrix = np.full((len(pop1), len(pop2)), np.nan)
    weight_matrix[conn.i[:], conn.j[:]] = conn.w[:]

    num_conn = []
    for ii in range(weight_matrix.shape[0]):
        num_conn.append(int(np.count_nonzero(weight_matrix[ii, :]) - np.count_nonzero(np.isnan(weight_matrix[ii, :]))))

    return weight_matrix, num_conn


def savelookuptable(namefile: str, peaktrain: dict, neugroup: NeuronGroup):
    """first column : i (brian) = ID neuron\n
    second column : position id in spatial configuration within the population = ID position"""

    with open(namefile, "a+") as f:
        f.write('ID neuron\tID position\n')
        for key, value in peaktrain.items():
            f.write(str(key) + '\t' + str(neugroup.index[key]+1) + '\n')
        f.close()


def savepeaktrain(file: str, peaktrain: dict, neugroup: NeuronGroup, n_exc: int, n_inh: int, num_layer: int):
    """
    Saves the peak trains to txt files\n
    The names of the files are sorted on the position of the neuron\n
    Layer after layer
    """
    # nota: qui gli indici sono già cambiati in modo che siano: PYsx-INsx-PYdx-INd

    for key, value in peaktrain.items():
        if key < n_exc:
            nf = file + '_ex'
        else:
            nf = file + '_in'

        nn = (n_exc + n_inh)*num_layer + 1

        if (neugroup.index[key] + nn) < 10:
            namefile = nf + '_000' + str(neugroup.index[key] + nn)
        elif 10 <= (neugroup.index[key] + nn) < 100:
            namefile = nf + '_00' + str(neugroup.index[key] + nn)
        elif 100 <= (neugroup.index[key] + nn) < 1000:
            namefile = nf + '_0' + str(neugroup.index[key] + nn)
        else:
            namefile = nf + '_' + str(neugroup.index[key] + nn)
        with open(namefile, "a+") as f:
            line = str(value).lstrip('[')
            newline = line.rstrip('] s')
            f.write(newline + '\n')
            f.close()


def connecting_int(pop1: NeuronGroup, pop2: NeuronGroup, syn_effect: str, sigma: Unit, v_ax: float, p_max: float,
                   correcting_value):
    # def connecting_int(pop1: NeuronGroup, pop2: NeuronGroup, syn_effect: str, sigma: Unit, v_ax: float, p_max: float,
    #                    correcting_value, seed_val: int):



    """
    Connecting two populations within the compartment\n
    exc_pop ---     exc, pop, 'I_synEX+= alpha_ex * w', sigma, V_ax_P, pMAX, 1000 * mV\n
    inh_exc ---     inh, exc, 'I_synIN-= alpha_in * w', sigma, V_ax_I, pMAX, 1000 * mV\n

    :param pop1: NeuronGroup
    :param pop2: NeuronGroup
    :param syn_effect: str
    :param sigma: unit
    :param v_ax: float
    :param p_max: float
    :param correcting_value
    # :param seed_val seed value
    :return: syn: the generated synapses
    """
    # np.random.seed(seed_val)

    syn = Synapses(pop1, pop2, 'w : volt', on_pre=syn_effect)
    syn.connect(condition='i!=j',
                p='p_max*exp(-((x_pre - x_post)**2 + (y_pre - y_post)**2)/(2*sigma**2))')
    syn.delay = '(sqrt(((x_pre- x_post)/umeter)**2 + ((y_pre - y_post)/umeter)**2))/v_ax*ms'
    syn.w_ = 'correcting_value * (1/((sigma/umeter)*sqrt(2*pi)))*exp(-((x_pre - x_post)**2 + (y_pre - y_post)**2)/(' \
             '2*sigma**2)) '
    return syn


def connecting_3D(pop1: NeuronGroup, pop2: NeuronGroup, syn_effect: str, sigma_plane: Unit, sigma_3D: Unit,
                  v_ax: float, p_max: float, correcting_value):
    # def connecting_ext(pop1: NeuronGroup, pop2: NeuronGroup, syn_effect: str, sigma_plane: Unit, sigma_3D: Unit,
    #                    v_ax: float, p_max: float, correcting_value, seed_val: int):
    """
    Connecting two populations within the compartment\n
    exc_pop ---     exc, pop, 'I_synEX+=alpha_ex * w', sigma, sigma_3D, V_ax_P, pMAXext, 1000 * mV\n
    inh_exc ---     inh, exc, 'I_synIN-=alpha_in * w', sigma, sigma_3D, V_ax_I, pMAXext, 1000 * mV\n

    :param pop1: NeuronGroup
    :param pop2: NeuronGroup
    :param syn_effect: str
    :param sigma_plane: unit
    :param sigma_3D: unit
    :param v_ax: float
    :param p_max: float
    :param correcting_value
    # :param seed_val seed value
    :return: syn: the generated synapses
    """
    # np.random.seed(seed_val)

    syn = Synapses(pop1, pop2, 'w : volt', on_pre=syn_effect)
    syn.connect(condition='index_up_pre==True',
                p='p_max*exp(-((x_pre - x_post)**2/(2*sigma_plane**2) + (y_pre - y_post)**2/(2*sigma_plane**2) + \
                      (z_pre - z_post)**2/(2*sigma_3D**2)))')
    syn.delay = 'sqrt(((x_pre - x_post)/umeter)**2 + ((y_pre - y_post)/umeter)**2 +((z_pre - z_post)/umeter)**2 )   \
                 /v_ax*ms'
    syn.w_ = 'correcting_value * (1 / sqrt(2*pi*(sigma_plane/umeter)**2*(sigma_3D/umeter))) * \
              exp(-((x_pre - x_post) ** 2 / (2 * sigma_plane ** 2) + (y_pre - y_post) ** 2 / (2 * sigma_plane ** 2) + \
              (z_pre - z_post) ** 2 / (2 * sigma_3D ** 2)))'
    return syn


def connecting_ext(pop1: NeuronGroup, pop2: NeuronGroup, syn_effect: str, sigma_x: Unit, sigma_y: Unit,
                   v_ax: float, p_max: float, correcting_value: float):
    """
    Connecting two populations pertaining to different compartments (simplified)\n

    :param pop1: NeuronGroup
    :param pop2: NeuronGroup
    :param syn_effect: str
    :param sigma_x: unit
    :param sigma_y: unit
    :param v_ax: float
    :param p_max: float
    :param correcting_value: float
    :return: syn: the generated synapses
    """

    syn = Synapses(pop1, pop2, 'w : 1', on_pre=syn_effect)
    syn.connect(condition='abs((y_pre - y_post)/umeter) < max_dist_y',
                p='p_max*exp(-((x_pre - x_post)**2/(2*sigma_x**2) + (y_pre - y_post)**2/(2*sigma_y**2)))')
    syn.delay = 'sqrt(((x_pre - x_post)/umeter)**2 + ((y_pre - y_post)/umeter)**2)/v_ax*ms'
    syn.w_ = 'correcting_value * (1 / sqrt(2*pi*(sigma_x/umeter)*(sigma_y/umeter))) * \
              exp(-((x_pre - x_post) ** 2 / (2 * sigma_x ** 2) + (y_pre - y_post) ** 2 / (2 * sigma_y ** 2)))'
    return syn


def max_conn_scale_free(a: float, degree_min: int, degree_max: int, n_cells: int, n_points: int, seed_val: int):
    """
    Creation of a list of maximum connections values,one for each neuron, with a power law distribution,
    indication of a scale-free topology\n
    :parameter a power law slope (-0.6)
    :parameter degree_min lower limit in degree distribution (1)
    :parameter degree_max upper limit in degree distribution (100)
    :parameter n_cells number of cells per layer
    :parameter n_points number of degree values
    :parameter seed_val seed value
    :return: max_conn_int - list of number of maximum outgoing connections, one for each neuron
    """
    np.random.seed(seed_val)

    xmin = np.log10(degree_min)
    xmax = np.log10(degree_max)

    def custom_distribution(x, rumore):
        return 10 ** ((a + rumore) * np.log10(x))


    degree_values = np.logspace(xmin, xmax, n_points)
    probabilities = np.array([custom_distribution(x, np.random.uniform(0, 0.05)) for x in degree_values])  # distrete version
    probabilities *= 1 / np.sum(probabilities)  # normalisation to have probability

    # plt.plot(degree_values, probabilities, '.r')
    # plt.xscale('log')
    # plt.yscale('log')
    # plt.show()

    max_conn_int = []
    for ii in range(n_points):
        max_conn_int = max_conn_int + [int(degree_values[ii])] * math.ceil((probabilities[ii] * n_cells))
        # occurrences = probabilities[ii] * n_cells
        # print('num conn: ' + str(degree_values[ii]) + ' --> ' + str(int(degree_values[ii])) + '\t' + 'ripetizioni: ' + str(
        #     (probabilities[ii] * n_cells)) + ' --> ' + str(int((probabilities[ii] * n_cells))))

    # print(len(max_conn_int))
    # print(max_conn_int)
    np.random.shuffle(max_conn_int)
    # print(max_conn_int)
    if len(max_conn_int) > n_cells:
        del max_conn_int[-(len(max_conn_int) - n_cells):]
    elif len(max_conn_int) < n_cells:
        max_conn_int = max_conn_int + [int(degree_values[0])] * (n_cells - len(max_conn_int))
    # print(max_conn_int)
    return max_conn_int
