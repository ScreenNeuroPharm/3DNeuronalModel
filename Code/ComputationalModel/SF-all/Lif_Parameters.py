"""
PARAMETERS: Documentation
This module specifies all the parameters necessary for a neuronal papopulation made up of:
- Excitatory cells
- Inhibitory cells
"""
from brian2 import *
from ExternalFunctions import max_conn_scale_free

###########################################
#       User defined parameters
###########################################
val = 6
seed(val)
pMAX = 0.2
pMAXext = pMAX/3

num_layers = 1
perc_inh_cell = 20
perc_conn_int = 8       # percentuale di connessioni massime permesse all'interno della popolazione rispetto al numero di neuroni totale
perc_conn_ext = 8       # percentuale di connessioni massime permesse all'interno della popolazione rispetto al numero di neuroni totale
upward_conn = True
downward_conn = True
###########################################


###########################################
#       Spatial Layout
###########################################
dist_pop = 50 * umeter
grid_dist = 26 * umeter

n_cells = 1091
rows, cols = 33, 34
# rows, cols = 36, 28
# n_cells = rows * cols

n_inh = int(n_cells*perc_inh_cell/100)
n_exc = n_cells - n_inh

# max_conn_int = int(n_cells*perc_conn_int/100) + np.random.randint(-10, 10, size=n_cells)
a = -0.6
lower_degree_limit = 50
upper_degree_limit = 200
number_degree_values = 100
max_conn_int = max_conn_scale_free(a, lower_degree_limit, upper_degree_limit, n_cells, number_degree_values, val)
max_conn_ext = int(n_cells*perc_conn_ext/100) + np.random.randint(-10, 10, size=n_cells)

sigma = 1000*umeter
sigma_3D = sigma
###########################################


###########################################
#   LIF model variables
###########################################
tau = 10 * ms       # membrane time constant
v_r = -70 * mV      # reset potential
El = -70 * mV
v_th = -50 * mV     # threshold potential
refractory_period = 5*ms
###########################################


###########################################
#   HPC AXONAL DELAYS
###########################################
dx = 10e-3         # 10 microns = 10*10^-3 mm
V_ax_P = 0.5       # PY axonal conductance velocity (mm/msec)
# V_ax_I = 0.1       # IN axonal conductance velocity (mm/msec)
V_ax_I = 0.5       # IN axonal conductance velocity (mm/msec)
###########################################


###########################################
#   SYNAPTIC Parameters
###########################################
alpha_ex = 12
alpha_in = 16
###########################################


###########################################
###########################################
if __name__ == '__main__':
    print("Lif_Parameters - END")
else:
    print("Finished Lif_Parameters")
###########################################
###########################################


