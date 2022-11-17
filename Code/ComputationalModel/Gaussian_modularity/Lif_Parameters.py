"""
PARAMETERS: Documentation
This module specifies all the parameters necessary for a neuronal papopulation made up of:
- Excitatory cells
- Inhibitory cells
"""
from brian2 import *

cm_fig_conv = 1/2.54
rc('font', size=6)
###########################################
#       User defined parameters
###########################################
val = 6
seed(val)
pMAX = 0.2
pMAXext = pMAX/3

num_layers = 4
perc_inh_cell = 20
perc_conn_int = 8       # percentuale di connessioni massime permesse all'interno della popolazione rispetto al numero di neuroni totale
perc_conn_ext_exc = 8       # CAMBIATA
perc_conn_ext_inh = 8       # CAMBIATA # percentuale di connessioni massime permesse all'interno della popolazione rispetto al numero di neuroni totale
upward_conn = True
downward_conn = True
###########################################


###########################################
#       Spatial Layout
###########################################
dist_pop = 50 * umeter
dist_modules = 250 * umeter
grid_dist = 30 * umeter
max_dist_y = (50*2) + (10*3)    # umeter, the distance of 3 uchannels

n_cells = 1091      # number of cells in each layer and each compartment
rows, cols = 33, 34
# rows, cols = 36, 28
# n_cells = rows * cols
x_trasl= (grid_dist * cols + dist_modules) / 2

n_inh = int(n_cells*perc_inh_cell/100)
n_exc = n_cells - n_inh

max_conn_int = int(n_cells*perc_conn_int/100) + np.random.randint(-10, 10, size=n_cells)
max_conn_ext_exc = int(n_cells*perc_conn_ext_exc/100) + np.random.randint(-10, 10, size=n_cells)    # CAMBIATA
max_conn_ext_inh = int(n_cells*perc_conn_ext_inh/100) + np.random.randint(-10, 10, size=n_cells)    # CAMBIATA

sigma = 1000*umeter
sigma_3D = sigma
sigma_2D_x = max_dist_y * umeter + dist_modules
sigma_2D_y = max_dist_y * umeter
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


