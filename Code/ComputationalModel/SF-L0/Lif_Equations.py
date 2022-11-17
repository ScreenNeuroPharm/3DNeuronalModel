"""
EQUATIONS: Documentation

"""
from brian2 import Equations, units

###########################################
# CELLS Equations
###########################################
eqs_LIF = Equations('''dvm/dt = ((El - vm) + I_synEX + I_synIN + I_noise)/tau : volt (unless refractory)
              I_noise : volt
              x : meter
              y : meter
              z : meter
              index : integer (constant)
          ''')
###########################################


###########################################
#   SYNAPTIC Equations
###########################################
eqs_LIF += Equations('''dI_synEX/dt = (-I_synEX)/tau : volt          # post-synaptic exc. 
                        dI_synIN/dt = (-I_synIN)/tau : volt          # post-synaptic inh. 
                    ''')
###########################################

###########################################
###########################################
if __name__ == '__main__':
    print("LIF_Equations - END")
else:
    print("Finished LIF_Equations")
###########################################
###########################################
