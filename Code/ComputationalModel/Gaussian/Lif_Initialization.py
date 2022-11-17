from Lif_Connections import *

#######################################
#       INITIALIZATION
#######################################
for pop_ll in layers:
    pop_ll.vm = El
    pop_ll.I_synEX = 0 * mV
    pop_ll.I_synIN = 0 * mV
    pop_ll.I_noise = 0 * mV
#######################################


#####################################
#####################################
if __name__ == '__main__':
    print("Initialization 1 Hp population - END")
else:
    print("Finished Initialization 1 Hp population")
#####################################
#####################################
