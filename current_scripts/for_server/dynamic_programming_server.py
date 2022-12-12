import current_scripts.utilities as utilities
import current_scripts.params as params
#region Packages

import numpy as np
import os
import pickle

#endregion

#region Data import
def data_list2(input):
    """
    Utility function for data import from 4x4 landscapes

    :param input: int
    :return: list
    """
    path = "/home/simonjean/data/all_data/"
    files = os.listdir(path)

    filename = open(path+files[input], 'rb')
    list_nz = pickle.load(filename)
    filename.close()
    list_nz = list(list_nz)
    return list_nz

#endregion

def lowlev_dynprog_cut_mod(land):
    # Load the environment
    data_list2(input)
    land = tuple(land)
    # Set up the output variables
    # compute the new land and associated value
    new_land = utilities.fuel_dyn(land, params.pot_fire_budget)
    value = utilities.high_fuel_con_vect(new_land) + utilities.high_fuel_connect(land) + \
            1000000 * (utilities.high_biod_con_vect(new_land) <= params.biodiv)
    # find min
    a = value.min(0)
    b = value.argmin(0)

    store = [tuple(x) for x in new_land[np.asarray(b, int)][0].tolist()]
    listed_values = [[land], store, np.asarray(b)[0].tolist(), np.asarray(a)[0].tolist()]
    list_values = [item for sublist in listed_values for item in sublist]

    return list_values
