'''
Features the functions used to realize the dynamic programmin part of Mouysset & Jean, 2023.

Contains :
- data import functions for 3x3 landscapes and 4x4 landscapes : data_list2, data_list
- Unit level dynamic programming for each landscape, in 3x3 and 4x4 : lowlev3_dynprog, lowlev_dynprog
'''

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
    path = "C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/all_data/"
    files = os.listdir(path)

    filename = open(path+files[input], 'rb')
    list_nz = pickle.load(filename)
    filename.close()
    list_nz = list(list_nz)
    return list_nz

def data_list3(input):
    """
    Utility function for data import from 3x3 landscapes
    :param input: int
    :return: list
    """
    path = "C:/Users/jean/PycharmProjects/connectivity_dilemma/data3/all_data/"

    files = os.listdir(path)
    filename = open(path+files[input], 'rb')
    list_nz = pickle.load(filename)
    filename.close()
    if input != 0:
        list_nz = [item for sublist in list_nz for item in sublist]
    return list_nz
#endregion


def lowlev3_dynprog(land):
    """
    Function that computes the value function from Mouysset & Jean, 2023 on 3x3 landscapes, and returns the optimal landscape successions, the min
    value and the optimal prescribed burns

    :param land: np.array
    :return: list - lands, prescribed burns, values
    """
    # Load the environment
    #data_list3(input)
    land = tuple(land)
    # Set up the output variables
    # compute the new land and associated value
    new_land = utilities.fuel_dyn(land, params.pot_fire_budget)
    value = utilities.high_fuel_con_vect(new_land) + utilities.high_fuel_connect(land) + 1000000 * (utilities.high_biod_con_vect(new_land) <= params.biodiv3)
    # find min
    a = value.min(0)
    b = value.argmin(0)

    store = [tuple(x) for x in new_land[np.asarray(b, int)][0].tolist()]
    listed_values= [[land], store, np.asarray(b)[0].tolist(), np.asarray(a)[0].tolist()]
    list_values = [item for sublist in listed_values for item in sublist]

    return list_values

def lowlev_dynprog(land):
    """
    Function that computes the value function from Mouysset & Jean, 2023 on 4x4 landscapes, and returns the optimal landscape successions, the min
    value and the optimal prescribed burns

    :param land: np.array
    :return: list - lands, prescribed burns, values
    """
    # Load the environment
    #data_list2(input)
    land = tuple(land)
    # Set up the output variables
    # compute the new land and associated value
    new_land = utilities.fuel_dyn(land, params.pot_fire_budget)
    value = utilities.high_fuel_con_vect(new_land) + utilities.high_fuel_connect(land) + 1000000 * (utilities.high_biod_con_vect(new_land) <= params.biodiv)
    # find min
    a = value.min(0)
    b = value.argmin(0)

    store = [tuple(x) for x in new_land[np.asarray(b, int)][0].tolist()]
    listed_values = [[land], store, np.asarray(b)[0].tolist(), np.asarray(a)[0].tolist()]
    list_values = [item for sublist in listed_values for item in sublist]

    return list_values
