#region Packages
import statistics

import numpy as np
import pandas as pd
import time
import ast
import csv
import statistics as stats
import random as rd
import re as re
import math
import multiprocessing as mp
from playsound import playsound
import os
import json
import pickle
import itertools
import matplotlib.pyplot as plt
import seaborn as sns
#endregion

#region Define parameters
# region Parameters that can vary


#endregion
exec(open('current_scripts/params.py').read())

#endregion

#region Functions to generate landscapes


#region Landscape 3x3
def to_matrix(essai,size):
    """
    This function translates the array of length size into a sqrt(size)xsqrt(size) matrix
    """
    if size==3:
        return(np.array([list(essai[0:3]),
                         list(essai[3:6]),
                         list(essai[6:9])]))
    elif size==4:
        essai2 = np.array([list(essai[0:4]),
                           list(essai[4:8]),
                           list(essai[8:12]),
                           list(essai[12:16])])
        return (essai2)
    elif size==5:
        essai2 = np.array([list(essai[0:5]),
                           list(essai[5:10]),
                           list(essai[10:15]),
                           list(essai[15:20]),
                           list(essai[20:25])])
        return (essai2)
def to_matrix3(essai):
    """
    Transform vector landscape to matrix
    :param essai:
    :return: matrix of landscape
    """
    return(np.array([list(essai[0:3]),
              list(essai[3:6]),
              list(essai[6:9])]))
def equi_landscape3(essai):
    """
    This function computes the unique matrices equivalent to the current 3x3 landscape using :
    - the initial matrix and 3 consecutive 90° rotations
    - A matrix transformed using horizontal symmetry and 3 consecutive 90° rotations
    - A matrix transformed using vertical symmetry
    """
    # 3 transformations
    m2 = np.rot90(essai)
    m3 = np.rot90(m2)
    m4 = np.rot90(m3)
    # Horizontal symmetry
    hor_sym = [essai[2],essai[1],essai[0]]

    # 3 clockwise transformations
    m6 = np.rot90(hor_sym)
    m7 = np.rot90(m6)
    m8 = np.rot90(m7)

    evalu = (tuple(list(itertools.chain(*essai))),
            tuple(list(itertools.chain(*m2))),
            tuple(list(itertools.chain(*m3))),
            tuple(list(itertools.chain(*m4))),
            tuple(list(itertools.chain(*hor_sym))),
            tuple(list(itertools.chain(*m6))),
            tuple(list(itertools.chain(*m7))),
            tuple(list(itertools.chain(*m8))))
    return (set(evalu))
def low_lev_land3(nb2,non_z):
    """
    For a specified n_zero, decentralizes the computation inside unique_land4_nz
    """
    unique_landscapes_nb2 = set()
    size = 9
    non_zero = [1, 2]
    which = np.array(list(itertools.combinations(range(size), non_z)))
    grid = np.zeros((len(which), size), dtype="int8")
    grid[np.arange(len(which))[None].T, which] = 1

    for i in np.array(range(len(grid))):
        grid_arr = np.array(grid[i])
        index_replace = np.where(grid_arr == 1)
        candidat_lists = list(itertools.product(non_zero, repeat=non_z))

        nb2_list = [x for x in candidat_lists if sum(x) == (non_z-nb2)+nb2*2]
        # nb2_list = [x for x in candidat_lists if sum(x) == len(candidat_lists[1]) + nb2]
        # unique_landscapes_sub = list()
        # unique_landscapes_sub : unique landscapes for a given number of 2
        for j in nb2_list:
            # Want to keep only lists with a certain score so that we can compare further than non 0, use the number of 2
            candidate = grid_arr
            for k in np.array(range(len(j))):
                candidate[index_replace[0][k]] = j[k]

            mat_test = equi_landscape(to_matrix(candidate,3),3)

            if len(mat_test.intersection(unique_landscapes_nb2)) == 0:
                unique_landscapes_nb2.add(tuple(candidate))
    return(unique_landscapes_nb2)
#endregion
#region Landscape 4x4
def to_matrix4(essai):
    """
    This function translates the array of length 16 into a 4x4 matrix
    """
    essai2 = np.array([list(essai[0:4]),
                       list(essai[4:8]),
                       list(essai[8:12]),
                       list(essai[12:16])],dtype=object)
    return(essai2)
def equi_landscape4(essai):
    """
    This function computes the unique matrices equivalent to the current 4x4 landscape using :
    - the initial matrix and 3 consecutive 90° rotations
    - A matrix transformed using horizontal symmetry and 3 consecutive 90° rotations
    - A matrix transformed using vertical symmetry
    """
    # 3 transformations
    m2 = np.rot90(essai)
    m3 = np.rot90(m2)
    m4 = np.rot90(m3)
    # Horizontal symmetry
    hor_sym = [essai[3],essai[2],essai[1],essai[0]]

    # 3 clockwise transformations
    m6 = np.rot90(hor_sym)
    m7 = np.rot90(m6)
    m8 = np.rot90(m7)
    evalu = (tuple(list(itertools.chain(*m2))),
             tuple(list(itertools.chain(*m3))),
             tuple(list(itertools.chain(*m4))),
             tuple(list(itertools.chain(*hor_sym))),
             tuple(list(itertools.chain(*m6))),
             tuple(list(itertools.chain(*m7))),
             tuple(list(itertools.chain(*m8))))
    return (set(evalu))
def low_lev_land4(nb2,non_z):
    """
    Returns the unique landscapes from a combination of non zero elements and a number of twos.
    :param nb2: int
        Number of desired twos
    :param non_z: int
        Number of desired non zeros
    :return: set
        Set of unique landscapes of size 4x4

    """
    unique_landscapes_nb2 = set()
    size = 16
    non_zero = [1, 2]
    which = np.array(list(itertools.combinations(range(size), non_z)))
    grid = np.zeros((len(which), size), dtype="int8")
    grid[np.arange(len(which))[None].T, which] = 1

    for i in np.array(range(len(grid))):
        grid_arr = np.array(grid[i])
        index_replace = np.where(grid_arr == 1)
        candidat_lists = list(itertools.product(non_zero, repeat=non_z))

        nb2_list = [x for x in candidat_lists if sum(x) == (non_z-nb2)+nb2*2]
        # nb2_list = [x for x in candidat_lists if sum(x) == len(candidat_lists[1]) + nb2]
        # unique_landscapes_sub = list()
        # unique_landscapes_sub : unique landscapes for a given number of 2
        for j in nb2_list:
            # Want to keep only lists with a certain score so that we can compare further than non 0, use the number of 2
            candidate = grid_arr
            for k in np.array(range(len(j))):
                candidate[index_replace[0][k]] = j[k]

            mat_test = equi_landscape(to_matrix(candidate,4),4)

            if len(mat_test.intersection(unique_landscapes_nb2))== 0:
                unique_landscapes_nb2.add(tuple(candidate))
                #save
    #a_file = open("C:/Users/jean/PycharmProjects/connectivity_dilemma/data5/land5_" + str(non_z) + "_" + str(nb2) + ".pkl",
    #              "wb")
    #pickle.dump(results[k], a_file)
    #a_file.close()
    return(unique_landscapes_nb2)

#endregion
#region Landscape 5x5
def to_matrix5(essai):
    """
    This function translates the array of length 16 into a 4x4 matrix
    """
    essai2 = np.array([list(essai[0:5]),
                       list(essai[5:10]),
                       list(essai[10:15]),
                       list(essai[15:20]),
                       list(essai[20:25])],dtype=int)
    return(essai2)
def equi_landscape5(essai):
    """
    This function computes the unique matrices equivalent to the current 4x4 landscape using :
    - the initial matrix and 3 consecutive 90° rotations
    - A matrix transformed using horizontal symmetry and 3 consecutive 90° rotations
    - A matrix transformed using vertical symmetry
    """
    # 3 transformations
    m2 = np.rot90(essai)
    m3 = np.rot90(m2)
    m4 = np.rot90(m3)
    # Horizontal symmetry
    hor_sym = [essai[4],essai[3],essai[2],essai[1],essai[0]]

    # 3 clockwise transformations
    m6 = np.rot90(hor_sym)
    m7 = np.rot90(m6)
    m8 = np.rot90(m7)
    evalu = (tuple(list(itertools.chain(*m2))),
             tuple(list(itertools.chain(*m3))),
             tuple(list(itertools.chain(*m4))),
             tuple(list(itertools.chain(*hor_sym))),
             tuple(list(itertools.chain(*m6))),
             tuple(list(itertools.chain(*m7))),
             tuple(list(itertools.chain(*m8))))
    return (set(evalu))
#endregion
#endregion

#region Dynamics and objective functions
def connectivity_mat(size=3):
    """
    Returns the connectivity
    """
    if (size == 3):
        return (np.array([[1, 1, 0, 1, 1, 0, 0, 0, 0],
                          [1, 1, 1, 1, 1, 1, 0, 0, 0],
                          [0, 1, 1, 0, 1, 1, 0, 0, 0],
                          [1, 1, 0, 1, 1, 0, 1, 1, 0],
                          [1, 1, 1, 1, 1, 1, 1, 1, 1],
                          [0, 1, 1, 0, 1, 1, 0, 1, 1],
                          [0, 0, 0, 1, 1, 0, 1, 1, 0],
                          [0, 0, 0, 1, 1, 1, 1, 1, 1],
                          [0, 0, 0, 0, 1, 1, 0, 1, 1]])
                    )
    elif (size == 4):
        return (
                np.array([[1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                          [1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                          [0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                          [0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                          [1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0],
                          [1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0],
                          [0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0],
                          [0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0],
                          [0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0],
                          [0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0],
                          [0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1],
                          [0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1],
                          [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0],
                          [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0],
                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1],
                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1]])
            )
    elif size == 5:
        return (
                np.array([[1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # 1
                          [1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # 2
                          [0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # 3
                          [0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # 4
                          [0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # 5
                          [1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # 6
                          [1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # 7
                          [0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # 8
                          [0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # 9
                          [0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # 10
                          [0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],  # 11
                          [0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0],  # 12
                          [0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0],  # 13
                          [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0],  # 14
                          [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0],  # 15
                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0],  # 16
                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0],  # 17
                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0],  # 18
                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1],  # 19
                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1],  # 20
                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0],  # 21
                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0],  # 22
                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0],  # 23
                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1],  # 24
                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1]])
            )
    else:
        print("No connectivity matrix associated with this size")

def fuel_dyn(landscape, presc_burn):
    post = (np.array(landscape) + 1)*(1-np.array(presc_burn))
    land_post = np.minimum(post, 2)
    return(land_post)
def high_fuel(landscape):
    """
      This function returns a Boolean vector of fuel content larger than fire risk threshold for each cell.
    """
    return (np.array(landscape) >= d)
def high_fuel_connect(landscape,size=4):

    """
    Computes the fuel connectivity score for a given landscape.
    """
    return(high_fuel(landscape).dot(connectivity_mat(size)).dot(high_fuel(landscape)[:,np.newaxis]))
def high_biodiv(landscape):
    """
    This function returns a Boolean vector of fuel content larger than the biodiversity value threshold for each cell.
    """
    return (np.array(landscape) >= m_seuil)
def high_fuel_con_vect(land, size=4):
    return(np.transpose(np.matrix(np.diagonal((land>=d).dot(connectivity_mat(size)).dot(np.transpose((land>=d)))))))
def high_biod_con_vect(new_land, size=4):
   a= np.mat(np.diagonal((new_land>=m).dot(connectivity_mat(size)).dot(np.transpose((new_land>=m))))).T
   return(a)

def connectivity_mat3(size=3):
    """
    Returns the connectivity
    """
    if (size == 3):
        return (np.array([[1, 1, 0, 1, 1, 0, 0, 0, 0],
                          [1, 1, 1, 1, 1, 1, 0, 0, 0],
                          [0, 1, 1, 0, 1, 1, 0, 0, 0],
                          [1, 1, 0, 1, 1, 0, 1, 1, 0],
                          [1, 1, 1, 1, 1, 1, 1, 1, 1],
                          [0, 1, 1, 0, 1, 1, 0, 1, 1],
                          [0, 0, 0, 1, 1, 0, 1, 1, 0],
                          [0, 0, 0, 1, 1, 1, 1, 1, 1],
                          [0, 0, 0, 0, 1, 1, 0, 1, 1]])
                    )
    elif (size == 4):
        return (
                np.array([[1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                          [1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                          [0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                          [0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                          [1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0],
                          [1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0],
                          [0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0],
                          [0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0],
                          [0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0],
                          [0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0],
                          [0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1],
                          [0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1],
                          [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0],
                          [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0],
                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1],
                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1]])
            )
    elif size == 5:
        return (
                np.array([[1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # 1
                          [1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # 2
                          [0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # 3
                          [0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # 4
                          [0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # 5
                          [1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # 6
                          [1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # 7
                          [0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # 8
                          [0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # 9
                          [0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # 10
                          [0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],  # 11
                          [0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0],  # 12
                          [0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0],  # 13
                          [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0],  # 14
                          [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0],  # 15
                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0],  # 16
                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0],  # 17
                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0],  # 18
                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1],  # 19
                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1],  # 20
                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0],  # 21
                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0],  # 22
                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0],  # 23
                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1],  # 24
                          [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1]])
            )
    else:
        print("No connectivity matrix associated with this size")
def high_fuel_con_vect3(land, size=3):
    return(np.transpose(np.matrix(np.diagonal((land>=d).dot(connectivity_mat(size)).dot(np.transpose((land>=d)))))))
def high_fuel_connect3(landscape,size=3):

    """
    Computes the fuel connectivity score for a given landscape.
    """
    return(high_fuel(landscape).dot(connectivity_mat(size)).dot(high_fuel(landscape)[:,np.newaxis]))
def high_biod_con_vect3(new_land, size=3):
   a= np.mat(np.diagonal((new_land>=m).dot(connectivity_mat(size)).dot(np.transpose((new_land>=m))))).T
   return(a)
#endregion

#region Functions for dynamic programming
#region 3x3 landscape
def data_list3(input):
    path = "C:/Users/jean/PycharmProjects/connectivity_dilemma/data3/all_data/"

    files = os.listdir(path)
    filename = open(path+files[input], 'rb')
    list_nz = pickle.load(filename)
    filename.close()
    if input!=0:
        list_nz = [item for sublist in list_nz for item in sublist]
    return (list_nz)
def lowlev3_dynprog_cut_mod(land):
    # Load the environment
    #data_list3(input)
    land = tuple(land)
    # Set up the output variables
    # compute the new land and associated value
    new_land = fuel_dyn(land, pot_fire_budget)
    value = high_fuel_con_vect3(new_land) + high_fuel_connect3(land) + 1000000 * (high_biod_con_vect3(new_land) <= biodiv3)
    # find min
    a = value.min(0)
    b = value.argmin(0)

    store = [ tuple(x) for x in new_land[np.asarray(b, int)][0].tolist()]
    listed_values= [[land],store,np.asarray(b)[0].tolist(), np.asarray(a)[0].tolist() ]
    list_values = [item for sublist in listed_values for item in sublist]

    return list_values
#endregion

#region 4x4 landscape
def data_list2(input):
    path = "C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/all_data/"
    files = os.listdir(path)

    filename = open(path+files[input], 'rb')
    list_nz = pickle.load(filename)
    filename.close()
    list_nz = list(list_nz)
    return (list_nz)
def lowlev_dynprog_cut_mod(land, input):
    # Load the environment
    data_list2(input)
    land = tuple(land)
    # Set up the output variables
    # compute the new land and associated value
    new_land = fuel_dyn(land, pot_fire_budget)
    value = high_fuel_con_vect(new_land) + high_fuel_connect(land) + 1000000 * (high_biod_con_vect(new_land) <= biodiv)
    # find min
    a = value.min(0)
    b = value.argmin(0)

    store = [ tuple(x) for x in new_land[np.asarray(b, int)][0].tolist()]
    listed_values= [[land],store,np.asarray(b)[0].tolist(), np.asarray(a)[0].tolist() ]
    list_values = [item for sublist in listed_values for item in sublist]

    return list_values

#endregion
#endregion

#region Functions for matching
def nonzero(land):
    return(sum(land>0))
def nb2(land):
    return(sum(land==2))
def match_heavy_b(k):
    file = open(params.path+files[k])
    nums = re.findall(r'\d+',files[k])
    csvreader=csv.reader(file)
    header=next(csvreader)

    rows=[]
    for row in csvreader:
        row_add = row[1:42]
        rows.append(row_add)
    landscapes = [item for sublist in rows for item in sublist]

    singular_land = [np.array(ast.literal_eval(i),int) for i in list(set(landscapes))]


    path2 = "C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/all_data"
    files_data = os.listdir(path2)

    twos = set([matcher.nb2(i) for i in singular_land])
    nonzeros = set([matcher.nonzero(i) for i in singular_land])

    possible = list(itertools.product(nonzeros,twos))

    store = {}
    for i in range(len(possible)):
        store[possible[i]]=[j for j in singular_land if (nonzero(j)==possible[i][0] and nb2(j)==possible[i][1])]
    store2 = {k: v for k, v in store.items() if v}
    del(store)


#region Matching & replacing procedure
    debut = time.time()
    replacement = {}

    for j in list(store2.keys()):
        debut2 = time.time()
        if os.path.exists(path2+"/land4_" + str(j[0]) + "_" + str(j[1]) + ".pkl"):
            filename = open("C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/all_data/land4_" + str(j[0]) + "_" + str(
                j[1]) + ".pkl", 'rb')
            set_to_check = pickle.load(filename)
            filename.close()
            value = []
            to_rep=[]
            for v in store2[j]:
                set_verif = utilities.equi_landscape(utilities.to_matrix(v,4),4)
                set_verif.add(tuple(v))
                list_verif = [*set_verif.intersection(set_to_check)]
            #new = list(list_verif[0])
                replacement[str(tuple(v))] = str(list_verif[0])

        elif os.path.exists(path2+"/land4_" + str(j[0]) + "_" + str(j[1]) + "cut_0.pkl"):
            set_to_check_merg = []
            to_load = [x for x in files_data if x.startswith('land4_' + str(j[0]) + "_" + str(j[1]))]

            for x in to_load:
                filename = open("C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/all_data/"+x, 'rb')
                set_to_check = pickle.load(filename)
                filename.close()
                set_to_check =list(map(tuple, set_to_check))
                set_to_check_merg.append(set_to_check)

            set_to_check_merg = [item for sublist in set_to_check_merg for item in sublist]
            set_to_check_merg=set(set_to_check_merg)

            for v in store2[j]:
                set_verif = utilities.equi_landscape(utilities.to_matrix(v,4),4)
                set_verif.add(tuple(v))
                list_verif = [*set_verif.intersection(set_to_check_merg)]
            # new = list(list_verif[0])

                replacement[str(tuple(v))] = str(list_verif[0])

            del(set_to_check)
        else:
            print('ERROR : no file')
        #print("Step " +str(j)+" took "+str(time.time()-debut2) +" seconds")


    for i in list(range(len(rows))):
        rows[i] = [replacement.get(n, n) for n in rows[i]]

    tester3 = pd.DataFrame(rows,dtype=str)

    if len(nums) > 4:
        tester3.to_csv(path + "/keys_succession/equiv_results_" + nums[2] + "_" + nums[3] + "_cut_" + nums[4] + '.csv',
                  index=False)
    else:
        tester3.to_csv(path + "/keys_succession/equiv_results_" + nums[2] + "_" + nums[3] + ".csv", index=False)

    return(print('Dataset ' + str(files[k]) +" took " + str(time.time()-debut) + ' seconds'))
def match_heavy_3(k):
    file = open(path + files[k])
    nums = re.findall(r'\d+', files[k])
    csvreader = csv.reader(file)
    header = next(csvreader)

    rows = []
    for row in csvreader:
        row_add = row[0:26]
        rows.append(row_add)
    landscapes = [item for sublist in rows for item in sublist]

    singular_land = [np.array(ast.literal_eval(i), int) for i in list(set(landscapes))]
    # endregion

    # Make sublists

    # def nonzero(land):
    #    return(sum(land>0))
    # def nb2(land):
    #    return(sum(land==2))

    path2 = "C:/Users/jean/PycharmProjects/connectivity_dilemma/data3/all_data"
    files_data = os.listdir(path2)

    twos = set([nb2(i) for i in singular_land])
    nonzeros = set([nonzero(i) for i in singular_land])

    possible = list(itertools.product(nonzeros, twos))

    store = {}
    for i in range(len(possible)):
        store[possible[i]] = [j for j in singular_land if (nonzero(j) == possible[i][0] and nb2(j) == possible[i][1])]
    store2 = {k: v for k, v in store.items() if v}
    del (store)

    # region Matching & replacing procedure
    debut = time.time()
    replacement = {}

    for j in list(store2.keys()):
        debut2 = time.time()
        if os.path.exists(path2 + "/land3_" + str(j[0]) + ".pkl"):
            filename = open(
                "C:/Users/jean/PycharmProjects/connectivity_dilemma/data3/all_data/land3_" + str(j[0]) + ".pkl", 'rb')
            set_to_check2 = pickle.load(filename)
            filename.close()
            if k == 0 and j == (0, 0):
                set_to_check = set_to_check2
            else:
                set_to_check = set()
                for subset in set_to_check2:
                    set_to_check.update(subset)
            value = []
            to_rep = []
            for v in store2[j]:
                set_verif = equi_landscape3(to_matrix3(v))
                set_verif.add(tuple(v))
                list_verif = [*set_verif.intersection(set_to_check)]
                # new = list(list_verif[0])
                replacement[str(tuple(v))] = str(list_verif[0])

    for i in list(range(len(rows))):
        rows[i] = [replacement.get(n, n) for n in rows[i]]

    tester3 = pd.DataFrame(rows, dtype=str)
    tester3.to_csv(path + "/keys_succession/equiv_results_" + nums[2] + ".csv", index=False)

    return (print('Dataset ' + str(files[k]) + " took " + str(time.time() - debut) + ' seconds'))

#endregion

#region Functions for successions
def equi_landscape4_recover(essai):
    """
    This function computes the unique matrices equivalent to the current 4x4 landscape using :
    - the initial matrix and 3 consecutive 90° rotations
    - A matrix transformed using horizontal symmetry and 3 consecutive 90° rotations
    - A matrix transformed using vertical symmetry
    """
    # 3 transformations
    m2 = np.rot90(essai)
    m3 = np.rot90(m2)
    m4 = np.rot90(m3)
    # Horizontal symmetry
    hor_sym = [essai[3],essai[2],essai[1],essai[0]]

    # 3 clockwise transformations
    m6 = np.rot90(hor_sym)
    m7 = np.rot90(m6)
    m8 = np.rot90(m7)
    evalu = [tuple(list(itertools.chain(*m2))),
             tuple(list(itertools.chain(*m3))),
             tuple(list(itertools.chain(*m4))),
             tuple(list(itertools.chain(*hor_sym))),
             tuple(list(itertools.chain(*m6))),
             tuple(list(itertools.chain(*m7))),
             tuple(list(itertools.chain(*m8))),
             tuple(list(itertools.chain(*essai)))]
    return evalu
def equi_landscape3_recover(essai):
    """
    This function computes the unique matrices equivalent to the current 4x4 landscape using :
    - the initial matrix and 3 consecutive 90° rotations
    - A matrix transformed using horizontal symmetry and 3 consecutive 90° rotations
    - A matrix transformed using vertical symmetry
    """
    # 3 transformations
    m2 = np.rot90(essai)
    m3 = np.rot90(m2)
    m4 = np.rot90(m3)
    # Horizontal symmetry
    hor_sym = [essai[2],essai[1],essai[0]]

    # 3 clockwise transformations
    m6 = np.rot90(hor_sym)
    m7 = np.rot90(m6)
    m8 = np.rot90(m7)
    evalu = [tuple(list(itertools.chain(*m2))),
             tuple(list(itertools.chain(*m3))),
             tuple(list(itertools.chain(*m4))),
             tuple(list(itertools.chain(*hor_sym))),
             tuple(list(itertools.chain(*m6))),
             tuple(list(itertools.chain(*m7))),
             tuple(list(itertools.chain(*m8))),
             tuple(list(itertools.chain(*essai)))]
    return evalu

def recover(land, transfo):
    land = to_matrix(land,4)
    if transfo ==0:
        m2 = np.rot90(land,axes=(1,0))
        return(tuple(list(itertools.chain(*m2))))
    if transfo ==1:
        m3 = np.rot90(np.rot90(land,axes=(1,0)),axes=(1,0))
        return(tuple(list(itertools.chain(*m3))))
    if transfo ==2:
        m4 = np.rot90(np.rot90(np.rot90(land,axes=(1,0)),axes=(1,0)),axes=(1,0))
        return(tuple(list(itertools.chain(*m4))))
    if transfo ==3:
        m5 = [land[3],land[2],land[1],land[0]]
        return(tuple(list(itertools.chain(*m5))))
    if transfo ==4:
        m6 = np.rot90(land ,axes=(1,0))
        m6 = [m6[3],m6[2],m6[1],m6[0]]
        return(tuple(list(itertools.chain(*m6))))
    if transfo ==5:
        m7 = np.rot90(np.rot90(land,axes=(1,0)),axes=(1,0))
        m7 = [m7[3],m7[2],m7[1],m7[0]]
        return(tuple(list(itertools.chain(*m7))))
    if transfo ==6:
        m8 = np.rot90(np.rot90(np.rot90(land,axes=(1,0)),axes=(1,0)),axes=(1,0))
        m8 = [m8[3],m8[2],m8[1],m8[0]]
        return(tuple(list(itertools.chain(*m8))))
    if transfo ==7:
        return(tuple(list(itertools.chain(*land))))
    if transfo=="None":
        return(tuple(list(itertools.chain(*land))))
def recover3(land, transfo):
    land = to_matrix(land,3)
    if transfo ==0:
        m2 = np.rot90(land,axes=(1,0))
        return(tuple(list(itertools.chain(*m2))))
    if transfo ==1:
        m3 = np.rot90(np.rot90(land,axes=(1,0)),axes=(1,0))
        return(tuple(list(itertools.chain(*m3))))
    if transfo ==2:
        m4 = np.rot90(np.rot90(np.rot90(land,axes=(1,0)),axes=(1,0)),axes=(1,0))
        return(tuple(list(itertools.chain(*m4))))
    if transfo ==3:
        m5 = [land[2],land[1],land[0]]
        return(tuple(list(itertools.chain(*m5))))
    if transfo ==4:
        m6 = np.rot90(land ,axes=(1,0))
        m6 = [m6[2],m6[1],m6[0]]
        return(tuple(list(itertools.chain(*m6))))
    if transfo ==5:
        m7 = np.rot90(np.rot90(land,axes=(1,0)),axes=(1,0))
        m7 = [m7[2],m7[1],m7[0]]
        return(tuple(list(itertools.chain(*m7))))
    if transfo ==6:
        m8 = np.rot90(np.rot90(np.rot90(land,axes=(1,0)),axes=(1,0)),axes=(1,0))
        m8 = [m8[2],m8[1],m8[0]]
        return(tuple(list(itertools.chain(*m8))))
    if transfo ==7:
        return(tuple(list(itertools.chain(*land))))
    if transfo=="None":
        return(tuple(list(itertools.chain(*land))))

def recover_full(land,transfo,index):
    """
    Uses a list of transformations and the index of the land we're trying to recover.
    :param land:
    :param transfo:
    :param index:
    :return:
    """
    land3 = land
    if index ==0:
        return(land3)
    if index ==1:
        return(recover(land3,transfo[1]))
    if index ==2:
        return(recover(recover(land3,transfo[2]),transfo[1]))
    if index ==3:
        return(recover(recover(recover(land3,transfo[3]),transfo[2]),transfo[1]))
    if index ==4:
        return(recover(recover(recover(recover(land3,transfo[4]),transfo[3]),transfo[2]),transfo[1]))
    if index ==5:
        return(recover(recover(recover(recover(recover(land3,transfo[5]),transfo[4]),transfo[3]),transfo[2]),transfo[1]))
    if index ==6:
        return(recover(recover(recover(recover(recover(recover(land3,transfo[6]),transfo[5]),transfo[4]),transfo[3]),transfo[2]),transfo[1]))
    if index ==7 :
        return(recover(recover(recover(recover(recover(recover(recover(land3,transfo[7]),transfo[6]),transfo[5]),transfo[4]),transfo[3]),transfo[2]),transfo[1]))
    if index ==8:
        return(recover(recover(recover(recover(recover(recover(recover(recover(land3,transfo[8]),transfo[7]),transfo[6]),transfo[5]),transfo[4]),transfo[3]),transfo[2]),transfo[1]))
    if index ==9:
        return(recover(recover(recover(recover(recover(recover(recover(recover(recover(land3,transfo[9]),transfo[8]),transfo[7]),transfo[6]),transfo[5]),transfo[4]),transfo[3]),transfo[2]),transfo[1]))
    if index ==10:
        return(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(land3,transfo[10]),transfo[9]),transfo[8]),transfo[7]),transfo[6]),transfo[5]),transfo[4]),transfo[3]),transfo[2]),transfo[1]))
    if index ==11:
        return(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(land3,transfo[11]),transfo[10]),transfo[9]),transfo[8]),transfo[7]),transfo[6]),transfo[5]),transfo[4]),transfo[3]),transfo[2]),transfo[1]))
    if index ==12:
        return(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(land3,transfo[12]),transfo[11]),transfo[10]),transfo[9]),transfo[8]),transfo[7]),transfo[6]),transfo[5]),transfo[4]),transfo[3]),transfo[2]),transfo[1]))
    if index ==13:
        return(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(land3,transfo[13]),transfo[12]),transfo[11]),transfo[10]),transfo[9]),transfo[8]),transfo[7]),transfo[6]),transfo[5]),transfo[4]),transfo[3]),transfo[2]),transfo[1]))
    if index ==14:
        return(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(land3,transfo[14]),transfo[13]),transfo[12]),transfo[11]),transfo[10]),transfo[9]),transfo[8]),transfo[7]),transfo[6]),transfo[5]),transfo[4]),transfo[3]),transfo[2]),transfo[1]))
    if index ==15:
        return(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(land3,transfo[15]),transfo[14]),transfo[13]),transfo[12]),transfo[11]),transfo[10]),transfo[9]),transfo[8]),transfo[7]),transfo[6]),transfo[5]),transfo[4]),transfo[3]),transfo[2]),transfo[1]))
    if index ==16:
        return(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(land3,transfo[16]),transfo[15]),transfo[14]),transfo[13]),transfo[12]),transfo[11]),transfo[10]),transfo[9]),transfo[8]),transfo[7]),transfo[6]),transfo[5]),transfo[4]),transfo[3]),transfo[2]),transfo[1]))
    if index == 17:
        return(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(land3,transfo[17]),transfo[16]),transfo[15]),transfo[14]),transfo[13]),transfo[12]),transfo[11]),transfo[10]),transfo[9]),transfo[8]),transfo[7]),transfo[6]),transfo[5]),transfo[4]),transfo[3]),transfo[2]),transfo[1]))
    if index == 18:
        return(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(land3,transfo[18]),transfo[17]),transfo[16]),transfo[15]),transfo[14]),transfo[13]),transfo[12]),transfo[11]),transfo[10]),transfo[9]),transfo[8]),transfo[7]),transfo[6]),transfo[5]),transfo[4]),transfo[3]),transfo[2]),transfo[1]))
    if index == 19:
        return(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(land3,transfo[19]),transfo[18]),transfo[17]),transfo[16]),transfo[15]),transfo[14]),transfo[13]),transfo[12]),transfo[11]),transfo[10]),transfo[9]),transfo[8]),transfo[7]),transfo[6]),transfo[5]),transfo[4]),transfo[3]),transfo[2]),transfo[1]))
def recover_full2(land,transfo,index):
    """
    Uses a list of transformations and the index of the land we're trying to recover.
    :param land:
    :param transfo:
    :param index:
    :return:
    """
    land3 = land
    if index ==0:
        return(land3)
    if index ==1:
        return(recover(land3,transfo[0]))
    if index ==2:
        return(recover(recover(land3,transfo[1]),transfo[0]))
    if index ==3:
        return(recover(recover(recover(land3,transfo[2]),transfo[1]),transfo[0]))
    if index ==4:
        return(recover(recover(recover(recover(land3,transfo[3]),transfo[2]),transfo[1]),transfo[0]))
    if index ==5:
        return(recover(recover(recover(recover(recover(land3,transfo[4]),transfo[3]),transfo[2]),transfo[1]),transfo[0]))
    if index ==6:
        return(recover(recover(recover(recover(recover(recover(land3,transfo[5]),transfo[4]),transfo[3]),transfo[2]),transfo[1]),transfo[0]))
    if index ==7 :
        return(recover(recover(recover(recover(recover(recover(recover(land3,transfo[6]),transfo[5]),transfo[4]),transfo[3]),transfo[2]),transfo[1]),transfo[0]))
    if index ==8:
        return(recover(recover(recover(recover(recover(recover(recover(recover(land3,transfo[7]),transfo[6]),transfo[5]),transfo[4]),transfo[3]),transfo[2]),transfo[1]),transfo[0]))
    if index ==9:
        return(recover(recover(recover(recover(recover(recover(recover(recover(recover(land3,transfo[8]),transfo[7]),transfo[6]),transfo[5]),transfo[4]),transfo[3]),transfo[2]),transfo[1]),transfo[0]))
    if index ==10:
        return(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(land3,transfo[9]),transfo[8]),transfo[7]),transfo[6]),transfo[5]),transfo[4]),transfo[3]),transfo[2]),transfo[1]),transfo[0]))
    if index ==11:
        return(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(land3,transfo[10]),transfo[9]),transfo[8]),transfo[7]),transfo[6]),transfo[5]),transfo[4]),transfo[3]),transfo[2]),transfo[1]),transfo[0]))
    if index ==12:
        return(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(land3,transfo[11]),transfo[10]),transfo[9]),transfo[8]),transfo[7]),transfo[6]),transfo[5]),transfo[4]),transfo[3]),transfo[2]),transfo[1]),transfo[0]))
    if index ==13:
        return(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(land3,transfo[12]),transfo[11]),transfo[10]),transfo[9]),transfo[8]),transfo[7]),transfo[6]),transfo[5]),transfo[4]),transfo[3]),transfo[2]),transfo[1]),transfo[0]))
    if index ==14:
        return(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(land3,transfo[13]),transfo[12]),transfo[11]),transfo[10]),transfo[9]),transfo[8]),transfo[7]),transfo[6]),transfo[5]),transfo[4]),transfo[3]),transfo[2]),transfo[1]),transfo[0]))
    if index ==15:
        return(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(land3,transfo[14]),transfo[13]),transfo[12]),transfo[11]),transfo[10]),transfo[9]),transfo[8]),transfo[7]),transfo[6]),transfo[5]),transfo[4]),transfo[3]),transfo[2]),transfo[1]),transfo[0]))
    if index ==16:
        return(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(land3,transfo[15]),transfo[14]),transfo[13]),transfo[12]),transfo[11]),transfo[10]),transfo[9]),transfo[8]),transfo[7]),transfo[6]),transfo[5]),transfo[4]),transfo[3]),transfo[2]),transfo[1]),transfo[0]))
    if index == 17:
        return(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(land3,transfo[16]),transfo[15]),transfo[14]),transfo[13]),transfo[12]),transfo[11]),transfo[10]),transfo[9]),transfo[8]),transfo[7]),transfo[6]),transfo[5]),transfo[4]),transfo[3]),transfo[2]),transfo[1]),transfo[0]))
    if index == 18:
        return(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(land3,transfo[17]),transfo[16]),transfo[15]),transfo[14]),transfo[13]),transfo[12]),transfo[11]),transfo[10]),transfo[9]),transfo[8]),transfo[7]),transfo[6]),transfo[5]),transfo[4]),transfo[3]),transfo[2]),transfo[1]),transfo[0]))
    if index == 19:
        return(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(recover(land3,transfo[18]),transfo[17]),transfo[16]),transfo[15]),transfo[14]),transfo[13]),transfo[12]),transfo[11]),transfo[10]),transfo[9]),transfo[8]),transfo[7]),transfo[6]),transfo[5]),transfo[4]),transfo[3]),transfo[2]),transfo[1]),transfo[0]))
def recover_full3(land,transfo,index):
    """
    Uses a list of transformations and the index of the land we're trying to recover.
    :param land:
    :param transfo:
    :param index:
    :return:
    """
    land3 = land
    if index ==0:
        return(land3)
    if index ==1:
        return(recover3(land3,transfo[0]))
    if index ==2:
        return(recover3(recover3(land3,transfo[1]),transfo[0]))
    if index ==3:
        return(recover3(recover3(recover3(land3,transfo[2]),transfo[1]),transfo[0]))
    if index ==4:
        return(recover3(recover3(recover3(recover3(land3,transfo[3]),transfo[2]),transfo[1]),transfo[0]))
    if index ==5:
        return(recover3(recover3(recover3(recover3(recover3(land3,transfo[4]),transfo[3]),transfo[2]),transfo[1]),transfo[0]))
    if index ==6:
        return(recover3(recover3(recover3(recover3(recover3(recover3(land3,transfo[5]),transfo[4]),transfo[3]),transfo[2]),transfo[1]),transfo[0]))
    if index ==7 :
        return(recover3(recover3(recover3(recover3(recover3(recover3(recover3(land3,transfo[6]),transfo[5]),transfo[4]),transfo[3]),transfo[2]),transfo[1]),transfo[0]))
    if index ==8:
        return(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(land3,transfo[7]),transfo[6]),transfo[5]),transfo[4]),transfo[3]),transfo[2]),transfo[1]),transfo[0]))
    if index ==9:
        return(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(land3,transfo[8]),transfo[7]),transfo[6]),transfo[5]),transfo[4]),transfo[3]),transfo[2]),transfo[1]),transfo[0]))
    if index ==10:
        return(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(land3,transfo[9]),transfo[8]),transfo[7]),transfo[6]),transfo[5]),transfo[4]),transfo[3]),transfo[2]),transfo[1]),transfo[0]))
    if index ==11:
        return(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(land3,transfo[10]),transfo[9]),transfo[8]),transfo[7]),transfo[6]),transfo[5]),transfo[4]),transfo[3]),transfo[2]),transfo[1]),transfo[0]))
    if index ==12:
        return(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(land3,transfo[11]),transfo[10]),transfo[9]),transfo[8]),transfo[7]),transfo[6]),transfo[5]),transfo[4]),transfo[3]),transfo[2]),transfo[1]),transfo[0]))
    if index ==13:
        return(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(land3,transfo[12]),transfo[11]),transfo[10]),transfo[9]),transfo[8]),transfo[7]),transfo[6]),transfo[5]),transfo[4]),transfo[3]),transfo[2]),transfo[1]),transfo[0]))
    if index ==14:
        return(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(land3,transfo[13]),transfo[12]),transfo[11]),transfo[10]),transfo[9]),transfo[8]),transfo[7]),transfo[6]),transfo[5]),transfo[4]),transfo[3]),transfo[2]),transfo[1]),transfo[0]))
    if index ==15:
        return(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(land3,transfo[14]),transfo[13]),transfo[12]),transfo[11]),transfo[10]),transfo[9]),transfo[8]),transfo[7]),transfo[6]),transfo[5]),transfo[4]),transfo[3]),transfo[2]),transfo[1]),transfo[0]))
    if index ==16:
        return(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(land3,transfo[15]),transfo[14]),transfo[13]),transfo[12]),transfo[11]),transfo[10]),transfo[9]),transfo[8]),transfo[7]),transfo[6]),transfo[5]),transfo[4]),transfo[3]),transfo[2]),transfo[1]),transfo[0]))
    if index == 17:
        return(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(land3,transfo[16]),transfo[15]),transfo[14]),transfo[13]),transfo[12]),transfo[11]),transfo[10]),transfo[9]),transfo[8]),transfo[7]),transfo[6]),transfo[5]),transfo[4]),transfo[3]),transfo[2]),transfo[1]),transfo[0]))
    if index == 18:
        return(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(land3,transfo[17]),transfo[16]),transfo[15]),transfo[14]),transfo[13]),transfo[12]),transfo[11]),transfo[10]),transfo[9]),transfo[8]),transfo[7]),transfo[6]),transfo[5]),transfo[4]),transfo[3]),transfo[2]),transfo[1]),transfo[0]))
    if index == 19:
        return(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(recover3(land3,transfo[18]),transfo[17]),transfo[16]),transfo[15]),transfo[14]),transfo[13]),transfo[12]),transfo[11]),transfo[10]),transfo[9]),transfo[8]),transfo[7]),transfo[6]),transfo[5]),transfo[4]),transfo[3]),transfo[2]),transfo[1]),transfo[0]))
def succession_dataframe(data_path, biodiv_index,time_param='Yes',follow_up="No"):
    budget = globals()["budget"]
    debut = time.time()
    #region Load the current data

    file = open("C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_" + str(budget) + "/"+data_path)
    csvreader = csv.reader(file)
    header = next(csvreader)

    rows_data = []
    for row in csvreader:
        row_add = [row[0],row[int(biodiv_index/2)], row[int(82+biodiv_index/2)]]
        rows_data.append(row_add)
    #endregion
    #region Import data to scan through - dict format and list format
    path = "C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_" + str(budget) + "/"
    files = os.listdir(path)
    files.remove('keys_succession')
    files.remove('successions_matched')
    try:
        files.remove("updated_list.csv")
    except ValueError:
        pass

    files = [x for x in files if
             ast.literal_eval(re.findall(r'\d+', str(x))[2]) in list(range(16 - budget, 17))]
    datas = {}
    datas[data_path]=rows_data
    for x in files:
        file = open("C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_" + str(budget) + "/" +x)
        csvreader = csv.reader(file)
        header = next(csvreader)

        rows_dict = []
        for row in csvreader:
            row_add = [row[0],row[int(biodiv_index/2)], row[int(82+biodiv_index/2)]]
            rows_dict.append(row_add)
        datas[x] = rows_dict
    # endregion

    # region Import equivalent data - dict format and list formats
    path = "C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_" + str(budget) + "/keys_succession/"
    files = os.listdir(path)
    files = [x for x in files if
             ast.literal_eval(re.findall(r'\d+', str(x))[0]) in list(range(16 - budget, 17))]
    datas_equiv = {}
    for x in files:
        # region Load data in list format
        file = open(path + x)
        csvreader = csv.reader(file)
        header = next(csvreader)

        rows_dict = []
        for row in csvreader:
            row_add = row[int(biodiv_index/2-1)]
            rows_dict.append(row_add)
        datas_equiv[x] = rows_dict
        #endregion
    # endregion
    #region Add current data to dict repositories

    match_equiv = re.findall(r'\d+', str(data_path))
    if len(match_equiv)>4:
        file = open("C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_" + str(budget) + "/keys_succession/equiv_results_"+str(match_equiv[2])+"_"+str(match_equiv[3])+"_cut_"+str(match_equiv[4])+".csv")
        csvreader = csv.reader(file)
        header = next(csvreader)

        rows_dict = []
        for row in csvreader:
            row_add = row[int(biodiv_index/2-1)]
            rows_dict.append(row_add)
        datas_equiv[('equiv_results_'+str(match_equiv[2])+"_"+str(match_equiv[3])+"_cut_"+str(match_equiv[4])+'.csv')]=rows_dict
    else:
        file = open("C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_" + str(budget) + "/keys_succession/equiv_results_"+str(match_equiv[2])+"_"+str(match_equiv[3])+".csv")
        csvreader = csv.reader(file)
        header = next(csvreader)

        rows_dict = []
        for row in csvreader:
            row_add = row[int(biodiv_index/2-1)]
            rows_dict.append(row_add)
        datas_equiv[('equiv_results_'+str(match_equiv[2])+"_"+str(match_equiv[3])+'.csv')]=rows_dict

    #endregion

    rows_out = []

    for p in list(range(len(rows_data))):
    # input
        land = str(rows_data[p][0])
    # output
        succession_list = [land]
        transformation = ['None']
        equiv_key = [land]
        value = [np.nan]
    #Meta parameters
        path = "C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_" + str(budget) + "/"
        step = 0
        while step <= 20:
            check = nonzero(np.array(ast.literal_eval(equiv_key[-1])))
            check2 = nb2(np.array(ast.literal_eval(equiv_key[-1])))
            #print(check, check2)
    # use relevant result data from loaded data :
            namer = ("land4_budget_" + str(budget) + "_" + str(check) + "_" + str(check2) + ".csv")
            namer2 = ("land4_budget_" + str(budget) + "_" + str(check) + "_" + str(check2) + "_cut_0.csv")

            if os.path.exists(path + namer):
                if namer in list(datas.keys()):
                    d2 = datas[namer]
                else:
                    file = open(path + "land4_budget_" + str(budget) + "_" + str(check) + "_" + str(check2) + ".csv")
                    csvreader = csv.reader(file)
                    header = next(csvreader)

                    rows_dict = []
                    for row in csvreader:
                        row_add = [row[0],row[int(biodiv_index/2)], row[int(82+biodiv_index/2)]]
                        rows_dict.append(row_add)
                    datas["land4_budget_" + str(budget) + "_" + str(check) + "_" + str(check2) + ".csv"] = rows_dict
                    d2 = rows_dict

                for g in list(range(len(d2))):
                    if d2[g][0]==equiv_key[-1]:
                        succession_index=g

                name2 = "equiv_results_" + str(check) + "_" + str(check2) + ".csv"

            elif os.path.exists(path + namer2):
                if namer2 in list(datas.keys()):
                    to_load = [x for x in list(datas.keys()) if x.startswith('land4_budget_' + str(budget) + "_" + str(check) + "_" + str(check2))]
                else :
                    path = "C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_" + str(
                        budget) + "/"
                    files = os.listdir(path)
                    files.remove('keys_succession')
                    files.remove('successions_matched')
                    to_load = [x for x in files if
                               x.startswith('land4_budget_' + str(budget) + "_" + str(check) + "_" + str(check2))]
                    for y in to_load:
                        file = open(path + y)
                        csvreader = csv.reader(file)
                        header = next(csvreader)

                        rows_dict = []
                        for row in csvreader:
                            row_add = [row[0], row[int(biodiv_index / 2)], row[int(82 + biodiv_index / 2)]]
                            rows_dict.append(row_add)
                        datas[y] = rows_dict

                for x in to_load:
                    d2 = datas[x]
                    for g in list(range(len(d2))):
                        if d2[g][0]==equiv_key[-1]:
                            succession_index = g
                            indexes = re.findall(r'\d+', str(x))
                            keep = x
                            name2 = "equiv_results_" + str(check) + "_" + str(check2) + "_cut_" + str(
                                indexes[4]) + ".csv"

                d2=datas[keep]
            succession_key = d2[succession_index][1]
            value_column = int(82 + biodiv_index / 2)

            if name2 in list(datas_equiv.keys()):
                d_equiv=datas_equiv[name2]
            else:
                file = open(path + "keys_succession/"+name2)
                csvreader = csv.reader(file)
                header = next(csvreader)

                rows_dict = []
                for row in csvreader:
                    row_add = row[int(biodiv_index/2-1)]
                    rows_dict.append(row_add)
                d_equiv = rows_dict
                datas_equiv[name2] = d_equiv

            #region Computations
            candidate_succ = d_equiv[succession_index]
            mat_succession = to_matrix(ast.literal_eval(succession_key),4)
            transfo = equi_landscape_recover(mat_succession,4).index(ast.literal_eval(candidate_succ))
            #endregion
            #region  list updates
            value.append(d2[succession_index][2])
            succession_list.append(succession_key)
            equiv_key.append(candidate_succ)
            transformation.append(transfo)
            #endregion

            step += 1
        #region Final processing : recovery of optimal succession with chain transformation
        c = list(map(ast.literal_eval, succession_list))
        final_succession = []
        for i in list(range(20)):
            e = transformation[0: i ]
            # deleted +1
            final_succession.append(recover_full2(c[i], transfo=e, index=i))

        final_succession = final_succession + value
        rows_out.append(final_succession)
        if follow_up == 'Yes':
            print(str(100*len(rows_out) / len(rows_data))+" %")
        #endregion

    data_path_saved = data_path.replace('.csv', '_biodiv_index_' + str(biodiv_index) + '.csv')
    data_return = pd.DataFrame(rows_out)
    data_return.to_csv(path+'successions_matched/'+data_path_saved)


    if time_param == 'Yes':
        return(data_path + " is done and took "+str(time.time() - debut)+' seconds' )
    else :
        return(data_path + " is done")


def succession_dataframeCHECK(data_path, biodiv_index,time_param='Yes',budget=1,follow_up="No"):

    debut = time.time()
    #region Load the current data

    file = open("C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_" + str(budget) + "/"+data_path)
    csvreader = csv.reader(file)
    header = next(csvreader)

    rows_data = []
    for row in csvreader:
        row_add = [row[0],row[int(biodiv_index/2)], row[int(82+biodiv_index/2)]]
        rows_data.append(row_add)
    #endregion
    #region Import data to scan through - dict format and list format
    path = "C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_" + str(budget) + "/"
    files = os.listdir(path)
    files.remove('keys_succession')
    files.remove('successions_matched')
    files = [x for x in files if
             ast.literal_eval(re.findall(r'\d+', str(x))[2]) in list(range(16 - budget, 17))]
    datas = {}
    datas[data_path]=rows_data
    for x in files:
        file = open("C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_" + str(budget) + "/" +x)
        csvreader = csv.reader(file)
        header = next(csvreader)

        rows_dict = []
        for row in csvreader:
            row_add = [row[0],row[int(biodiv_index/2)], row[int(82+biodiv_index/2)]]
            rows_dict.append(row_add)
        datas[x] = rows_dict
    # endregion

    # region Import equivalent data - dict format and list formats
    path = "C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_" + str(budget) + "/keys_succession/"
    files = os.listdir(path)
    files = [x for x in files if
             ast.literal_eval(re.findall(r'\d+', str(x))[0]) in list(range(16 - budget, 17))]
    datas_equiv = {}
    for x in files:
        # region Load data in list format
        file = open(path + x)
        csvreader = csv.reader(file)
        header = next(csvreader)

        rows_dict = []
        for row in csvreader:
            row_add = row[int(biodiv_index/2-1)]
            rows_dict.append(row_add)
        datas_equiv[x] = rows_dict
        #endregion
    # endregion
    #region Add current data to dict repositories

    match_equiv = re.findall(r'\d+', str(data_path))
    if len(match_equiv)>4:
        file = open("C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_" + str(budget) + "/keys_succession/equiv_results_"+str(match_equiv[2])+"_"+str(match_equiv[3])+"_cut_"+str(match_equiv[4])+".csv")
        csvreader = csv.reader(file)
        header = next(csvreader)

        rows_dict = []
        for row in csvreader:
            row_add = row[int(biodiv_index/2-1)]
            rows_dict.append(row_add)
        datas_equiv[('equiv_results_'+str(match_equiv[2])+"_"+str(match_equiv[3])+"_cut_"+str(match_equiv[4])+'.csv')]=rows_dict
    else:
        file = open("C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_" + str(budget) + "/keys_succession/equiv_results_"+str(match_equiv[2])+"_"+str(match_equiv[3])+".csv")
        csvreader = csv.reader(file)
        header = next(csvreader)

        rows_dict = []
        for row in csvreader:
            row_add = row[int(biodiv_index/2-1)]
            rows_dict.append(row_add)
        datas_equiv[('equiv_results_'+str(match_equiv[2])+"_"+str(match_equiv[3])+'.csv')]=rows_dict

    #endregion

    rows_out = []

    for p in list(range(len(rows_data))):
    # input
        land = str(rows_data[p][0])
    # output
        succession_list = [land]
        transformation = ['None']
        equiv_key = [land]
        value = [np.nan]
    #Meta parameters
        path = "C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_" + str(budget) + "/"
        step = 0
        while step <= 20:
            check = nonzero(np.array(ast.literal_eval(equiv_key[-1])))
            check2 = nb2(np.array(ast.literal_eval(equiv_key[-1])))
            #print(check, check2)
    # use relevant result data from loaded data :
            namer = ("land4_budget_" + str(budget) + "_" + str(check) + "_" + str(check2) + ".csv")
            namer2 = ("land4_budget_" + str(budget) + "_" + str(check) + "_" + str(check2) + "_cut_0.csv")

            if os.path.exists(path + namer):
                if namer in list(datas.keys()):
                    d2 = datas[namer]
                else:
                    file = open(path + "land4_budget_" + str(budget) + "_" + str(check) + "_" + str(check2) + ".csv")
                    csvreader = csv.reader(file)
                    header = next(csvreader)

                    rows_dict = []
                    for row in csvreader:
                        row_add = [row[0],row[int(biodiv_index/2)], row[int(82+biodiv_index/2)]]
                        rows_dict.append(row_add)
                    datas["land4_budget_" + str(budget) + "_" + str(check) + "_" + str(check2) + ".csv"] = rows_dict
                    d2 = rows_dict

                for g in list(range(len(d2))):
                    if d2[g][0]==equiv_key[-1]:
                        succession_index=g

                name2 = "equiv_results_" + str(check) + "_" + str(check2) + ".csv"

            elif os.path.exists(path + namer2):
                if namer2 in list(datas.keys()):
                    to_load = [x for x in list(datas.keys()) if x.startswith('land4_budget_' + str(budget) + "_" + str(check) + "_" + str(check2))]
                else :
                    path = "C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_" + str(
                        budget) + "/"
                    files = os.listdir(path)
                    files.remove('keys_succession')
                    files.remove('successions_matched')
                    to_load = [x for x in files if
                               x.startswith('land4_budget_' + str(budget) + "_" + str(check) + "_" + str(check2))]
                    for y in to_load:
                        file = open(path + y)
                        csvreader = csv.reader(file)
                        header = next(csvreader)

                        rows_dict = []
                        for row in csvreader:
                            row_add = [row[0], row[int(biodiv_index / 2)], row[int(82 + biodiv_index / 2)]]
                            rows_dict.append(row_add)
                        datas[y] = rows_dict

                for x in to_load:
                    d2 = datas[x]
                    for g in list(range(len(d2))):
                        if d2[g][0]==equiv_key[-1]:
                            succession_index = g
                            indexes = re.findall(r'\d+', str(x))
                            keep = x
                            name2 = "equiv_results_" + str(check) + "_" + str(check2) + "_cut_" + str(
                                indexes[4]) + ".csv"

                d2=datas[keep]
            succession_key = d2[succession_index][1]
            value_column = int(82 + biodiv_index / 2)

            if name2 in list(datas_equiv.keys()):
                d_equiv=datas_equiv[name2]
            else:
                file = open(path + "keys_succession/"+name2)
                csvreader = csv.reader(file)
                header = next(csvreader)

                rows_dict = []
                for row in csvreader:
                    row_add = row[int(biodiv_index/2-1)]
                    rows_dict.append(row_add)
                d_equiv = rows_dict
                datas_equiv[name2] = d_equiv

            #region Computations
            candidate_succ = d_equiv[succession_index]
            mat_succession = to_matrix(ast.literal_eval(succession_key),4)
            transfo = equi_landscape_recover(mat_succession,4).index(ast.literal_eval(candidate_succ))
            #endregion
            #region  list updates
            value.append(d2[succession_index][2])
            succession_list.append(succession_key)
            equiv_key.append(candidate_succ)
            transformation.append(transfo)
            #endregion

            step += 1
        #region Final processing : recovery of optimal succession with chain transformation
        c = list(map(ast.literal_eval, succession_list))
        final_succession = []
        for i in list(range(20)):
            e = transformation[0: i ]
            # deleted +1
            final_succession.append(recover_full2(c[i], transfo=e, index=i))

        final_succession = final_succession + value
        rows_out.append(final_succession)
        if follow_up == 'Yes':
            print(str(100*len(rows_out) / len(rows_data))+" %")
        #endregion

    data_path_saved = data_path.replace('.csv', '_biodiv_index_' + str(biodiv_index) + '.csv')
    data_return = pd.DataFrame(rows_out)
    data_return.to_csv(path+'successions_matched/NEW'+data_path_saved)


    if time_param == 'Yes':
        return(data_path + " is done and took "+str(time.time() - debut)+' seconds' )
    else :
        return(data_path + " is done")
def succession_dataframe3(data_path, biodiv_index,time_param='Yes',follow_up="No"):
    budget = globals()['budget']
    debut = time.time()
    #region Load the current data

    file = open("C:/Users/jean/PycharmProjects/connectivity_dilemma/data3/results/budget_" + str(budget) + "/"+data_path)
    csvreader = csv.reader(file)
    header = next(csvreader)

    rows_data = []
    for row in csvreader:
        row_add = [row[0],row[int(biodiv_index/2)], row[int(51+biodiv_index/2)]]
        rows_data.append(row_add)
    #endregion
    #region Import data to scan through - dict format and list format
    path = "C:/Users/jean/PycharmProjects/connectivity_dilemma/data3/results/budget_" + str(budget) + "/"
    files = os.listdir(path)
    files.remove('keys_succession')
    files.remove('successions_matched')
    files = [x for x in files if
             ast.literal_eval(re.findall(r'\d+', str(x))[2]) in list(range(9 - budget, 10))]
    datas = {}
    datas[data_path]=rows_data
    for x in files:
        file = open("C:/Users/jean/PycharmProjects/connectivity_dilemma/data3/results/budget_" + str(budget) + "/" +x)
        csvreader = csv.reader(file)
        header = next(csvreader)

        rows_dict = []
        for row in csvreader:
            row_add = [row[0],row[int(biodiv_index/2)], row[int(51+biodiv_index/2)]]
            rows_dict.append(row_add)
        datas[x] = rows_dict
    # endregion

    # region Import equivalent data - dict format and list formats
    path = "C:/Users/jean/PycharmProjects/connectivity_dilemma/data3/results/budget_" + str(budget) + "/keys_succession/"
    files = os.listdir(path)
    files = [x for x in files if
             ast.literal_eval(re.findall(r'\d+', str(x))[0]) in list(range(9 - budget, 10))]
    datas_equiv = {}
    for x in files:
        # region Load data in list format
        file = open(path + x)
        csvreader = csv.reader(file)
        header = next(csvreader)

        rows_dict = []
        for row in csvreader:
            row_add = row[int(biodiv_index/2)]
            rows_dict.append(row_add)
        datas_equiv[x] = rows_dict
        #endregion
    # endregion
    #region Add current data to dict repositories

    match_equiv = re.findall(r'\d+', str(data_path))
    file = open("C:/Users/jean/PycharmProjects/connectivity_dilemma/data3/results/budget_" + str(budget) + "/keys_succession/equiv_results_"+str(match_equiv[2])+".csv")
    csvreader = csv.reader(file)
    header = next(csvreader)

    rows_dict = []
    for row in csvreader:
        row_add = row[int(biodiv_index/2)]
        rows_dict.append(row_add)
    datas_equiv[('equiv_results_'+str(match_equiv[2])+'.csv')]=rows_dict

    #endregion

    rows_out = []

    for p in list(range(len(rows_data))):
    # input
        land = str(rows_data[p][0])
    # output
        succession_list = [land]
        transformation = ['None']
        equiv_key = [land]
        value = [np.nan]
    #Meta parameters
        path = "C:/Users/jean/PycharmProjects/connectivity_dilemma/data3/results/budget_" + str(budget) + "/"
        step = 0
        while step <= 20:
            check = nonzero(np.array(ast.literal_eval(equiv_key[-1])))
            #print(check, check2)
    # use relevant result data from loaded data :
            namer = ("land3_budget_" + str(budget) + "_" + str(check) + ".csv")

            if os.path.exists(path + namer):
                if namer in list(datas.keys()):
                    d2 = datas[namer]
                else:
                    file = open(path + "land3_budget_" + str(budget) + "_" + str(check)+ ".csv")
                    csvreader = csv.reader(file)
                    header = next(csvreader)

                    rows_dict = []
                    for row in csvreader:
                        row_add = [row[0],row[int(biodiv_index/2)], row[int(51+biodiv_index/2)]]
                        rows_dict.append(row_add)
                    datas["land3_budget_" + str(budget) + "_" + str(check) + ".csv"] = rows_dict
                    d2 = rows_dict

                for g in list(range(len(d2))):
                    if d2[g][0]==equiv_key[-1]:
                        succession_index=g

                name2 = "equiv_results_" + str(check) + ".csv"

            succession_key = d2[succession_index][1]
            value_column = int(51 + biodiv_index / 2)

            if name2 in list(datas_equiv.keys()):
                d_equiv=datas_equiv[name2]
            else:
                file = open(path + "keys_succession/"+name2)
                csvreader = csv.reader(file)
                header = next(csvreader)

                rows_dict = []
                for row in csvreader:
                    row_add = row[int(biodiv_index/2)]
                    rows_dict.append(row_add)
                d_equiv = rows_dict
                datas_equiv[name2] = d_equiv

            #region Computations
            candidate_succ = d_equiv[succession_index]
            mat_succession = to_matrix(ast.literal_eval(succession_key),3)
            transfo = equi_landscape_recover(mat_succession,3).index(ast.literal_eval(candidate_succ))
            #endregion
            #region  list updates
            value.append(d2[succession_index][2])
            succession_list.append(succession_key)
            equiv_key.append(candidate_succ)
            transformation.append(transfo)
            #endregion

            step += 1
        #region Final processing : recovery of optimal succession with chain transformation
        c = list(map(ast.literal_eval, succession_list))
        final_succession = []
        for i in list(range(20)):
            e = transformation[0: i ]
            # deleted +1
            final_succession.append(recover_full3(c[i], transfo=e, index=i))

        final_succession = final_succession + value
        rows_out.append(final_succession)
        if follow_up == 'Yes':
            print(str(100*len(rows_out) / len(rows_data))+" %")
        #endregion

    data_path_saved = data_path.replace('.csv', '_biodiv_index_' + str(biodiv_index) + '.csv')
    data_return = pd.DataFrame(rows_out)
    data_return.to_csv(path+'successions_matched/'+data_path_saved)


    if time_param == 'Yes':
        return(data_path + " is done and took "+str(time.time() - debut)+' seconds' )
    else :
        return(data_path + " is done")
#endregion


#region Potential fires


#endregion

def data_list2_np(input):
    path = "C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/all_data/"
    files = os.listdir(path)

    filename = open(path+files[input], 'rb')
    list_nz = pickle.load(filename)
    filename.close()
    list_nz = list(list_nz)
    list_nz = [np.asarray(x) for x in list_nz]
    return (list_nz)

def lowlev_dynprog_cut_mod_np(land, input):
    # Load the environment
    data_list2_np(input)
    # Set up the output variables
    # compute the new land and associated value
    new_land = fuel_dyn(land, pot_fire_budget)
    value = high_fuel_con_vect(new_land) + high_fuel_connect(land) + 1000000 * (high_biod_con_vect(new_land) <= biodiv)
    # find min
    a = value.min(0)
    b = value.argmin(0)

    store = [ x for x in new_land[np.asarray(b, int)][0].tolist()]
    listed_values= [[land],store,np.asarray(b)[0].tolist(), np.asarray(a)[0].tolist() ]
    list_values = [item for sublist in listed_values for item in sublist]

    return list_values
#################################################################################
class Graph:

    def __init__(self, row, col, g):
        self.ROW = row
        self.COL = col
        self.graph = g

    # A function to check if a given cell
    # (row, col) can be included in DFS
    def isSafe(self, i, j, visited):
        # row number is in range, column number
        # is in range and value is 1
        # and not yet visited
        return (i >= 0 and i < self.ROW and
                j >= 0 and j < self.COL and
                not visited[i][j] and self.graph[i][j])

    # A utility function to do DFS for a 2D
    # boolean matrix. It only considers
    # the 8 neighbours as adjacent vertices

    def DFS(self, i, j, visited):

        # These arrays are used to get row and
        # column numbers of 8 neighbours
        # of a given cell
        rowNbr = [-1, -1, -1, 0, 0, 1, 1, 1]
        colNbr = [-1, 0, 1, -1, 1, -1, 0, 1]

        # Mark this cell as visited
        visited[i][j] = True

        # Recur for all connected neighbours
        for k in range(8):
            if self.isSafe(i + rowNbr[k], j + colNbr[k], visited):
                self.DFS(i + rowNbr[k], j + colNbr[k], visited)

    # The main function that returns
    # count of islands in a given boolean
    # 2D matrix

    def countIslands(self):
        # Make a bool array to mark visited cells.
        # Initially all cells are unvisited
        visited = [[False for j in range(self.COL)]for i in range(self.ROW)]

        # Initialize count as 0 and traverse
        # through the all cells of
        # given matrix
        count = 0
        area = [0]
        zone_l=[]
        for i in range(self.ROW):
            for j in range(self.COL):
                # If a cell with value 1 is not visited yet,
                # then new island found
                if visited[i][j] == False and self.graph[i][j] == 1:
                    # Visit all cells in this island
                    # and increment island count

                    self.DFS(i, j, visited)
                    count += 1
                    area.append(sum([sum(x) for x in visited])-area[-1])
                    zone_l.append((i,j))
        area.pop(0)

        return count, area,zone_l

class QItem:
    def __init__(self, row, col, dist):
        self.row = row
        self.col = col
        self.dist = dist

    def __repr__(self):
        return f"QItem({self.row}, {self.col}, {self.dist})"


def minDistance(grid,sourcex,sourcey,destx,desty):
    source = QItem(sourcex, sourcey, 0)


    # To maintain location visit status
    visited = [[False for _ in range(len(grid[0]))]
               for _ in range(len(grid))]

    # applying BFS on matrix cells starting from source
    queue = []
    queue.append(source)
    visited[source.row][source.col] = True
    while len(queue) != 0:
        source = queue.pop(0)

        # Destination found;
        if(source.row==destx and source.col==desty):
            return(source.dist)
            #print("finito")

        # moving up left
        if isValid(source.row - 1, source.col-1, grid, visited):
            queue.append(QItem(source.row - 1, source.col-1, source.dist + math.sqrt(2)))
            visited[source.row - 1][source.col-1] = True

        # moving up right
        if isValid(source.row - 1, source.col+1, grid, visited):
            queue.append(QItem(source.row - 1, source.col+1, source.dist + math.sqrt(2)))
            visited[source.row - 1][source.col+1] = True

        # moving down right
        if isValid(source.row + 1, source.col + 1, grid, visited):
            queue.append(QItem(source.row + 1, source.col + 1, source.dist + math.sqrt(2)))
            visited[source.row + 1][source.col + 1] = True
        # moving down left
        if isValid(source.row + 1, source.col -1, grid, visited):
            queue.append(QItem(source.row - 1, source.col + 1, source.dist + math.sqrt(2)))
            visited[source.row + 1][source.col - 1] = True

        # moving up
        if isValid(source.row - 1, source.col, grid, visited):
            queue.append(QItem(source.row - 1, source.col, source.dist + 1))
            visited[source.row - 1][source.col] = True

        # moving down
        if isValid(source.row + 1, source.col, grid, visited):
            queue.append(QItem(source.row + 1, source.col, source.dist + 1))
            visited[source.row + 1][source.col] = True

        # moving left
        if isValid(source.row, source.col - 1, grid, visited):
            queue.append(QItem(source.row, source.col - 1, source.dist + 1))
            visited[source.row][source.col - 1] = True

        # moving right
        if isValid(source.row, source.col + 1, grid, visited):
            queue.append(QItem(source.row, source.col + 1, source.dist + 1))
            visited[source.row][source.col + 1] = True

    return 1000
    #dealt with 0 if no path exists, but it is problematic ; if all patches have the same size
    # and no

# checking where move is valid or not
def isValid(x, y, grid, visited):
    if ((x >= 0 and y >= 0) and
            (x < len(grid) and y < len(grid[0])) and
            (grid[x][y] != 0) and (visited[x][y] == False)):
        return True
    return False




class Land:
    def __init__(self,shape):
        self.shape= shape
        self.size = int(math.sqrt(len(shape)))
        self.mat = to_matrix(shape,size)
        self.biod = (shape>=1).dot(connectivity_mat(self.size)).dot(np.transpose((shape>=1)))
        self.biod_relative = self.biod/(np.repeat(1,len(self.shape)).dot(connectivity_mat(self.size)).dot(np.transpose(np.repeat(1,len(self.shape)))))
        self.fuel = (shape>=2).dot(connectivity_mat(self.size)).dot(np.transpose((shape>=2)))
        self.fuel_relative = self.fuel/(np.repeat(1,len(self.shape)).dot(connectivity_mat(self.size)).dot(np.transpose(np.repeat(1,len(self.shape)))))
        self.node_biod = sum(shape>0)
        self.node_fuel = sum(shape==2)
        self.zero = sum(shape==0)
        self.habitat_area_ratio = sum(self.shape>=1)/(self.size**2)
        self.fuel_area_ratio = sum(self.shape==2)/(self.size**2)
    def view(self):
        ax = plt.axes()
        ax.set_title('Initial biodiv = '+str(self.biod))
        sns.heatmap(self.mat,ax = ax, linewidth=0.5, cmap="Greens", vmin=0, vmax=2, center=1)
        plt.show()
    def view_succession(self,index):
        if len(index)==1:
            ax = plt.axes()
            ax.set_title('Biodiv constraint = ' + str(index[0])+" with budget = "+str(budget))
            y = to_matrix(self.succession_land[index[0]], self.size)
            sns.heatmap(y, ax=ax, linewidth=0.5, cmap="Greens", vmin=0, vmax=2, center=1)
            plt.show()
        else:
            number_x = 0
            rows = len(index)//2
            for i in range(rows):
                for j in range(2):
                    ax = plt.subplot2grid((rows, 2), (i, j))
                    ax.title.set_text("Biodiv constraint =" + str(index[number_x]) + ' with budget ='+str(budget))
                    y = to_matrix(self.succession_land[index[number_x]], self.size)
                    ax = sns.heatmap(y, linewidth=0.5, cmap="Greens", vmin=0, vmax=2, center=1)
                    number_x += 1
            plt.show()
    def equivalent(self):
        m2 = np.rot90(self.mat)
        m3 = np.rot90(m2)
        m4 = np.rot90(m3)
        if self.size==5:
             # Horizontal symmetry
            hor_sym = [self.mat[4], self.mat[3], self.mat[2], self.mat[1], self.mat[0]]
        elif self.size==4:
            hor_sym=[self.mat[3], self.mat[2], self.mat[1], self.mat[0]]
        elif self.size==3:
            hor_sym = [self.mat[2], self.mat[1], self.mat[0]]
            # 3 clockwise transformations
        else:
            ValueError("Size is not supported")
        m6 = np.rot90(hor_sym)
        m7 = np.rot90(m6)
        m8 = np.rot90(m7)
        evalu = (tuple(list(itertools.chain(*m2))),
                tuple(list(itertools.chain(*m3))),
                tuple(list(itertools.chain(*m4))),
                tuple(list(itertools.chain(*hor_sym))),
                tuple(list(itertools.chain(*m6))),
                tuple(list(itertools.chain(*m7))),
                tuple(list(itertools.chain(*m8))))
        return set(evalu)
    def fuel_dynamics(self,presc_burn):
        post = (np.array(self.shape) + 1) * (1 - np.array(presc_burn))
        land_post = np.minimum(post, 2)
        self.fuel_dynamics=land_post
    def succession(self,constraint_biodiv):
        if len(self.fuel_dynamics)>1:
            biod =np.mat(np.diagonal((self.fuel_dynamics >= 1).dot(connectivity_mat(self.size)).dot(np.transpose((self.fuel_dynamics >= 1))) - (self.fuel_dynamics >= 1).sum(
                    axis=1))).T
            fuel = np.transpose(np.matrix(np.diagonal(
                (self.fuel_dynamics >= 2).dot(connectivity_mat(self.size)).dot(np.transpose((self.fuel_dynamics >= 2))) - (self.fuel_dynamics >= 2).sum(axis=1))))

            value = fuel + self.fuel + 1000*(biod<=constraint_biodiv)
            a = value.min(0)
            b = value.argmin(0)

            self.succession_land = {k:v for (k,v) in zip(constraint_biodiv, [np.array(x) for x in self.fuel_dynamics[np.asarray(b, int)][0].tolist()])}
            self.succession_value ={k:v for (k,v) in zip(constraint_biodiv, np.asarray(a)[0].tolist())}
            return self.succession_land, self.succession_value

        else :
            print("Problem")
    # Functions for characterization
    def test_algo(self,k):
        if k>=0:
            grid = (self.mat >= k)
        elif k==0:
            grid = (self.mat==k)
        grid = grid.astype(int)

        all_coordinates = []
        for i in range(4):
            all_coordinates.append(list(zip([i] * 4, range(4))))
        all_coordinates = [item for sublist in all_coordinates for item in sublist]

        paths = {}
        for i in all_coordinates:
            for j in all_coordinates:
                paths[(i,j)]=minDistance(grid, i[0], i[1], j[0], j[1])

        b = [key for key in list(paths.keys()) if paths[key] == 1000]

        for key in b:
            if paths[(key[1], key[0])] != 1000:
                paths[(key[0], key[1])] = paths[(key[1], key[0])]
        return(paths)
    def shortest_path(self,sourcex,sourcey,destx,desty):
        grid = (self.mat>=1)
        grid = grid.astype(int)
        a = minDistance(grid,sourcex,sourcey,destx,desty)
        b = minDistance(grid,destx,desty,sourcey,sourcex)
        if a!=b:
            return(min(a,b))
        else:
            return(a)

    def IIC_with_nc(self,var="biod"):
        if var=="biod":
            paths = self.test_algo(1)
        elif var=="fuel":
            paths = self.test_algo(2)
        invert = [1 / (1 + x) for x in list(paths.values())]
        # replace 1/1001 or 1/101 by 0
        iic = sum(invert) / (self.size) ** 4
        return(iic)
    def IIC_without_nc(self,var="biod"):
        if var=="biod":
            paths = self.test_algo(1)
        elif var=="fuel":
            paths = self.test_algo(2)
        screen1 = [x for x in list(paths.values()) if x<1000]
        invert = [1 / (1 + x) for x in screen1]
        # replace 1/1001 or 1/101 by 0
        iic = sum(invert) / (self.size) ** 4
        return(iic)
    def CPL(self,var="biod"):
        if var == "biod":
            paths = self.test_algo(1)
        elif var=="fuel":
            paths = self.test_algo(2)
        # need to remove all pairs such that (i=j) from paths
        paths2 = [ paths[k] for k in list(paths.keys()) if k[0]!=k[1]]
        return stats.mean(paths2)
    def CPL_without_nc(self,var="biod"):
        if var == "biod":
            paths = self.test_algo(1)
        elif var=="fuel":
            paths = self.test_algo(2)
        # need to remove all pairs such that (i=j) from paths
        paths2 = [ paths[k] for k in list(paths.keys()) if k[0]!=k[1]]
        paths2 = [x for x in paths2 if x<1000]
        try:
            result = stats.mean(paths2)
        except statistics.StatisticsError:
            result = float("nan")
        return result
    def diameter(self,var="biod"):
        if var=="zero":
            paths = self.test_algo(0)
        if var=="biod":
            paths = self.test_algo(1)
        elif var=="fuel":
            paths = self.test_algo(2)
        paths2 = [paths[x] for x in paths.keys() if paths[x]<1000]
        return max(paths2)
    def landscape_shape_index(self,var="biod"):
        if var=='biod':
            return 0.25*self.perimeter(var="biod")/self.size
        elif var=='fuel':
            return 0.25*self.perimeter(var="fuel")/self.size
    def diagonal(self,direction,var="biod"):
        if var=="biod":
            paths = self.test_algo(1)
        elif var == "fuel":
            paths = self.test_algo(2)

        if direction == 'NO-SE':
            return(paths[((0,0),(self.size-1,self.size-1))])
        elif direction=='SO-NE':
            return(paths[((size-1,0),(0,size-1))])
    # Description
    def number_empty_neighbors(self,i,j,var="biod"):
        size = self.size - 1
        if var=="biod":
            threshold=1
        elif var=="fuel":
            threshold=2

        if i>0:
            up = (self.mat[i-1][j]<threshold)
        elif i==0:
            up=1

        if i<size:
            down = (self.mat[i+1][j]<threshold)
        elif i==size:
            down=1

        if j>0:
            left = (self.mat[i][j-1]<threshold)
        elif j==0:
            left=1

        if j<size:
            right= (self.mat[i][j+1]<threshold)
        elif j==size:
            right=1

        return(sum([up,down,right,left]))
    def number_zero_adjacent(self,i,j, var="biod"):
        size = self.size - 1

        if var=="biod":
            threshold=1
        elif var=="fuel":
            threshold=2

        if i<0 or j<0 or i>size or j>size:
            return(0)

        if i>0:
            up = int(self.mat[i-1][j]==0)
        elif i<=0:
            up=1

        if i<size:
            down = int(self.mat[i+1][j]==0)
        elif i>=size:
            down=1
        if j>0:
            left = int(self.mat[i][j-1]==0)
        elif j<=0:
            left=1
        if j<size:
            right= int(self.mat[i][j+1]==0)
        elif j>=size:
            right=1

        if i>0 and j<size:
            up_right = int(self.mat[i-1][j+1]==0)
        else:
            up_right=1

        if i>0 and j>0:
            up_left = int(self.mat[i-1][j-1]==0)
        else:
            up_left =1

        if i<size and j>0:
            down_left = int(self.mat[i+1][j-1]==0)
        else:
            down_left=1

        if i<size and j<size:
            down_right = int(self.mat[i+1][j+1]==0)
        else:
            down_right=1

        sum_neighbor = sum([up,down,right,left,up_right,up_left,down_right,down_left])
        return(sum_neighbor*(self.mat[i][j]>=threshold))
    def nb_patches(self,var="biod"):
        size=self.size
        nb = []
        if var=="biod":
            for i in range(size):
                for j in range(size):
                    nb.append(self.number_zero_adjacent(i,j,var="biod"))
        elif var=="fuel":
            for i in range(size):
                for j in range(size):
                    nb.append(self.number_zero_adjacent(i,j,var="fuel"))

        return(1+sum([x==8 for x in nb]))
        # Components
    def components(self,var = "biod"):
        if var == "biod":
            check = (self.mat>=1).astype(int)
        elif var=="fuel":
            check = (self.mat>=2).astype(int)
        graph = list([list(x) for x in check])
        row = len(graph)
        col = len(graph[0])

        g = Graph(row, col, graph)

        return g.countIslands()[0]
    def components_area(self,var='biod'):
        if var == "biod":
            check = (self.mat>=1).astype(int)
        elif var=="fuel":
            check = (self.mat>=2).astype(int)
        graph = list([list(x) for x in check])
        row = len(graph)
        col = len(graph[0])

        g = Graph(row, col, graph)
        result = g.countIslands()[1]
        if len(result)==1:
            result = int(result[0])
        elif len(result)==0:
            result = float("nan")
        return result
    def components_area_max(self,var='biod'):
        if var == "biod":
            check = (self.mat>=1).astype(int)
        elif var=="fuel":
            check = (self.mat>=2).astype(int)
        graph = list([list(x) for x in check])
        row = len(graph)
        col = len(graph[0])

        g = Graph(row, col, graph)

        result = g.countIslands()[1]
        if len(result) == 1:
            result = int(result[0])
        elif len(result) == 0:
            result = float("nan")
        else:
            result = max(result)
        return result
    def components_perimeter(self, var= "biod",components_var=False):
        vars = var
        if var == "biod":
            check = (self.mat>=1).astype(int)
        elif var=="fuel":
            check = (self.mat>=2).astype(int)

        graph = list([list(x) for x in check])
        row = len(graph)
        col = len(graph[0])

        g = Graph(row,col,graph)
        first_nodes = g.countIslands()[2]

        if len(first_nodes) <= 1:
            if var=="biod":
                return(self.perimeter(var="biod"))
            elif var=='fuel':
                return(self.perimeter(var="fuel"))
        else:
            components = []
            perimeter_components = []
            for nodes in first_nodes:
                queue = []
                queue.append(nodes)
                visited = np.ones((4, 4), dtype=bool) == False
                patches = []
                cols = [0, 0, -1, 1, -1, -1, 1, 1]
                rows = [-1, 1, 0, 0, 1, -1, 1, -1]
                while len(queue) > 0:
                    node = queue.pop(0)

                    for i, j in list(zip(rows, cols)):
                        if isValid(node[0] + i, node[1] + j, check, visited):
                            queue.append((node[0] + i, node[1] + j))
                            visited[node[0] + i, node[1] + j] = True
                            patches.append((node[0] + i, node[1] + j))
                components.append(patches)
                perimeter = 0
                for patch in patches:
                    perimeter += check[patch[0], patch[1]]* self.number_empty_neighbors(patch[0], patch[1], var=vars)

                perimeter_components.append(perimeter)
            if components_var== True:
                return(perimeter_components,components)
            else:
                return (perimeter_components)
    def components_shape_index(self,var="biod"):
        if var =="biod":
            if self.components()>1:
                area = max(self.components_area(var="biod"))
                candidate = self.components_area(var="biod").index(area)
                perimeter = self.components_perimeter(var="biod")[candidate]
            else:
                area = self.components_area(var="biod")
                perimeter = self.components_perimeter(var="biod")

        elif var=="fuel":
            if self.components("fuel")>1:
                area = max(self.components_area(var="fuel"))
                candidate = self.components_area(var="fuel").index(area)
                perimeter = self.components_perimeter(var="fuel")[candidate]
            else:
                area = self.components_area(var="fuel")
                perimeter = self.components_perimeter(var="fuel")
        return(0.25*perimeter/math.sqrt(area))

        # Overall graph
    def node_degree(self,i,j,var="biod"):
        if i<0 or j<0 or j>self.size-1 or i>self.size-1:
            return(0)
        if var=="biod":
            if self.mat[i,j]==0:
                return(0)
            else:
                return 8- self.number_zero_adjacent(i,j,var='biod')
        elif var=="fuel":
            if self.mat[i,j]==0:
                return(0)
            else:
                return 8- self.number_zero_adjacent(i,j,var='fuel')
    def connectivity_correlation(self,var="biod",option_plot="No"):
        store_node = []
        store_neighbors = []
        for i in range(self.size):
            for j in range(self.size):
                first_coord = [-1,1,0,0,-1,-1,1,1]
                second_coord = [0,0,-1,1,-1,1,-1,1]
                store = []
                if var=="biod":
                    for nb in range(8):
                        store.append(self.node_degree(i+first_coord[nb],j+second_coord[nb],var="biod"))
                    store_node.append(self.node_degree(i,j,var="biod"))
                elif var=="fuel":
                    for nb in range(8):
                        store.append(self.node_degree(i+first_coord[nb],j+second_coord[nb],var="fuel"))
                    store_node.append(self.node_degree(i,j,var="fuel"))

                positive_x=[x for x in store if x!=0]
                try:
                    store_neighbors.append(sum(positive_x)/len(positive_x))
                except stats.StatisticsError:
                    store_neighbors.append(0)
        try:
            coefficient = np.corrcoef(np.array(store_node),np.array(store_neighbors))
        except RuntimeWarning:
            coefficient = float("nan")
        if option_plot=="Yes":
            plt.scatter(store_node,store_neighbors)
            plt.ylabel("Average # of edges for 8 neighbors ")
            plt.xlabel("Number of edges of node")
            plt.show()
        else:
            return coefficient[0,1]
    def perimeter(self,var="biod"):
        perimeter = 0
        if var=="biod":
            for j in range(len(self.mat[0])):
                for i in range(len(self.mat[0])):
                    perimeter+= (self.mat[i][j]>0)*self.number_empty_neighbors(i,j)
        elif var=="fuel":
            for j in range(len(self.mat[0])):
                for i in range(len(self.mat[0])):
                    perimeter += (self.mat[i][j] > 1) * self.number_empty_neighbors(i, j, var='fuel')
        return perimeter

    def perimeter_to_area(self,var='biod'):
        if var=='biod':
            return self.perimeter(var="biod")/math.sqrt(sum(sum(self.mat==1)))
        elif var=='fuel':
            return self.perimeter(var="fuel")/math.sqrt(sum(sum(self.mat==2)))
    def coord_zero(self):
        coord_zero=[]
        for i in range(self.size):
            for j in range(self.size):
                if self.mat[i,j]==0:
                    coord_zero.append((i,j))
        return coord_zero

class Land_lite:
    def __init__(self,shape):
        self.shape= shape
        self.size = int(math.sqrt(len(shape)))
        self.mat = to_matrix(shape,size)
        self.biod = (shape>=1).dot(connectivity_mat(self.size)).dot(np.transpose((shape>=1)))-sum(shape>=1)
        #self.biod_relative = self.biod/(np.repeat(1,len(self.shape)).dot(connectivity_mat(self.size)).dot(np.transpose(np.repeat(1,len(self.shape))))-sum(np.repeat(1,len(self.shape))))
        self.fuel = (shape>=2).dot(connectivity_mat(self.size)).dot(np.transpose((shape>=2)))-sum(shape>=2)
        #self.fuel_relative = self.fuel/(np.repeat(1,len(self.shape)).dot(connectivity_mat(self.size)).dot(np.transpose(np.repeat(1,len(self.shape))))-sum(np.repeat(1,len(self.shape))))
        self.node_biod = sum(shape>0)
        self.node_fuel = sum(shape==2)
        self.zero = sum(shape==0)
        self.habitat_area_ratio = sum(self.shape>=1)/(self.size**2)
        self.fuel_area_ratio = sum(self.shape==2)/(self.size**2)
        self.test_algo_0 = self.test_algo(0)
        self.test_algo_1 = self.test_algo(1)
        self.test_algo_2 = self.test_algo(2)
        self.components_biod_tot = self.components()
        self.components_fuel_tot = self.components('fuel')
    # Functions for characterization
    def test_algo(self,k):
        if k>=0:
            grid = (self.mat >= k)
        elif k==0:
            grid = (self.mat==k)
        grid = grid.astype(int)

        all_coordinates = []
        for i in range(4):
            all_coordinates.append(list(zip([i] * 4, range(4))))
        all_coordinates = [item for sublist in all_coordinates for item in sublist]

        paths = {}
        for i in all_coordinates:
            for j in all_coordinates:
                paths[(i,j)]=minDistance(grid, i[0], i[1], j[0], j[1])

        b = [key for key in list(paths.keys()) if paths[key] == 1000]

        for key in b:
            if paths[(key[1], key[0])] != 1000:
                paths[(key[0], key[1])] = paths[(key[1], key[0])]
        return(paths)
    def shortest_path(self,sourcex,sourcey,destx,desty):
        grid = (self.mat>=1)
        grid = grid.astype(int)
        a = minDistance(grid,sourcex,sourcey,destx,desty)
        b = minDistance(grid,destx,desty,sourcey,sourcex)
        if a!=b:
            return(min(a,b))
        else:
            return(a)

    def IIC_with_nc(self,var="biod"):
        if var=="biod":
            paths = self.test_algo_1
        elif var=="fuel":
            paths = self.test_algo_2
        invert = [1 / (1 + x) for x in list(paths.values())]
        # replace 1/1001 or 1/101 by 0
        iic = sum(invert) / (self.size) ** 4
        return(iic)
    def IIC_without_nc(self,var="biod"):
        if var=="biod":
            paths = self.test_algo_1
        elif var=="fuel":
            paths = self.test_algo_2
        screen1 = [x for x in list(paths.values()) if x<1000]
        invert = [1 / (1 + x) for x in screen1]
        # replace 1/1001 or 1/101 by 0
        iic = sum(invert) / (self.size) ** 4
        return(iic)
    def CPL(self,var="biod"):
        if var == "biod":
            paths = self.test_algo_1
        elif var=="fuel":
            paths = self.test_algo_2
        # need to remove all pairs such that (i=j) from paths
        paths2 = [ paths[k] for k in list(paths.keys()) if k[0]!=k[1]]
        return stats.mean(paths2)
    def CPL_without_nc(self,var="biod"):
        if var == "biod":
            paths = self.test_algo_1
        elif var=="fuel":
            paths = self.test_algo_2
        # need to remove all pairs such that (i=j) from paths
        paths2 = [ paths[k] for k in list(paths.keys()) if k[0]!=k[1]]
        paths2 = [x for x in paths2 if x<1000]
        try:
            result = stats.mean(paths2)
        except statistics.StatisticsError:
            result = float("nan")
        return result

    def diameter(self,var="biod"):
        if var=="zero":
            paths = self.test_algo_0
        if var=="biod":
            paths = self.test_algo_1
        elif var=="fuel":
            paths = self.test_algo_2
        paths2 = [paths[x] for x in paths.keys() if paths[x]<1000]
        return max(paths2)


    def landscape_shape_index(self,var="biod"):
        if var=='biod':
            return 0.25*self.perimeter(var="biod")/self.size
        elif var=='fuel':
            return 0.25*self.perimeter(var="fuel")/self.size
    def landscape_perimeter_area(self,var="biod"):
        if var=='biod':
            return 0.25*self.perimeter(var="biod")/sum(self.shape>=1)
        elif var=='fuel':
            return 0.25*self.perimeter(var="fuel")/sum(self.shape==2)

    # Description
    # Utility functions : empty neighbors and zero adjacent
    def number_empty_neighbors(self,i,j,var="biod"):
        """
        4 direction empty neighbors
        :param i: Coordinate in grid
        :param j: Coordinate in grid
        :param var: "biod" or "fuel" depending on which metric is evaluated
        :return: number of empty 4 neighbors
        """
        size = self.size - 1
        if var=="biod":
            threshold=1
        elif var=="fuel":
            threshold=2

        if i>0:
            up = (self.mat[i-1][j]<threshold)
        elif i==0:
            up=1

        if i<size:
            down = (self.mat[i+1][j]<threshold)
        elif i==size:
            down=1

        if j>0:
            left = (self.mat[i][j-1]<threshold)
        elif j==0:
            left=1

        if j<size:
            right= (self.mat[i][j+1]<threshold)
        elif j==size:
            right=1

        return(sum([up,down,right,left]))
    def number_zero_adjacent(self,i,j, var="biod"):
        """
        Number of 8 neighbors empty neighbors of cell (i,j) in grid
        :param i: Coordinate in grid
        :param j: Coordinate in grid
        :param var: "biod" or "fuel" depending on which metric is evaluated
        :return:
        """
        size = self.size - 1

        if var=="biod":
            threshold=1
        elif var=="fuel":
            threshold=2

        if i<0 or j<0 or i>size or j>size:
            return(0)

        if i>0:
            up = int(self.mat[i-1][j]==0)
        elif i<=0:
            up=1

        if i<size:
            down = int(self.mat[i+1][j]==0)
        elif i>=size:
            down=1
        if j>0:
            left = int(self.mat[i][j-1]==0)
        elif j<=0:
            left=1
        if j<size:
            right= int(self.mat[i][j+1]==0)
        elif j>=size:
            right=1

        if i>0 and j<size:
            up_right = int(self.mat[i-1][j+1]==0)
        else:
            up_right=1

        if i>0 and j>0:
            up_left = int(self.mat[i-1][j-1]==0)
        else:
            up_left =1

        if i<size and j>0:
            down_left = int(self.mat[i+1][j-1]==0)
        else:
            down_left=1

        if i<size and j<size:
            down_right = int(self.mat[i+1][j+1]==0)
        else:
            down_right=1

        sum_neighbor = sum([up,down,right,left,up_right,up_left,down_right,down_left])
        return(sum_neighbor*(self.mat[i][j]>=threshold))
    def nb_patches(self,var="biod"):
        size=self.size
        nb = []
        if var=="biod":
            for i in range(size):
                for j in range(size):
                    nb.append(self.number_zero_adjacent(i,j,var="biod"))
        elif var=="fuel":
            for i in range(size):
                for j in range(size):
                    nb.append(self.number_zero_adjacent(i,j,var="fuel"))

        return(1+sum([x==8 for x in nb]))
        # Components
    def components(self,var = "biod"):
        if var == "biod":
            check = (self.mat>=1).astype(int)
        elif var=="fuel":
            check = (self.mat>=2).astype(int)
        graph = list([list(x) for x in check])
        row = len(graph)
        col = len(graph[0])

        g = Graph(row, col, graph)

        return g.countIslands()

    def components_perimeter(self, var= "biod",components_var=False):
        vars = var
        if var == "biod":
            first_nodes = self.components_biod_tot[2]
            check = (self.mat>=1).astype(int)
        elif var=="fuel":
            first_nodes = self.components_fuel_tot[2]
            check = (self.mat==2).astype(int)
        if len(first_nodes) <= 1:
            if var=="biod":
                return(self.perimeter(var="biod"))
            elif var=='fuel':
                return(self.perimeter(var="fuel"))
        else:
            components = []
            perimeter_components = []
            for nodes in first_nodes:
                queue = []
                queue.append(nodes)
                visited = np.ones((4, 4), dtype=bool) == False
                patches = []
                cols = [0, 0, -1, 1, -1, -1, 1, 1]
                rows = [-1, 1, 0, 0, 1, -1, 1, -1]
                while len(queue) > 0:
                    node = queue.pop(0)

                    for i, j in list(zip(rows, cols)):
                        if isValid(node[0] + i, node[1] + j, check, visited):
                            queue.append((node[0] + i, node[1] + j))
                            visited[node[0] + i, node[1] + j] = True
                            patches.append((node[0] + i, node[1] + j))
                components.append(patches)
                perimeter = 0
                for patch in patches:
                    perimeter += check[patch[0], patch[1]]* self.number_empty_neighbors(patch[0], patch[1], var=vars)

                perimeter_components.append(perimeter)
            if components_var== True:
                return(perimeter_components,components)
            else:
                return (perimeter_components)
    def components_shape_index(self,var="biod"):
        if var =="biod":
            if self.components_biod_tot[0]>1:
                area = max(self.components_biod_tot[1])
                candidate = self.components_biod_tot[1].index(area)
                try:
                    perimeter = self.components_perimeter(var="biod")[candidate]
                except IndexError:
                    perimeter =  self.components_perimeter(var="biod")
            else:
                area = self.components_biod_tot[1][0]
                perimeter = self.components_perimeter(var="biod")

        elif var=="fuel":
            if self.components_fuel_tot[0]>1:
                area = max(self.components_fuel_tot[1])
                candidate = self.components_fuel_tot[1].index(area)
                try:
                    perimeter = self.components_perimeter(var="fuel")[candidate]
                except IndexError:
                    perimeter =  self.components_perimeter(var="fuel")
            else:
                area = self.components_fuel_tot[1][0]
                perimeter = self.components_perimeter(var="fuel")
        return(0.25*perimeter/math.sqrt(area))

        # Overall graph

    def node_degree(self,i,j,var="biod"):
        if i<0 or j<0 or j>self.size-1 or i>self.size-1:
            return(0)
        if var=="biod":
            if self.mat[i,j]==0:
                return(0)
            else:
                return 8- self.number_zero_adjacent(i,j,var='biod')
        elif var=="fuel":
            if self.mat[i,j]==0:
                return(0)
            else:
                return 8- self.number_zero_adjacent(i,j,var='fuel')
    def connectivity_correlation(self,var="biod",option_plot="No"):
        store_node = []
        store_neighbors = []
        for i in range(self.size):
            for j in range(self.size):
                first_coord = [-1,1,0,0,-1,-1,1,1]
                second_coord = [0,0,-1,1,-1,1,-1,1]
                store = []
                if var=="biod":
                    for nb in range(8):
                        store.append(self.node_degree(i+first_coord[nb],j+second_coord[nb],var="biod"))
                    store_node.append(self.node_degree(i,j,var="biod"))
                elif var=="fuel":
                    for nb in range(8):
                        store.append(self.node_degree(i+first_coord[nb],j+second_coord[nb],var="fuel"))
                    store_node.append(self.node_degree(i,j,var="fuel"))

                positive_x=[x for x in store if x!=0]
                try:
                    store_neighbors.append(sum(positive_x)/len(positive_x))
                except stats.StatisticsError:
                    store_neighbors.append(0)
                except ZeroDivisionError:
                    store_neighbors.append(0)

        try:
            coefficient = np.corrcoef(np.array(store_node),np.array(store_neighbors))
        except RuntimeWarning:
            coefficient = float("nan")
        if option_plot=="Yes":
            plt.scatter(store_node,store_neighbors)
            plt.ylabel("Average # of edges for 8 neighbors ")
            plt.xlabel("Number of edges of node")
            plt.show()
        else:
            return coefficient[0,1]
    def perimeter(self,var="biod"):
        perimeter = 0
        if var=="biod":
            for j in range(len(self.mat[0])):
                for i in range(len(self.mat[0])):
                    perimeter+= (self.mat[i][j]>0)*self.number_empty_neighbors(i,j)
        elif var=="fuel":
            for j in range(len(self.mat[0])):
                for i in range(len(self.mat[0])):
                    perimeter += (self.mat[i][j] > 1) * self.number_empty_neighbors(i, j, var='fuel')
        return perimeter

    def perimeter_to_area(self,var='biod'):
        if var=='biod':
            return self.perimeter(var="biod")/math.sqrt(sum(sum(self.mat==1)))
        elif var=='fuel':
            return self.perimeter(var="fuel")/math.sqrt(sum(sum(self.mat==2)))

def convergence(successions):
    """
    Returns the time of convergence and length of cycle if applicable
    :param successions: data series
    :return: a list of tuples (t_convergence, cycle_period)
    """
    lands = [Land(shape=np.array(ast.literal_eval(x))) for x in successions]
    horizon = list(range(len(successions)))
    horizon.remove(0)

    for t in horizon:
        if len(lands[t].equivalent().intersection(set([tuple(lands[t - 1].shape)]))) > 0:
            convergence_time=t
            break

    for t in range(convergence_time+ 1, horizon[-1], 1):
        if tuple(lands[t].shape) == tuple(lands[convergence_time].shape):
            cycle = t - convergence_time
            break
    return convergence_time,cycle



def statistics_dataframe(data_path):
    index = re.findall(r'\d+',data_path)
    budget_here=index[2]
    path="C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_"+str(budget_here)+"/successions_matched/"
    data = pd.read_csv(path+data_path).drop(columns=['Unnamed: 0'])

    data = data.iloc[:, 0:20]
    names = list(range(20))
    data.columns = names

    data["index"] = list(range(len(data)))
    # data['index_biod'] = [12]*len(data)
    # could be better to do in the end, with initial numbers as well.

    checker = data.melt(id_vars=["index"])
    checker["variable"] = checker['variable'].astype(int)
    checker = checker.sort_values(["index", "variable"])
    checker = checker.reset_index(drop=True)

    values = list(checker.value)
    values_unique = list(set(values))

    # region Setting up of data storage dictionnaries
    nodes_biod = {}
    nodes_fuel = {}

    biod_score = {}
    fuel_score = {}

    land_perimeter_biod = {}
    land_perimeter_fuel = {}

    land_habitat_area_ratio = {}
    land_fuel_area_ratio = {}

    land_shape_index_biod = {}
    land_shape_index_fuel = {}

    land_diameter_biod = {}
    land_diameter_fuel = {}

    land_connectivity_correlation_biod = {}
    land_connectivity_correlation_fuel = {}

    land_IIC_nc_biod = {}
    land_IIC_nc_fuel = {}

    land_IIC_biod = {}
    land_IIC_fuel = {}

    land_CPL_nc_biod = {}
    land_CPL_nc_fuel = {}

    land_CPL_biod = {}
    land_CPL_fuel = {}

    components_number_biod = {}
    components_number_fuel = {}

    components_area_biod_max = {}
    components_area_fuel_max = {}

    components_shape_index_biod = {}
    components_shape_index_fuel = {}
    # endregion

    # region Compute statistics
    debut = time.time()
    for shape in values_unique:
        land = Land_lite(np.array(ast.literal_eval(shape)))

        nodes_biod[shape] = land.node_biod
        nodes_fuel[shape] = land.node_fuel

        biod_score[shape] = land.biod
        fuel_score[shape] = land.fuel

        land_perimeter_biod[shape] = land.perimeter()
        land_perimeter_fuel[shape] = land.perimeter(var="fuel")

        land_habitat_area_ratio[shape] = land.habitat_area_ratio
        land_fuel_area_ratio[shape] = land.fuel_area_ratio

        land_shape_index_biod[shape] = land.landscape_shape_index()
        land_shape_index_fuel[shape] = land.landscape_shape_index(var="fuel")

        land_diameter_biod[shape] = land.diameter()
        land_diameter_fuel[shape] = land.diameter(var='fuel')

        land_connectivity_correlation_biod[shape] = land.connectivity_correlation()
        land_connectivity_correlation_fuel[shape] = land.connectivity_correlation(var="fuel")

        land_IIC_nc_biod[shape] = land.IIC_with_nc()
        land_IIC_nc_fuel[shape] = land.IIC_with_nc(var="fuel")

        land_IIC_biod[shape] = land.IIC_without_nc()
        land_IIC_fuel[shape] = land.IIC_without_nc(var="fuel")

        land_CPL_nc_biod[shape] = land.CPL()
        land_CPL_nc_fuel[shape] = land.CPL(var='fuel')

        land_CPL_biod[shape] = land.CPL_without_nc()
        land_CPL_fuel[shape] = land.CPL_without_nc(var="fuel")

        components_number_biod[shape] = land.components_biod_tot[0]
        components_number_fuel[shape] = land.components_fuel_tot[0]

        components_area_biod_max[shape] = max(land.components_biod_tot[1])
        components_area_fuel_max[shape] = max(land.components_fuel_tot[1])

        components_shape_index_biod[shape] = land.components_shape_index()
        components_shape_index_fuel[shape] = land.components_shape_index("fuel")
        # print(values_unique.index(shape)/len(values_unique))

    checker["nodes_biod"] = [nodes_biod.get(n, n) for n in values]
    checker["nodes_fuel"] = [nodes_fuel.get(n, n) for n in values]
    checker["nodes_zero"] = [16 - checker.nodes_biod.loc[x] for x in range(len(checker))]
    checker["score_biod"] = [biod_score.get(n, n) for n in values]
    checker["score_fuel"] = [fuel_score.get(n, n) for n in values]

    checker["land_perimeter_biod"] = [land_perimeter_biod.get(n, n) for n in values]
    checker['land_perimeter_fuel'] = [land_perimeter_fuel.get(n, n) for n in values]
    checker["land_shape_index_biod"] = [land_shape_index_biod.get(n, n) for n in values]
    checker["land_shape_index_fuel"] = [land_shape_index_fuel.get(n, n) for n in values]

    checker["land_diameter_biod"] = [land_diameter_biod.get(n, n) for n in values]
    checker["land_diameter_fuel"] = [land_diameter_fuel.get(n, n) for n in values]

    checker["land_connectivity_correlation_biod"] = [land_connectivity_correlation_biod.get(n, n) for n in values]
    checker["land_connectivity_correlation_fuel"] = [land_connectivity_correlation_fuel.get(n, n) for n in values]

    checker["land_IIC_nc_biod"] = [land_IIC_nc_biod.get(n, n) for n in values]
    checker["land_IIC_nc_fuel"] = [land_IIC_nc_fuel.get(n, n) for n in values]

    checker["land_IIC_biod"] = [land_IIC_biod.get(n, n) for n in values]
    checker["land_IIC_fuel"] = [land_IIC_fuel.get(n, n) for n in values]

    checker["land_CPL_nc_biod"] = [land_CPL_nc_biod.get(n, n) for n in values]
    checker["land_CPL_nc_fuel"] = [land_CPL_nc_fuel.get(n, n) for n in values]

    checker["land_CPL_biod"] = [land_CPL_biod.get(n, n) for n in values]
    checker["land_CPL_fuel"] = [land_CPL_fuel.get(n, n) for n in values]

    checker["components_number_biod"] = [components_number_biod.get(n, n) for n in values]
    checker["components_number_fuel"] = [components_number_fuel.get(n, n) for n in values]

    checker["components_area_max_biod"] = [components_area_biod_max.get(n, n) for n in values]
    checker["components_area_max_fuel"] = [components_area_fuel_max.get(n, n) for n in values]

    checker["components_shape_index_biod"] = [components_shape_index_biod.get(n, n) for n in values]
    checker["components_shape_index_fuel"] = [components_shape_index_fuel.get(n, n) for n in values]

    convergence_time = []
    convergence_period = []
    indices = list(range(0, len(checker) + 1, 20))
    indices.remove(0)
    for x in indices:
        succession = checker.value.loc[x - 20:x - 1]
        result = convergence(succession)
        convergence_time_s = result[0]
        convergence_period_s = result[1]
        convergence_time.extend([convergence_time_s] * 20)
        convergence_period.extend([convergence_period_s] * 20)
    del result

    checker['convergence_time'] = convergence_time
    checker['convergence_period'] = convergence_period

    checker['index_biod'] = [index[-1]] * len(checker)
    checker['initial_nb_1'] = [index[3]-index[4]] * len(checker)
    checker['initial_nb_2'] = [index[4]] * len(checker)
    checker['budget']=[budget_here]*len(checker)

    checker.to_csv("C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results/budget_"+str(budget_here)+"/statistics/stats_"+data_path)
    return(checker)
# checking where move is valid or not
