#region Packages
import numpy as np
import pandas as pd
import time
import ast
import csv
import random as rd
import re as re
import math
import multiprocessing as mp
import os
import json
import pickle
import itertools
#endregion

#region Define parameters
# region Parameters that can vary


budget = 1
#endregion
#region Underlying parameters


# Values
size    = 4
R       = size**2
# paths
path = "/data/public_data/connectivity_jean/budget_" + str(budget) + "/"
files = os.listdir(path)
files.remove('keys_succession')
files.remove('successions_matched')

# Biodiversity and fire thresholds
d_seuil = 2
d       = d_seuil*np.ones(R)
m_seuil = 1
m       = m_seuil*np.ones(R)
A_max   = 2
nb_age  = A_max+1
biodiv  = np.arange(2,84,2)
biodiv3 = np.arange(0,50,2)


#endregion

#endregion

#region Functions to generate landscapes
def cost(size):
    if size==4:
        return(np.array([[1,1,1,1],
                        [1,1,1,1],
                        [1,1,1,1],
                        [1,1,1,1]]))
    elif size==3:
        return(np.array([[1, 1, 1],
                         [1, 1, 1],
                         [1, 1, 1]]))
    elif size==5:
        return (np.array([[1, 1, 1, 1, 1],
                          [1, 1, 1, 1, 1],
                          [1, 1, 1, 1, 1],
                          [1, 1, 1, 1, 1],
                          [1, 1, 1, 1, 1]]))

costs    = cost(size)

def nonzero(land):
    return(sum(land>0))
def nb2(land):
    return(sum(land==2))
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
    return(np.array([list(essai[0:3]),
              list(essai[3:6]),
              list(essai[6:9])]))
def equi_landscape3(essai):
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

            mat_test = equi_landscape3(to_matrix3(candidate))

            if len(mat_test.intersection(unique_landscapes_nb2)) == 0:
                unique_landscapes_nb2.add(tuple(candidate))
    return(unique_landscapes_nb2)
#endregion
#region Landscape 4x4
def to_matrix4(essai):
    essai2 = np.array([list(essai[0:4]),
                       list(essai[4:8]),
                       list(essai[8:12]),
                       list(essai[12:16])],dtype=object)
    return(essai2)
def equi_landscape4(essai):
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

            mat_test = equi_landscape4(to_matrix4(candidate))

            if len(mat_test.intersection(unique_landscapes_nb2))== 0:
                unique_landscapes_nb2.add(tuple(candidate))

    return(unique_landscapes_nb2)

#endregion
#region Landscape 5x5
def to_matrix5(essai):
    essai2 = np.array([list(essai[0:5]),
                       list(essai[5:10]),
                       list(essai[10:15]),
                       list(essai[15:20]),
                       list(essai[20:25])],dtype=int)
    return(essai2)
def equi_landscape5(essai):
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
    return (np.array(landscape) >= d)
def high_fuel_connect(landscape,size=4):
    return(high_fuel(landscape).dot(connectivity_mat(size)).dot(high_fuel(landscape)[:,np.newaxis])-sum(high_fuel(landscape)))
def high_biodiv(landscape):
    return (np.array(landscape) >= m_seuil)
def high_fuel_con_vect(land, size=4):
    return(np.transpose(np.matrix(np.diagonal((land>=d).dot(connectivity_mat(size)).dot(np.transpose((land>=d)))-(land>=d).sum(axis=1)))))
def high_biod_con_vect(new_land, size=4):
   a= np.mat(np.diagonal((new_land>=m).dot(connectivity_mat(size)).dot(np.transpose((new_land>=m)))-(new_land>=m).sum(axis=1))).T
   return(a)

def connectivity_mat3(size=3):
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
    return(np.transpose(np.matrix(np.diagonal((land>=d).dot(connectivity_mat(size)).dot(np.transpose((land>=d)))-(land>=d).sum(axis=1)))))
def high_fuel_connect3(landscape,size=3):
    return(high_fuel(landscape).dot(connectivity_mat(size)).dot(high_fuel(landscape)[:,np.newaxis])-sum(high_fuel(landscape)))
def high_biod_con_vect3(new_land, size=3):
   a= np.mat(np.diagonal((new_land>=m).dot(connectivity_mat(size)).dot(np.transpose((new_land>=m)))-(new_land>=m).sum(axis=1))).T
   return(a)
#endregion

#region Functions for dynamic programming
#region 3x3 landscape

#endregion


#region Functions for successions
def equi_landscape4_recover(essai):
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
    land = to_matrix4(land)
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
        m7 = np.rot90(np.rot90([land[3],land[2],land[1],land[0]],axes=(1,0)),axes=(1,0))
        return(tuple(list(itertools.chain(*m7))))
    if transfo ==6:
        m8 = np.rot90(np.rot90(np.rot90([land[3],land[2],land[1],land[0]],axes=(1,0)),axes=(1,0)),axes=(1,0))
        return(tuple(list(itertools.chain(*m8))))
    if transfo ==7:
        return(tuple(list(itertools.chain(*land))))
    if transfo=="None":
        return(tuple(list(itertools.chain(*land))))
def recover3(land, transfo):
    land = to_matrix3(land)
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
        m7 = np.rot90(np.rot90([land[2],land[1],land[0]],axes=(1,0)),axes=(1,0))
        return(tuple(list(itertools.chain(*m7))))
    if transfo ==6:
        m8 = np.rot90(np.rot90(np.rot90([land[2],land[1],land[0]],axes=(1,0)),axes=(1,0)),axes=(1,0))
        return(tuple(list(itertools.chain(*m8))))
    if transfo ==7:
        return(tuple(list(itertools.chain(*land))))
    if transfo=="None":
        return(tuple(list(itertools.chain(*land))))

def recover_full(land,transfo,index):
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
def succession_dataframe(data_path, biodiv_index,time_param='Yes',budget=1,follow_up="No"):

    debut = time.time()
    #region Load the current data

    file = open("/data/public_data/connectivity_jean/budget_"+ str(budget) + "/"+data_path)
    csvreader = csv.reader(file)
    header = next(csvreader)

    rows_data = []
    for row in csvreader:
        row_add = [row[0],row[int(biodiv_index/2)], row[int(82+biodiv_index/2)]]
        rows_data.append(row_add)
    #endregion
    #region Import data to scan through - dict format and list format
    path = "/data/public_data/connectivity_jean/budget_" + str(budget) + "/"
    files = os.listdir(path)
    files.remove('keys_succession')
    files.remove('successions_matched')
    files = [x for x in files if
             ast.literal_eval(re.findall(r'\d+', str(x))[2]) in list(range(16 - budget, 17))]
    datas = {}
    datas[data_path]=rows_data
    for x in files:
        file = open("/data/public_data/connectivity_jean/budget_" + str(budget) + "/" +x)
        csvreader = csv.reader(file)
        header = next(csvreader)

        rows_dict = []
        for row in csvreader:
            row_add = [row[0],row[int(biodiv_index/2)], row[int(82+biodiv_index/2)]]
            rows_dict.append(row_add)
        datas[x] = rows_dict
    # endregion

    # region Import equivalent data - dict format and list formats
    path = "/data/public_data/connectivity_jean/budget_" + str(budget) + "/keys_succession/"
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
        file = open("/data/public_data/connectivity_jean/budget_" + str(budget) + "/keys_succession/equiv_results_"+str(match_equiv[2])+"_"+str(match_equiv[3])+"_cut_"+str(match_equiv[4])+".csv")
        csvreader = csv.reader(file)
        header = next(csvreader)

        rows_dict = []
        for row in csvreader:
            row_add = row[int(biodiv_index/2-1)]
            rows_dict.append(row_add)
        datas_equiv[('equiv_results_'+str(match_equiv[2])+"_"+str(match_equiv[3])+"_cut_"+str(match_equiv[4])+'.csv')]=rows_dict
    else:
        file = open("/data/public_data/connectivity_jean/budget_" + str(budget) + "/keys_succession/equiv_results_"+str(match_equiv[2])+"_"+str(match_equiv[3])+".csv")
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
        path = "/data/public_data/connectivity_jean/budget_" + str(budget) + "/"
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
                    path = "/data/public_data/connectivity_jean/budget_" + str(
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
            mat_succession = to_matrix4(ast.literal_eval(succession_key))
            transfo = equi_landscape4_recover(mat_succession).index(ast.literal_eval(candidate_succ))
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
#endregion


#region Potential fires

pot_fire_value  = (0,1)
pot_fire        = list(itertools.product(pot_fire_value,repeat=R))
pot_fire_budget = [x for x in pot_fire if sum(sum(to_matrix(x,size)*costs))<=budget]
#endregion
class Land:
    def __init__(self,shape):
        self.shape= shape
        self.size = int(math.sqrt(len(shape)))
        self.mat = to_matrix(shape,size)
        self.biod = (shape>=1).dot(connectivity_mat(size)).dot(np.transpose((shape>=1)))-sum(shape>=1)
        self.fuel = (shape>=2).dot(connectivity_mat(size)).dot(np.transpose((shape>=2)))-sum(shape>=2)
        self.nonz = sum(shape>0)
        self.nb2 = sum(shape==2)
        self.zero = sum(shape==0)
        self.habitat_area_ratio = sum(self.shape>=1)/(self.size**2)


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
            print("prout")

    def shortest_path(self,sourcex,sourcey,destx,desty):
        grid = (self.mat>=1)
        grid = grid.astype(int)
        a = minDistance(grid,sourcex,sourcey,destx,desty)
        b = minDistance(grid,destx,desty,sourcey,sourcex)
        if a!=b:
            return(min(a,b))
        else:
            return(a)


    def IIC_fuel(self):
        paths = self.test_algo(2)
        invert = [1/(1+x) for x in list(paths.values())]
        iic = sum(invert)/(self.size)**4
        return(iic)
    def IIC_biod(self):
        paths = self.test_algo(1)

        invert = [1/(1+x) for x in list(paths.values())]
        iic = sum(invert)/(self.size)**4
        return(iic)
    def characteristic_path_length_biod(self):
        paths = self.test_algo(1)
        # need to remove all pairs such that (i=j) from paths
        paths2 = [ paths[k] for k in list(paths.keys()) if k[0]!=k[1]]
        return stats.mean(paths2)
    def characteristic_path_length_fire(self):
        paths = self.test_algo(2)
        paths2 = [ paths[k] for k in list(paths.keys()) if k[0]!=k[1]]


        return stats.mean(paths2)


    def test_algo(self,k):

        grid = (self.mat >= k)
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

    def diameter_biod(self):
        paths = self.test_algo(1)

        paths2 = [paths[x] for x in paths.keys() if paths[x]<1000]
        return max(paths2)
    def diameter_fire(self):
        paths = self.test_algo(2)

        paths2 = [paths[x] for x in paths.keys() if paths[x]<1000]
        return max(paths2)

    def diagonal_biod(self,direction):
        paths = self.test_algo(1)
        if direction == 'NO-SE':
            return(paths[((0,0),(self.size-1,self.size-1))])
        elif direction=='SO-NE':
            return(paths[((size-1,0),(0,size-1))])

    def adjacency_biod(self):
        paths = self.test_algo(1)

def final_approach(file,index):
    debut2 = time.time()
    path = "/home/jean/apprentissage/data/"
    filename = open(path+file, 'rb')

    list_nz = pickle.load(filename)
    filename.close()
    list_nz = list(list_nz)
    list_nz = [np.array(x) for x in list_nz]

    rows = []
    for land in list_nz:
        landscape = Land(shape=land)
        landscape.fuel_dynamics(pot_fire_budget)
        landscape.succession([index])
        storage_succession = [landscape.succession([index])[0][index]]
        storage_value = [landscape.succession([index])[1][index]]
        for t in list(range(20)):
            landscape_next = Land(shape=storage_succession[-1])
            landscape_next.fuel_dynamics(pot_fire_budget)
            storage_succession.append(landscape_next.succession([index])[0][index])
            storage_value.append(landscape_next.succession([index])[1][index])

        storage_succession.extend(storage_value)
        rows.append(storage_succession)

        file2 = file.replace(".pkl",'')
    with open('/home/jean/apprentissage/results/budget_'+str(budget)+'/'+file2+"_biod_index_" + str(
            index) + '.csv', 'w', encoding='UTF8', newline='') as f:
        writer = csv.writer(f)
        # write multiple rows
        writer.writerows(rows)
        return("Process " + file +" for index "+str(index)+" took "+ str(time.time()-debut2)+" seconds")


def match_heavy_b(k):
    file = open(path+files[k])
    nums = re.findall(r'\d+',files[k])
    csvreader=csv.reader(file)
    header=next(csvreader)

    rows=[]
    for row in csvreader:
        row_add = row[1:42]
        rows.append(row_add)
    landscapes = [item for sublist in rows for item in sublist]

    singular_land = [np.array(ast.literal_eval(i),int) for i in list(set(landscapes))]
#endregion

# Make sublists

    #def nonzero(land):
    #    return(sum(land>0))
    #def nb2(land):
    #    return(sum(land==2))

    path2 = "/data/public_data/connectivity_jean/all_data"
    files_data = os.listdir(path2)

    twos = set([nb2(i) for i in singular_land])
    nonzeros = set([nonzero(i) for i in singular_land])

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
            filename = open("/data/public_data/connectivity_jean/all_data/land4_" + str(j[0]) + "_" + str(
                j[1]) + ".pkl", 'rb')
            set_to_check = pickle.load(filename)
            filename.close()
            value = []
            to_rep=[]
            for v in store2[j]:
                set_verif = equi_landscape4(to_matrix4(v))
                set_verif.add(tuple(v))
                list_verif = [*set_verif.intersection(set_to_check)]
            #new = list(list_verif[0])
                replacement[str(tuple(v))] = str(list_verif[0])

        elif os.path.exists(path2+"/land4_" + str(j[0]) + "_" + str(j[1]) + "cut_0.pkl"):
            set_to_check_merg = []
            to_load = [x for x in files_data if x.startswith('land4_' + str(j[0]) + "_" + str(j[1]))]

            for x in to_load:
                filename = open("/data/public_data/connectivity_jean/all_data/"+x, 'rb')
                set_to_check = pickle.load(filename)
                filename.close()
                set_to_check =list(map(tuple, set_to_check))
                set_to_check_merg.append(set_to_check)

            set_to_check_merg = [item for sublist in set_to_check_merg for item in sublist]
            set_to_check_merg=set(set_to_check_merg)

            for v in store2[j]:
                set_verif = equi_landscape4(to_matrix4(v))
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
