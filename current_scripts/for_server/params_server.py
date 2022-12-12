"""
Module for parameter set-up for Mouysset & Jean, 2023.

Parameters :
-----------
- Budget
- Landscape size (n) in nxn matrix
- Biodiversity values range for various cases
- Potential prescribed burns given budget.

Functions :
-------------------
- to_matrix returns matrix of landscape
- costs : returns homogeneous cost grid

"""

#region Underlying parameters
import os
import numpy as np
import itertools

# Values
budget = 1
size = 4
R = size**2
# paths
if size == 4:
    path = "/home/simonjean/data/budget_" + str(budget) + "/"
    files = os.listdir(path)
    files.remove('keys_succession')
    files.remove('successions_matched')


# Biodiversity and fire thresholds
d_seuil = 2
d = d_seuil*np.ones(R)
m_seuil = 1
m = m_seuil*np.ones(R)
A_max = 2
nb_age = A_max+1
biodiv = np.arange(2, 100, 2)
biodiv3 = np.arange(0, 50, 2)

def to_matrix(land):
    """
    This function translates the array of length size into a sqrt(size)xsqrt(size) matrix
    :param land: np.array or tuple
    :return: matrix
    """
    if size == 3:
        return(np.array([list(land[0:3]),
                         list(land[3:6]),
                         list(land[6:9])]))
    elif size == 4:
        essai2 = np.array([list(land[0:4]),
                           list(land[4:8]),
                           list(land[8:12]),
                           list(land[12:16])])
        return essai2
    elif size == 5:
        essai2 = np.array([list(land[0:5]),
                           list(land[5:10]),
                           list(land[10:15]),
                           list(land[15:20]),
                           list(land[20:25])])
        return essai2
def cost():
    if size == 4:
        return(np.array([[1, 1, 1, 1],
                        [1, 1, 1, 1],
                        [1, 1, 1, 1],
                        [1, 1, 1, 1]]))
    elif size == 3:
        return(np.array([[1, 1, 1],
                         [1, 1, 1],
                         [1, 1, 1]]))
    elif size == 5:
        return (np.array([[1, 1, 1, 1, 1],
                          [1, 1, 1, 1, 1],
                          [1, 1, 1, 1, 1],
                          [1, 1, 1, 1, 1],
                          [1, 1, 1, 1, 1]]))

costs = cost()

pot_fire_value = (0, 1)
pot_fire = list(itertools.product(pot_fire_value, repeat=R))
pot_fire_budget = [x for x in pot_fire if sum(sum(to_matrix(x)*costs)) <= budget]

#endregion