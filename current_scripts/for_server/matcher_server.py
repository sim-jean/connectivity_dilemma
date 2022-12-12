import modules.utilities as utilities
import modules.params as params
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
import os
import json
import pickle
import itertools
import matplotlib.pyplot as plt
import seaborn as sns
#endregion

def nonzero(land):
    return sum(land > 0)

def nb2(land):
    return sum(land == 2)

def match_heavy_b(k):
    '''
    Matches results to their equivalent landscape in the set of unique landscapes for 4x4 landscapes
    :param k: int
    :return: dataset in .csv
    '''
    file = open(params.path+params.files[k])
    nums = re.findall(r'\d+', params.files[k])
    csvreader = csv.reader(file)
    header = next(csvreader)

    rows = []
    for row in csvreader:
        row_add = row[1:42]
        rows.append(row_add)
    landscapes = [item for sublist in rows for item in sublist]

    singular_land = [np.array(ast.literal_eval(i), int) for i in list(set(landscapes))]


    path2 = "/home/simonjean/data/all_data"
    files_data = os.listdir(path2)

    twos = set([nb2(i) for i in singular_land])
    nonzeros = set([nonzero(i) for i in singular_land])

    possible = list(itertools.product(nonzeros, twos))

    store = {}
    for i in range(len(possible)):
        store[possible[i]] = [j for j in singular_land if (nonzero(j) == possible[i][0] and nb2(j) == possible[i][1])]
    store2 = {k: v for k, v in store.items() if v}
    del store


#region Matching & replacing procedure
    debut = time.time()
    replacement = {}

    for j in list(store2.keys()):
        if os.path.exists(path2 + "/land4_" + str(j[0]) + "_" + str(j[1]) + ".pkl"):
            filename = open("/home/simonjean/data/all_data/land4_" + str(j[0]) + "_" + str(
                j[1]) + ".pkl", 'rb')
            set_to_check = pickle.load(filename)
            filename.close()

            for v in store2[j]:
                set_verif = utilities.equi_landscape(utilities.to_matrix(v))
                set_verif.add(tuple(v))
                list_verif = [*set_verif.intersection(set_to_check)]
            #new = list(list_verif[0])
                replacement[str(tuple(v))] = str(list_verif[0])

        elif os.path.exists(path2 + "/land4_" + str(j[0]) + "_" + str(j[1]) + "cut_0.pkl"):
            set_to_check_merg = []
            to_load = [x for x in files_data if x.startswith('land4_' + str(j[0]) + "_" + str(j[1]))]

            for x in to_load:
                filename = open("/home/simonjean/data/all_data/"+x, 'rb')
                set_to_check = pickle.load(filename)
                filename.close()
                set_to_check = list(map(tuple, set_to_check))
                set_to_check_merg.append(set_to_check)

            set_to_check_merg = [item for sublist in set_to_check_merg for item in sublist]
            set_to_check_merg = set(set_to_check_merg)

            for v in store2[j]:
                set_verif = utilities.equi_landscape(utilities.to_matrix(v))
                set_verif.add(tuple(v))
                list_verif = [*set_verif.intersection(set_to_check_merg)]
            # new = list(list_verif[0])

                replacement[str(tuple(v))] = str(list_verif[0])

            del set_to_check
        else:
            print('ERROR : no file')


    for i in list(range(len(rows))):
        rows[i] = [replacement.get(n, n) for n in rows[i]]

    tester3 = pd.DataFrame(rows, dtype=str)

    if len(nums) > 4:
        tester3.to_csv(params.path + "/keys_succession/equiv_results_" + nums[2] + "_" + nums[3] + "_cut_" + nums[4] + '.csv',
                  index=False)
    else:
        tester3.to_csv(params.path + "/keys_succession/equiv_results_" + nums[2] + "_" + nums[3] + ".csv", index=False)

    return print('Dataset ' + str(params.files[k]) +" took " + str(time.time()-debut) + ' seconds')


