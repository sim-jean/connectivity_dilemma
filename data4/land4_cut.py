import current_scripts.params as params
import current_scripts.dynamic_programming as dynprog

import os
import re
import pickle
import time
import pandas as pd
import multiprocessing as mp

__name__ = "__main__"
__spec__ = None

if __name__ == "__main__":
    path = "C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/all_data/"
    files = os.listdir(path)

    def data_list2(input):
        path = "C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/all_data/"
        files = os.listdir(path)

        filename = open(path+files[input], 'rb')
        list_nz = pickle.load(filename)
        filename.close()
        list_nz = list(list_nz)
        return list_nz


    debut2 = time.time()
    print("Process is running:")
    list_names = [["land"], ["success_biod_" + str(i) for i in params.biodiv], ["fire_biod_" + str(i) for i in params.biodiv],
                  ["value_biod_" + str(i) for i in params.biodiv]]
    flat_list_names = [item for sublist in list_names for item in sublist]

    for i in list(range(len(files))):

        size = 4
        list_nz = dynprog.data_list2(i)

        start = time.time()
        with mp.Pool(processes=10) as pool:
            results = pool.map(dynprog.lowlev_dynprog, list_nz)
        pool.terminate()
        # save results

        d = pd.DataFrame(results, columns=list(flat_list_names))

        #Opt succession
        indexes = re.findall(r'\d+', str(files[i]))
        #indexes = re.findall(r'\d+', str(files))
        if len(indexes) > 3:
            d.to_csv("C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results_habitat_availability/budget_"+str(params.budget)+"/land4_budget_"
                     + str(params.budget) + "_" + str(indexes[1]) + '_' + str(indexes[2]) + '_cut_' +
                     str(indexes[3]) + ".csv", index=False)
            print(str(time.time() - start) + " seconds elapsed for " + str(indexes[2]) + " twos among " + str(
                indexes[1]) + " ones - cut " + str(indexes[3]))
        else:
            d.to_csv("C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/results_habitat_availability/budget_" + str(
                params.budget) + "/land4_budget_" + str(params.budget) + "_" + str(indexes[1]) + '_' + str(
                indexes[2]) + ".csv", index=False)
            print(str(time.time() - start) + " seconds elapsed for " + str(indexes[2]) + " twos among " + str(
                indexes[1]))

        print(str(100 * (len(files) - i - 1) / len(files)) + "% remaining")
    print("Overall :" + str(time.time() - debut2) + "seconds have elapsed")

#################################################################################
