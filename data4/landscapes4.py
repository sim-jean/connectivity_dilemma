# Program for the landscape connectivity dilemma - 4x4 landscape, Mouysset, 2022
import sys

globals().clear()

#region Import packages and functions
import numpy as np
import time
import multiprocessing as mp
import os
import pickle
import itertools
exec(open('current_scripts/params.py').read())
import current_scripts.params as params
import current_scripts.utilities as utilities
#endregion


#region Get unique landscapes as dictionaries and save them
__spec__ = None
__name__ = "__main__"


debut = time.time()
if __name__ == "__main__":
    size = 16
    list_non_z = list(range(size + 1))
    list_non_z.remove(0)
    store = dict()
    #Loop for number of non-zero elements in the landscape
    for i in list_non_z:
        all_two = list(range(i + 1))
        # Number of possible 2 (including 0)
        non_zz = [i]*(i + 1)
        # Number of non_zero elements
        start_time = time.time()

        with mp.Pool(processes=16) as pool:
            results = pool.starmap(utilities.low_lev_land4, list(zip(all_two, non_zz)))
        pool.terminate()
        print("Step :", i)
        print("took %s seconds" % (time.time() - start_time))

        #Save file
        for k in range(len(results)):
            a_file = open("C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/land4_" + str(i) + "_" + str(k) + ".pkl", "wb")
            pickle.dump(results[k], a_file)
            a_file.close()
print("Overall time : %s", (time.time() - debut))


results = list(itertools.product(range(1), repeat=16))


a_file = open("C:/Users/jean/PycharmProjects/connectivity_dilemma/data3/land4_0.pkl", "wb")
pickle.dump(results, a_file)
a_file.close()
#endregion

#region Cut large files
admissible = list(range(17))
admissible.remove(0)
list_cases = [list(itertools.product(list(range(i + 1)), [i])) for i in admissible]
list_cases = [item for sublist in list_cases for item in sublist]

# set size of files
path = "C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/to_cut/land4_"
min_size = os.path.getsize(path+str(7)+"_"+str(1)+".pkl")
files = os.listdir("C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/to_cut")

done_files = True


if done_files == True:
    for k in range(len(list_cases)):
        print(k/len(list_cases))
        a = list_cases[k][1]
        b = list_cases[k][0]
        if os.path.getsize(path + str(a) + "_" + str(b) + ".pkl") >= min_size:
            filename = open(path + str(a) + "_" + str(b) + ".pkl", 'rb')
            list_nz = pickle.load(filename)
            filename.close()
            list_nz = list(list_nz)

            n = round(os.path.getsize(path + str(a) + "_" + str(b) + ".pkl")/min_size)
            new_l = np.array_split(list_nz, n)

            for l in range(len(new_l)):
                new_path =  "C:/Users/jean/PycharmProjects/connectivity_dilemma/data4/cut/land4_"
                filename = open(new_path + str(a) + "_" + str(b) + "cut_"+str(l)+".pkl", 'wb')
                pickle.dump(new_l[l], filename)
                filename.close()
#endregion

filename = open(path + str(6) + "_" + str(3) + ".pkl", 'rb')
list_nz = pickle.load(filename)
filename.close()
list_nz = list(list_nz)

filename = open(path + str(6) + "_" + str(3) + "cut_0.pkl", 'rb')
list_nz2 = pickle.load(filename)
filename.close()
list_nz2 = list(list_nz2)