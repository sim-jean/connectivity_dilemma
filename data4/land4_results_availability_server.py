exec(open("/data/public_data/connectivity_jean/functions_project_server.py", encoding="utf-8", errors='ignore').read(), globals())

import current_scripts.params as params
import current_scripts.matcher as matcher
import os
import re
import time
import multiprocessing as mp
import pandas as pd

__name__ = "__main__"
__spec__ = None

if __name__ == "__main__":
    path = "/data/public_data/connectivity_jean/all_data/"
    files = os.listdir(path)

    #final_files = set(files) - set(files2)
    #final_files = list(final_files)



    debut2 = time.time()
    print("Process is running")
    success_biod = ["success_biod_" + str(i) for i in params.biodiv]
    fire_biod = ["fire_biod_" + str(i) for i in params.biodiv]
    value_biod = ["value_biod_" + str(i) for i in params.biodiv]
    list_names = [["land"], success_biod, fire_biod, value_biod]

    flat_list_names = [item for sublist in list_names for item in sublist]
    for i in list(range(len(files))):
        # proceed to normal routine
        # and list(range(
        size = 4
        list_nz = matcher.data_list2(i)

        start = time.time()

        with mp.Pool(processes=16) as pool:
            results = pool.starmap(matcher.lowlev_dynprog_cut_mod, list(zip(list_nz, [i] * len(list_nz))))
        pool.terminate()

        d = pd.DataFrame(results, columns=list(flat_list_names))

        # Save output
        indexes = re.findall(r'\d+', str(files[i]))
        if len(indexes) > 3:
            d.to_csv("/data/public_data/connectivity_jean/results_habitat_availability/budget_" + str(
                params.budget) + "/land4_budget_" + str(params.budget) + "_" + str(indexes[1]) + '_' + str(
                indexes[2]) + '_cut_' + str(indexes[3]) + ".csv", index=False)
            print(str(time.time() - start) + " elapsed for " + str(indexes[2]) + " twos among " + str(
                indexes[1]) + " ones - cut " + str(indexes[3]))
        else:
            d.to_csv("/data/public_data/connectivity_jean/results_habitat_availability/budget_" + str(params.budget) + "/land4_budget_" + str(params.budget) + "_" + str(indexes[1]) + '_' + str(
                indexes[2]) + ".csv", index=False)
            print(str(time.time() - start) + " elapsed for " + str(indexes[2]) + " twos among " + str(
                indexes[1]))

        print(str(100 * (len(files) - i - 1) / len(files)) + "% remaining")
    print("Overall :" + str(time.time() - debut2) + "seconds have elapsed")