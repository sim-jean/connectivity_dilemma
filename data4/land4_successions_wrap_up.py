import time
import multiprocessing as mp

import current_scripts.params as params
import current_scripts.successions as successions
import current_scripts.outputs as outputs

#find the maximum of biodiversity value
maxbiod = outputs.Land(params.m).biod

step = 10
__name__ = '__main__'
__spec__ = None

if __name__ == "__main__":
    start = time.time()
    for j in params.files:
        candidates = list(zip([j]*len(list(range(2, maxbiod + 1, step))), list(range(2, maxbiod + 1, step))))
        with mp.Pool(processes=12) as pool:
            results = pool.starmap(successions.succession_dataframe, candidates)
        pool.terminate()
    print("The process took " + str((time.time()-start)//60) + ' minutes and ' + str((time.time()-start) % 60) + ' seconds')