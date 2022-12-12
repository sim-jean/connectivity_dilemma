import time
import multiprocessing as mp

import current_scripts.params as params
import current_scripts.matcher as matcher


__spec__ = None
__name__ = "__main__"

if __name__ == "__main__":
    start = time.time()
    with mp.Pool(processes=14) as pool:
        results = pool.map(matcher.match_heavy, list(range(len(params.files))))
    pool.terminate()
    print("The process took " + str((time.time()-start)//60) + ' minutes and ' + str((time.time()-start) % 60) + ' seconds')

