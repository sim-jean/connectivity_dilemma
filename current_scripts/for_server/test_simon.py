import modules.params_server as params
import os
import multiprocessing as mp

print('Hello World!')

print('Successful import of parameters : budget is '
      + str(params.budget))
print('Size is : ' + str(params.size))



n_cores = os.getenv('SLURM_NTASKS', '1')  # env var is always a 'str'
n_cores = int(n_cores)                    # coerce to 'int'
print('Number of core is requested: ',  n_cores, '\n')