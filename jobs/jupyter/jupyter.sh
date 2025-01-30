#!/bin/bash -l

#SBATCH --job-name=ops-jupyter
#SBATCH --partition=nd_bioinformatics_cpu,cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=1-00:00:00
#SBATCH --output=/scratch/users/%u/slurm_jobs/%j_%x.log
#SBATCH --chdir /scratch/prj/ppn_als_longread/ont-benchmark
#SBATCH --signal=USR2

module load anaconda3
source activate jupyter

# Get unused remote socket
readonly IPADDRESS=$(hostname -I | tr ' ' '\n' | grep '10.211.4.')
readonly REMOTE_PORT=$(python -c 'import socket; s=socket.socket(); s.bind(("", 0)); print(s.getsockname()[1]); s.close()')

# Find an available local port
LOCAL_PORT=$(python -c '
import socket
port = 8888
while True:
    try:
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.bind(("localhost", port))
        s.close()
        print(port)
        break
    except socket.error:
        port += 1
')

cat 1>&2 <<END
1. Run the following SSH command on your local machine to create the tunnel:

   ssh -NL ${LOCAL_PORT}:${HOSTNAME}:${REMOTE_PORT} create

2. Point your web browser to http://localhost:${LOCAL_PORT}/lab?token=<add the token from the jupyter output below>

When done using the notebook, terminate the job by
issuing the following command on the login node:

      scancel -f ${SLURM_JOB_ID}

END

jupyter-lab --port=${REMOTE_PORT} --ip=${IPADDRESS} --no-browser

printf 'notebook exited' 1>&2
