import os 
import socket
import time

def check_port(port):
    import socket
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    try:
        sock.bind(("0.0.0.0", port))
        available = True
    except:
        available = False
    sock.close()
    return available

slurm_env = [
     'export XRD_RUNFORKHANDLER=1',
     'export XRD_STREAMTIMEOUT=10',
     f'source {os.environ["HOME"]}/.bashrc',
     f'conda activate myenv',
     f'export X509_USER_PROXY={os.environ["HOME"]}/x509up_u146312'
]

extra_args=[
     "--output=DASKlogs/dask_job_output_%j.out",
     "--error=DASKlogs/dask_job_output_%j.err",
     "--partition=submit",
     "--clusters=submit",
]

from dask_jobqueue import SLURMCluster
from distributed import Client
from dask.distributed import performance_report

#cluster = dask_jobqueue.SLURMCluster(cores=12, memory='24 GB', processes=1, interface='ib0')

while not check_port(6820):
    time.sleep(5)
    
cluster = SLURMCluster(
        queue='all',
        project="Hrare_Slurm",
        cores=20, #
        processes=1,
        memory="2 GB",
        #retries=10,
        walltime='00:30:00',
        scheduler_options={
              'port': 6820,
              'dashboard_address': 8000,
              'host': socket.gethostname()
        },
        job_extra=extra_args,
        env_extra=slurm_env,
)

#cluster.scale(1)
client = Client(cluster)

print(client)
