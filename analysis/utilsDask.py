import os 
import socket
import time

#def check_port(port):
#    import socket
#    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
#    try:
#        sock.bind(("0.0.0.0", port))
#        available = True
#    except:
#        available = False
#    sock.close()
#    return available

#while not check_port(6820):
#    time.sleep(5)

slurm_env = [
     'export XRD_RUNFORKHANDLER=1',
     'export XRD_STREAMTIMEOUT=10',
     f'source {os.getenv("HOME")}/.bashrc',
     f'conda activate myenv',
     f'export X509_USER_PROXY={os.getenv("HOME")}/x509up_u146312'
]

extra_args=[
     "--output=DASKlogs/dask_job_output_%j.out",
     "--error=DASKlogs/dask_job_output_%j.err",
#     "--partition=submit",
#     "--partition=submit-gpu1080",
#     "--clusters=submit",
]

from dask_jobqueue import SLURMCluster
from distributed import Client
from dask.distributed import performance_report

cluster = SLURMCluster(
#        queue='all',
#        project="Hrare_Slurm",
        cores=1,
        memory='2GB',
#        #retries=10,
#        walltime='00:30:00',
#        scheduler_options={
#              'port': 6820,
#              'dashboard_address': 8000,
#              'host': socket.gethostname()
#        },
    job_extra_directives=extra_args,
    job_script_prologue=slurm_env
)

cluster.adapt(maximum_jobs=30)
cluster.scale(10)
client_ = Client(cluster)

print(client_)
print(cluster.job_script())
