import os 
import socket
import time

def create_remote_connection():

    slurm_env = [
        'export DASK_DISTRIBUTED__COMM__ALLOWED_TRANSPORTS=["tcp://[::]:0"]',
        'export XRD_RUNFORKHANDLER=1',
        'export XRD_STREAMTIMEOUT=10',
        f'source {os.getenv("HOME")}/.bashrc',
        f'conda activate myenv',
    ]

    extra_args=[
        "--output=DASKlogs/dask_job_output_%j.out",
        "--error=DASKlogs/dask_job_output_%j.err",
        "--partition=submit",
    ]


    from distributed import Client
    from dask_jobqueue import SLURMCluster

    cluster = SLURMCluster(
        project="Hrare_Slurm",
        job_name="test1",
        cores=1,
        memory='10GB',
        scheduler_options={
            'dashboard_address': 8000,
            'host': socket.gethostname()
        },
        silence_logs="debug",
        job_extra_directives=extra_args,
        job_script_prologue=slurm_env
    )

    cluster.scale(2)
    client = Client(cluster)

    print(client)
    print(cluster.job_script())

    return client

def create_local_connection(n_workers):

    from distributed import Client
    from dask.distributed import LocalCluster

    cluster = LocalCluster(n_workers=n_workers, threads_per_worker=1, processes=True, memory_limit="5GiB")
    try:
        client = Client(cluster,timeout='2s')
    except TimeoutError:
        pass
    return client
