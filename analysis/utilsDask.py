import os 
import socket
import time

def create_Purdue_connection():

    print("inside create_remote_connection")
    import dask_gateway
    from dask_gateway import Gateway
    from distributed import Client

    os.environ["X509_USER_PROXY"] = "/home/dalfonso-cern/JULY_ex/x509up_u146312"
    print("done with environ")

    gateway = Gateway(
        "http://dask-gateway-k8s.geddes.rcac.purdue.edu/",
        proxy_address="traefik-dask-gateway-k8s.cms.geddes.rcac.purdue.edu:8786",
    )
    print(gateway)

    clusters = gateway.list_clusters()
    cluster_name = clusters[0].name
    print("cluster_name",cluster_name)

    options = gateway.cluster_options()
    print("options.worker_cores",options.worker_cores)
    print("options.worker_memory",options.worker_memory)
    print("options.environment",options.environment)

    # Create the cluster                                                                                                                                                                                
    cluster = gateway.new_cluster(
        conda_env = "/depot/cms/kernels/root632", # path to conda env
        worker_cores = 1,    # cores per worker
        worker_memory = 10,   # memory per worker in GB
        env = dict(os.environ), # pass environment as a dictionary
    )
    cluster.scale(10)
    print(cluster)

    client = Client(cluster)
    print(client)
    return client

def set_env():
    os.environ["X509_USER_PROXY"] = "/tmp/x509up_u7592"

def create_remote_DaskGateway():

    from dask_gateway import Gateway, GatewayCluster, BasicAuth

    gateway = Gateway(address="http://submit.mit.edu:6820",
                      proxy_address="http://submit.mit.edu:6821")

    clusters = gateway.list_clusters()
    for cl in clusters:
        cluster_name = cl.name
        print("cluster_name",cluster_name)

    # make a new cluster
    cluster = gateway.new_cluster(
#        env = "/work/submit/mariadlf/miniforge3/envs/myenvAF", # path to conda env
        worker_cores = 1,    # cores per worker
        worker_memory = 6,   # memory per worker in GB
    )

    cluster.scale(10)

    options = gateway.cluster_options()
    print("options.worker_cores",options.worker_cores)
    print("options.worker_memory",options.worker_memory)

    # get a client
    client = cluster.get_client()
    client.run(set_env)

    return client


def create_remote_Dask():

    slurm_env = [
        'export DASK_DISTRIBUTED__COMM__ALLOWED_TRANSPORTS=["tcp://[::]:0"]',
        'export XRD_RUNFORKHANDLER=1',
        'export XRD_STREAMTIMEOUT=10',
        'echo "Landed on $HOSTNAME"',
        f'source {os.getenv("HOME")}/.bashrc',
        f'conda activate myenvLuca',
    ]

    extra_args=[
        "--output=DASKlogs/dask_job_output_%j.out",
        "--error=DASKlogs/dask_job_output_%j.err",
        "--partition=submit,submit-gpu,submit-gpu-a30",
    ]


    from distributed import Client
    from dask_jobqueue import SLURMCluster
    import warnings

    cluster = SLURMCluster(
        project="Hrare_Slurm",
        job_name="test1",
        cores=1,
        memory='10GB',
        walltime='00:30:00',
        scheduler_options={
            'dashboard_address': 8000,
            'host': socket.gethostname()
        },
        silence_logs="debug",
        job_extra_directives=extra_args,
        job_script_prologue=slurm_env
    )

    cluster.scale(10)
    client = Client(cluster)

    return client

def create_local_connection(n_workers):

    from distributed import Client
    from dask.distributed import LocalCluster

    cluster = LocalCluster(n_workers=n_workers, threads_per_worker=1, processes=True, memory_limit="10GiB")
    try:
        client = Client(cluster,timeout='60s') # 10minutes
    except TimeoutError:
        pass
    return client
