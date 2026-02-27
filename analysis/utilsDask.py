import os 
import socket
import time
import sys
import ROOT

from utilsAna import loadCorrectionSet, loadUserCode, loadtmvahelper

AF="MIT-daskgateway"

#MIT
if AF=="MIT-daskgateway": myDir = '/home/submit/mariadlf/Hrare/CMSSW_10_6_27_new/src/Hrare/analysis/'
sys.path.insert(0, myDir)

doINTERACTIVE=True
doLocalCluster=False

def loadUserCodeGateway(fileName):
    print('loadUserCodeGateway for ', fileName)

    from pathlib import Path
    from dask.distributed import get_worker

    try:
        #this try is needed at CERN Swan-Dask
        # Get the local directory where the file is stored

        worker = get_worker()
        localdir = worker.local_directory

        # Define the path to the file
        lib_path = Path(localdir) / fileName
        #    ROOT.gInterpreter.Declare(f'#include "{lib_path}"')

        ROOT.gInterpreter.ProcessLine(f'#include "{lib_path}"')
    except:
        pass


def init():

    if doINTERACTIVE or doLocalCluster:
        loadUserCode()
        year = 2024
        loadCorrectionSet(year)
        loadtmvahelper()
    else:
        if AF=="MIT-daskgateway":
            loadUserCodeGateway("functions.h")
            loadUserCodeGateway("sfCorrLib.h")
            loadUserCodeGateway("tmva_helper_xml.h")
            loadUserCodeGateway("tmva_helper_xgb.h")
        else:
            loadUserCode()
            year = 2024
            loadCorrectionSet(year)
            #    readDataQuality(2018)

def makeRDF(files,year):

    if doINTERACTIVE:
        RDataFrame = ROOT.RDataFrame
        df = RDataFrame("Events", files)
        init()
    else:
        RDataFrame = ROOT.RDF.Experimental.Distributed.Dask.RDataFrame
        if doLocalCluster:
            from utilsDask import create_local_connection
            NPARTITIONS=10 # assuming one launch 7 Graphs (7*15 = 105 cores)
            nCLUSTER=10 #(like scale for remote)
            connection = create_local_connection(nCLUSTER)
            print(connection)
            # The `npartitions` optional argument tells the RDataFrame how many tasks are desired

            df = RDataFrame("Events", files, daskclient=connection, npartitions=NPARTITIONS)
            ROOT.RDF.Experimental.Distributed.initialize(init)  # this load the .h files

            df._headnode.backend.distribute_unique_paths(
                [
                    myDir+"./config/functions.h",
                    myDir+"./config/sfCorrLib.h",
                    myDir+"./config/tmva_helper_xml.h",
                    myDir+"./config/tmva_helper_xgb.h",
                    myDir+"helper_tmva.py",
                    myDir+"utilsAna.py",
                ]
            )
            print('DONE with LocalCluster')
        else:
            # this is set up for AF=="MIT-daskgateway"
            from utilsDask import create_DaskGateway_MIT
            client = create_DaskGateway_MIT()
            client.upload_file(myDir+"utilsAna.py")

            print(client)
            NPARTITIONS=5 # assuming one launch 10 Graphs (10*30 = 200 cores)

            #Create the RDF and initialize user lib

            df = RDataFrame("Events", files, daskclient=client, npartitions=NPARTITIONS)
            df._headnode.backend.distribute_unique_paths(
                [
                    myDir+"./config/functions.h",
                    myDir+"./config/sfCorrLib.h",
                    myDir+"./config/tmva_helper_xml.h",
                    myDir+"./config/tmva_helper_xgb.h",
                    #                        "../weights_mva_oct/ggH_rho/TMVAClassification_wp80_6vars_mh110-160_BDTG.weights.xml"
                ]
            )

            ROOT.RDF.Experimental.Distributed.initialize(init)

    return df

def close_DaskGateway_MIT(gateway):

    clusters = gateway.list_clusters()
    print(clusters)
    print(' -- inside close_DaskGateway_MIT')
    for cl in clusters:
        cluster_name = cl.name
        print("old cluster_name",cluster_name)
        cluster = gateway.connect(cluster_name)
        cluster.shutdown()

def create_DaskGateway_MIT():

    print('HELLO -- inside create_DaskGateway_MIT')

    from dask_gateway import Gateway, GatewayCluster, BasicAuth

    gateway = Gateway(address="http://submit.mit.edu:6820",
                      proxy_address="http://submit.mit.edu:6821")

    close_DaskGateway_MIT(gateway)

    options = gateway.cluster_options()
    options['environment'] = "/work/submit/mariadlf/miniforge3/envs/myenv"
    options['worker_memory'] = 8.0
    options['worker_cores'] = 4

    cluster = gateway.new_cluster(options)
    cluster.scale(64)

    # need to close all the old clusters first
    clusters = gateway.list_clusters()
    for cl in clusters:
        cluster_name = cl.name
        print("cluster_name",cluster_name)

    client = cluster.get_client()
#    print('HELLO -- total numer of cores (64x4 ? ) ',client.get_ncores())
    # Get the number of workers

    cluster_info = client.scheduler_info()
    num_workers = len(cluster_info['workers'])

    # Calculate the total number of cores
    total_cores = num_workers * options['worker_cores']
    print('num_workers=',num_workers,' total_cores=',total_cores)

    print(client)

    return client

#    print(cluster.scheduler_info)
#    cluster = gateway.connect(cluster_name)
#    cluster.shutdown()

def create_remote_Dask():

    print('setting up Dask + SLURM')
    slurm_env = [
        'export DASK_DISTRIBUTED__COMM__ALLOWED_TRANSPORTS=["tcp://[::]:0"]',
        'export XRD_RUNFORKHANDLER=1',
        'export XRD_STREAMTIMEOUT=10',
        'echo "Landed on $HOSTNAME"',
        #'export DASK_CONFIG=dask/dask.yaml',
        #~/.config/dask/dask.yaml
        f'source {os.getenv("HOME")}/.bashrc',
        f'conda activate myenvAF',
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
