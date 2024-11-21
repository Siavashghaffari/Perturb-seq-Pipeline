# A little bit of set-up
import pydsdb
import multiassayexperiment as mae
import singlecellexperiment as sce
import pandas as pd
from gpauth import GPAuth
import json
def Dataset_maker(DSID, DEV):
    """
    This function creates a test dataset for the test of the pipeline
    """
    # Retrieve the dataset
    dm = pydsdb.get_dataset(DSID, version=2, dev=DEV)
    # Get the name of our experiment
    your_experiment = list(dm.experiments.keys())[0]
    #Convert SingleCellExperiment object to AnnData object
    adata, adatas = dm.experiments[your_experiment].toAnnData(alts=True)
    # Create metadata parameters
    s = dm.metadata['.internal']["metadata"]
    json_acceptable_string = s.replace("'", "\"")
    d = json.loads(json_acceptable_string)
    title = d["dataset"]["title"]
    description = d["dataset"]["description"]
    name_space = [{"id": "GRCh38", "type": "genome"}]
    sources = [{"id": "Siavash-1234", "name": "Geo-ID"}]
    tech_name = "scRNA-seq"
    organism = 'human'
    author="SG"
    # prepare experiment level metadata¶
    sce_metadata = pydsdb.create_sce_metadata(
    description=description,
    name_space=name_space,
    organism=organism,
    sources=sources,
    technology_name=tech_name,
    title=title)
    ##### Added these lines to avoid getting a bug related to datatype:
    #####'object' dtype row values have more than one type (excluding 'None')
    if 'demux_type' in adata.obs.columns:
        adata.obs["demux_type"]=adata.obs["demux_type"].fillna("unknown")
    if 'assignment' in adata.obs.columns:
        adata.obs["assignment"]=adata.obs["assignment"].fillna("unknown")
    if 'demux_type_v2' in adata.obs.columns:
        adata.obs = adata.obs.drop('demux_type_v2', axis=1)
    if 'assignment_v2' in adata.obs.columns:
        adata.obs = adata.obs.drop('assignment_v2', axis=1)
    if 'havana_gene' in adata.var.columns:
        adata.var["havana_gene"]=adata.var["havana_gene"].fillna("unknown")
    if 'hgnc_id' in adata.var.columns:
        adata.var["hgnc_id"]=adata.var["hgnc_id"].fillna("unknown")
    if 'tag' in adata.var.columns:
        adata.var["tag"]=adata.var["tag"].fillna("unknown")
    if 'genomic_ranges' in adata.var.columns:
        adata.var = adata.var.drop('genomic_ranges', axis=1)
    # convert anndata to a SingleCellExperiment object
    tse_2 = sce.io.anndata.fromAnnData(adata)
    # attach experiment level metadata to SingleCellExperiment objects
    pydsdb.add_metadata(tse_2, sce_metadata)
    # Create a variable which is a dictionary for alternative experiments
    tse_alt=[]
    tse_alt_names =[]
    for i, (k, v) in enumerate(adatas.items()):
        tse_alt.append(sce.io.anndata.fromAnnData(v))
        tse_alt_names.append(k)
    # attach experiment level metadata to SingleCellExperiment objects of alt experiments
    for i in range(len(tse_alt)):
        pydsdb.add_metadata(tse_alt[i], sce_metadata)
    AltExps = dict(zip(tse_alt_names, tse_alt))
    # Creating singlecellexperiment object for upload preparation
    tse = sce.SingleCellExperiment(
    assays={"counts": adata.layers['counts'].T}, rowData=adata.var, colData=adata.obs,
    reducedDims={}, altExps=AltExps)
    # attach experiment level metadata to SingleCellExperiment objects
    pydsdb.add_metadata(tse, sce_metadata)
    # Convert SingleCellExperiment object to MAE object
    # Add new data to the MAE
    new_mae = mae.makeMAE({your_experiment: tse})
    # prepare dataset level metadata
    dataset_metadata = pydsdb.create_dataset_metadata(
    description=("Siavash test upload; this should be gone"),
    title="Siavash test upload",
    authors="SG")
    #add dataset metadata
    pydsdb.add_metadata(new_mae, dataset_metadata)
    # Initialize the upload
    upload = pydsdb.Upload(new_mae)
    # Set Permissions
    permissions = pydsdb.create_permissions_info()
    # authenticating interactively: Please insert your password
    auth = GPAuth(False)
    # Update the existing dataset
    dsid, version = upload.submit(auth,permissions, dev=True, test=True)
        
        
    return dsid, version