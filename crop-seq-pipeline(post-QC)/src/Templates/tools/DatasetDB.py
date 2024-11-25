## This script provides a convient way to work with pydsdb
## @author  Siavash Ghaffari

import pydsdb
import multiassayexperiment as mae
import singlecellexperiment as sce
import pandas as pd
from gpauth import GPAuth
#from multiassayexperiment import makeMAE

# to ignore SSL errors
from os import environ
environ["GP_DISABLE_SSL_VERIFICATION"] = "True"

class DATASET (object):
    """
    This class encapsulates all the methods and instances necessary for the DATASET download, upload and updtae
    """
    def __init__(self, DSID, DEV, **kwargs):
        self.DSID = DSID
        #self.Version = Version
        self.DEV = DEV
        
        # Unpack keyword arguments
        self.title = kwargs.pop('title', "Dataset metadata")
        self.description = kwargs.pop('description', "Dataset decription")
        self.name_space = kwargs.pop('name_space', [{"id": "GRCh38", "type": "genome"}])
        self.organism = kwargs.pop('organism', f'human')
        self.sources = kwargs.pop('sources', [{"id": "1234", "name": "Geo-ID"}])
        self.tech_name = kwargs.pop('tech_name', "scRNA-seq")
        self.author = kwargs.pop('author', "SG")
        
        self.sce_metadata = pydsdb.create_sce_metadata(
        description=self.description,
        name_space=self.name_space,
        organism=self.organism,
        sources=self.sources,
        technology_name=self.tech_name,
        title=self.title)

    def load_dataset(self, Version=1, your_experiment=None, Corr=True):
        dm = pydsdb.get_dataset(self.DSID, version=Version, dev=self.DEV)
        
        if your_experiment == None:
            your_experiment = list(dm.experiments.keys())[0]
        adata, adatas = dm.experiments[your_experiment].toAnnData(alts=True)
        # In a new version of AnnData adata.X is None, so compatible to this version 
        # we define it with layers['counts'] and and pass it to adata.X to be compatible with our scripts
        adata.X = adata.layers['counts']
        ##### Added these lines to avoid getting a bug related to datatype:
        #####'object' dtype row values have more than one type (excluding 'None')
        if Corr:
            if 'demux_type' in adata.obs.columns:
                adata.obs["demux_type"]=adata.obs["demux_type"].fillna("unknown")
            if 'assignment' in adata.obs.columns:
                adata.obs["assignment"]=adata.obs["assignment"].fillna("unknown")
            if 'cellline' in adata.obs.columns:
                adata.obs["cellline"]=adata.obs["cellline"].fillna("unknown")
            if 'timepoint' in adata.obs.columns:
                adata.obs["timepoint"]=adata.obs['timepoint'].fillna("unknown")
            if 'gene_symbol' in adata.obs.columns:
                adata.obs["gene_symbol"]=adata.obs["gene_symbol"].fillna("unknown")
            if 'class' in adata.obs.columns:
                adata.obs["class"]=adata.obs["class"].fillna("unknown")
            if 'HTO' in adata.obs.columns:
                adata.obs["HTO"]=adata.obs["HTO"].fillna("unknown")
            if 'NGS_ID' in adata.obs.columns:
                adata.obs["NGS_ID"]=adata.obs["NGS_ID"].fillna("unknown")
            if 'Biological_replicate' in adata.obs.columns:
                adata.obs["Biological_replicate"]=adata.obs["Biological_replicate"].fillna("unknown")
            if 'sublibrary' in adata.obs.columns:
                adata.obs["sublibrary"]=adata.obs["sublibrary"].fillna("unknown")
            if 'gRNA_library_MOI' in adata.obs.columns:
                adata.obs["gRNA_library_MOI"]=adata.obs["gRNA_library_MOI"].fillna("unknown")
            #if 'demux_type_v2' in adata.obs.columns:
            #    adata.obs = adata.obs.drop('demux_type_v2', axis=1)
            #if 'assignment_v2' in adata.obs.columns:
            #    adata.obs = adata.obs.drop('assignment_v2', axis=1)
            if 'demux_type_v2' in adata.obs.columns:
                adata.obs["demux_type_v2"]=adata.obs["demux_type_v2"].fillna("unknown")
            if 'assignment_v2' in adata.obs.columns:
                adata.obs["assignment_v2"]=adata.obs["assignment_v2"].fillna("unknown")
            if 'havana_gene' in adata.var.columns:
                adata.var["havana_gene"]=adata.var["havana_gene"].fillna("unknown")
            if 'hgnc_id' in adata.var.columns:
                adata.var["hgnc_id"]=adata.var["hgnc_id"].fillna("unknown")
            if 'tag' in adata.var.columns:
                adata.var["tag"]=adata.var["tag"].fillna("unknown")
            if 'genomic_ranges' in adata.var.columns:
                adata.var = adata.var.drop('genomic_ranges', axis=1)
        
        self.Coldata = dm.colData
        self.Samplemap = dm.sampleMap
        self.experiments = dm.experiments
        self.metadata = dm.metadata
        
        return (adata,adatas)
    
    def get_metadata_param (self):
        s = ds.metadata['.internal']["metadata"]
        json_acceptable_string = s.replace("'", "\"")
        d = json.loads(json_acceptable_string)
        title = d["dataset"]["title"]
        description = d["dataset"]["description"]
        name_space = [{"id": "GRCh38", "type": "genome"}]
        sources = [{"id": "Siavash-1234", "name": "Geo-ID"}]
        tech_name = "scRNA-seq"
        organism = 'human'
        author = "SG"
        return title, description, name_space, sources,tech_name, organism, author 

    def upload_dataset(self, new_experiment, adata_updated, adatas_updated, test):
        """
        This function will upload the updated dataset (new experimnet added) with a new dataset ID
        """
        # tse_2 is fisrt scenario in which there is no alternative experiments
        # convert anndata to a SingleCellExperiment object
        tse_2 = sce.io.anndata.fromAnnData(adata_updated)
        # attach experiment level metadata to SingleCellExperiment objects
        pydsdb.add_metadata(tse_2, self.sce_metadata)

        # Create a variable which is a dictionary for alternative experiments
        tse_alt=[]
        tse_alt_names =[]
        for i, (k, v) in enumerate(adatas_updated.items()):
            tse_alt.append(sce.io.anndata.fromAnnData(v))
            tse_alt_names.append(k)

        # attach experiment level metadata to SingleCellExperiment objects of alt experiments
        for i in range(len(tse_alt)):
            pydsdb.add_metadata(tse_alt[i], self.sce_metadata)

        AltExps = dict(zip(tse_alt_names, tse_alt))  

        # Creating singlecellexperiment object for upload preparation
        tse = sce.SingleCellExperiment(
        assays={"counts": adata_updated.layers['counts'].T}, rowData=adata_updated.var, colData=adata_updated.obs,
        reducedDims={}, altExps=AltExps
        )

        # attach experiment level metadata to SingleCellExperiment objects
        pydsdb.add_metadata(tse, self.sce_metadata)

        # Convert SingleCellExperiment object to MAE object
        # Add new data to the MAE
        maeobj = mae.makeMAE({new_experiment: tse})
        # these two are needed if new experiment(s) added
        column_data = pd.concat([maeobj.colData, self.Coldata ])
        sample_map = pd.concat([maeobj.sampleMap, self.Samplemap])
        
        # Create new MAE object to upload
        new_mae = mae.MultiAssayExperiment(
        {**self.experiments,New_experiment: tse},
        column_data,
        sample_map,
        self.metadata
        )

        # Initialize the upload
        upload = pydsdb.Upload(new_mae)
        # Set Permissions
        permissions = pydsdb.create_permissions_info()
        # authenticating interactively: Please insert your password
        auth = GPAuth(False)
        # Update the existing dataset
        dsid, version = upload.submit(auth,permissions, dev=self.DEV, test=test)
        return dsid, version

    def update_dataset(self, new_experiment, adata_updated, adatas_updated, test, Corr=True):
        """
        This function will updtae the updated dataset (new experimnet added) with a new dataset Id and a new version
        """
        # tse_2 is fisrt scenario in which there is no alternative experiments
        # convert anndata to a SingleCellExperiment object
        if Corr:
            dm = pydsdb.get_dataset(self.DSID, version=1, dev=self.DEV)
        tse_2 = sce.io.anndata.fromAnnData(adata_updated)
        # attach experiment level metadata to SingleCellExperiment objects
        pydsdb.add_metadata(tse_2, self.sce_metadata)

        # Create a variable which is a dictionary for alternative experiments
        tse_alt=[]
        tse_alt_names =[]
        for i, (k, v) in enumerate(adatas_updated.items()):
            tse_alt.append(sce.io.anndata.fromAnnData(v))
            tse_alt_names.append(k)

        # attach experiment level metadata to SingleCellExperiment objects of alt experiments
        for i in range(len(tse_alt)):
            pydsdb.add_metadata(tse_alt[i], self.sce_metadata)

        AltExps = dict(zip(tse_alt_names, tse_alt))  

        # Creating singlecellexperiment object for upload preparation
        tse = sce.SingleCellExperiment(
        assays={"counts": adata_updated.layers['counts'].T}, rowData=adata_updated.var, colData=adata_updated.obs,
        reducedDims={}, altExps=AltExps
        )

        # attach experiment level metadata to SingleCellExperiment objects
        pydsdb.add_metadata(tse, self.sce_metadata)

        # Convert SingleCellExperiment object to MAE object
        # Add new data to the MAE
        maeobj = mae.makeMAE({new_experiment: tse})
        # these two are needed if new experiment(s) added
        if Corr:
            column_data = pd.concat([maeobj.colData, dm.colData])
            sample_map = pd.concat([maeobj.sampleMap, dm.sampleMap])
        else:
            column_data = pd.concat([maeobj.colData, self.Coldata ])
            sample_map = pd.concat([maeobj.sampleMap, self.Samplemap])
        
        # avoid geeting datatype error
        if 'samples' in column_data:
            column_data = column_data.drop('samples', axis=1)
        if 'has.crispr' in column_data:
            column_data["has.crispr"]=column_data['has.crispr'].fillna(True)
        if 'has.hashing' in column_data:
            column_data["has.hashing"]=column_data['has.hashing'].fillna(True)            
            
        # Create new MAE object to upload
        if Corr:
            new_mae = mae.MultiAssayExperiment(
            {**dm.experiments,new_experiment: tse},
            column_data,
            sample_map,
            dm.metadata
            )
        else:
            new_mae = mae.MultiAssayExperiment(
            {**self.experiments,new_experiment: tse},
            column_data,
            sample_map,
            self.metadata
            )
        # Initialize the upload
        upload = pydsdb.Upload(new_mae)
        # Set Permissions
        permissions = pydsdb.create_permissions_info()
        # authenticating interactively: Please insert your password
        auth = GPAuth(False)
        # Update the existing dataset
        dsid, version = upload.submit(auth,permissions, self.DSID, dev=self.DEV, test=test)
        return dsid, version
    
    def update_another_dataset(self, test_DSID, test_dev, test_version, new_experiment, adata_updated, adatas_updated, test):
        """
        This function will updtae another dataset (new experimnet added) which is uusually a test dataset for testing purposes
        """
        dm = pydsdb.get_dataset(test_DSID, version=test_version, dev=test_dev)
        # tse_2 is fisrt scenario in which there is no alternative experiments
        # convert anndata to a SingleCellExperiment object
        tse_2 = sce.io.anndata.fromAnnData(adata_updated)
        # attach experiment level metadata to SingleCellExperiment objects
        pydsdb.add_metadata(tse_2, self.sce_metadata)

        # Create a variable which is a dictionary for alternative experiments
        tse_alt=[]
        tse_alt_names =[]
        for i, (k, v) in enumerate(adatas_updated.items()):
            tse_alt.append(sce.io.anndata.fromAnnData(v))
            tse_alt_names.append(k)

        # attach experiment level metadata to SingleCellExperiment objects of alt experiments
        for i in range(len(tse_alt)):
            pydsdb.add_metadata(tse_alt[i], self.sce_metadata)

        AltExps = dict(zip(tse_alt_names, tse_alt))  

        # Creating singlecellexperiment object for upload preparation
        tse = sce.SingleCellExperiment(
        assays={"counts": adata_updated.layers['counts'].T}, rowData=adata_updated.var, colData=adata_updated.obs,
        reducedDims={}, altExps=AltExps
        )

        # attach experiment level metadata to SingleCellExperiment objects
        pydsdb.add_metadata(tse, self.sce_metadata)

        # Convert SingleCellExperiment object to MAE object
        # Add new data to the MAE
        maeobj = mae.makeMAE({new_experiment: tse})
        # these two are needed if new experiment(s) added
        column_data = pd.concat([maeobj.colData, dm.colData ])
        sample_map = pd.concat([maeobj.sampleMap, dm.sampleMap])
        
        # Create new MAE object to upload
        new_mae = mae.MultiAssayExperiment(
        {**dm.experiments,new_experiment: tse},
        column_data,
        sample_map,
        dm.metadata
        )

        # Initialize the upload
        upload = pydsdb.Upload(new_mae)
        # Set Permissions
        permissions = pydsdb.create_permissions_info()
        # authenticating interactively: Please insert your password
        auth = GPAuth(False)
        # Update the existing dataset
        dsid, version = upload.submit(auth,permissions, test_DSID, dev=test_dev, test=test)
        return dsid, version
    
    
    
def upload(adata, new_experiment, dev, test, title, description, name_space, organism, sources, tech_name, author):
    """
    This function just upload a dataset which does niot exist before into datasetDB and get it registered
    with a new ID and version 1. This case is useful once we have no DSDBID to update
    """
    
    #Create metadata
    scr_metadata = pydsdb.create_scr_metadata(
    description=description,
    name_space=name_space,
    organism=organism,
    sources=sources,
    technology_name=tech_name,
    title=title,
    # other optional metadata
    )
    
    # Creating singlecellexperiment object for upload preparation
    sExpt = sce.fromAnnData(adata)
    pydsdb.add_metadata(sExpt, scr_metadata)
    
    # prepare dataset level metadata
    dataset_metadata = pydsdb.create_dataset_metadata(
    description=("SIA test upload; this should be gone"),
    title="ghaffars test upload",
    authors="ghaffars",
    )
    
    # Convert SingleCellExperiment object to MAE object
    obj = mae.makeMAE({new_experiment: sExpt})
    
    # attach dataset level metadata
    pydsdb.add_metadata(obj, dataset_metadata)
    
    # Initialize the upload
    upload = pydsdb.Upload(obj)
    
    # Set Permissions
    permissions = pydsdb.create_permissions_info()
    
    # authenticating interactively: Please insert your password
    auth = GPAuth(False)
    
    # Update the existing dataset
    dsid, version = upload.submit(auth,permissions, dev=DEV, test=False, mode="sts:boto3")
    
    return dsid, version