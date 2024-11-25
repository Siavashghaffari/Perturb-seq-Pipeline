## Import Libraries

from __future__ import division 

#single cell libraries
import scanpy as sc
import anndata

#general
import pandas as pd
import sys, os
import numpy as np
import itertools
import json
import subprocess
import papermill as pm
import warnings
#from tqdm import tqdm
warnings.simplefilter(action='ignore', category=FutureWarning) 

#DatasetDB
import pydsdb


class Features (object):
    """
    This class encapsulates all the feautres necessary for the single cell seq pipeline
    """
    def __init__(self, DatasetID, auth, **kwargs):
        
        self.DatasetID = DatasetID
        self.auth = auth
        
        # Unpack keyword arguments
        # General and metadata parameters
        self.DEV = kwargs.pop('DEV', False)
        self.test_DatasetID = kwargs.pop('test_DatasetID', '')
        self.DS_test = kwargs.pop('DS_test', False)
        self.title = kwargs.pop('title', "Dataset metadata")
        self.description = kwargs.pop('description', "Dataset decription")
        self.name_space = kwargs.pop('name_space', [{"id": "GRCh38", "type": "genome"}])
        self.sources = kwargs.pop('sources', [{"id": "1234", "name": "Geo-ID"}])
        self.tech_name = kwargs.pop('tech_name', "scRNA-seq")
        self.author = kwargs.pop('author', "SG")
        # Saving options
        self.Notebooks_HOME = kwargs.pop('Notebooks_HOME',".")
        self.Notebooks_Savepath = kwargs.pop('Notebooks_Savepath',"/gstore/home/ghaffars/Notebooks")
        self.Reports_format = kwargs.pop('Reports_format',"notebook")
        self.Dataset_HOME = kwargs.pop('Dataset_HOME',".")
        # Notebooks parameters
        self.gene_param = kwargs.pop('gene_param', dict({'organism':f'human',                                                                          'groupby': f'Sample', 'min_cells': 3, 'min_genes': 500, 'max_genes': None, 'pct_mt': None, 
        'total_counts': 500,'n_top_genes':2000, 'ribo': True, 'doublets': True}))
        self.hashing_param = kwargs.pop('hashing_param', dict({'topBarcodesToPlot':5,                                                                   'bottomBarcodesToPlot':5, 'fix_barcodes':False, 'valid_assignments':None}))        
        self.crispr_param = kwargs.pop('crispr_param', dict({'topBarcodesToPlot':5,                                                                     'bottomBarcodesToPlot':5, 'fix_barcodes':False, 'valid_assignments':None}))
        self.Embeddings_param = kwargs.pop('Embeddings_param', dict({'subsample_dic':{},                                                                    'subsample_type':'n'}))
        self.classifier_param = kwargs.pop('classifier_param', dict({'training_set':'NGS4390',                                                   'training_full':f'/gstore/project/singlecell_sm/crc_base_profiling/matrices/NGS4390.raw.qc.label.h5ad', 
            'dLevel':'leiden_label_cell', 'study_annotation':'Sample',
        'min_genes':500, 'probability_cutoff':0.7, 'margin_cutoff':1, 'use_DSDB':False}))
        self.Fisher_param = kwargs.pop('Fisher_param', dict({'genes':[],'controls':["NTC"], 
            'x_category':"Untreatedclus", 'y_category' : "DOXclus",
        'remove_gRNA_doublets': True, 'permutations' : 10000}))
        self.ED_param = kwargs.pop('ED_param', dict({'pVAL':0.05, 'n':500, 'embeddings_key':"X_pca"}))
        
        # Setting templates notebooks
        self.crispr_hashing_harmonization = os.path.join(self.Notebooks_HOME, "DatasetCrisprHashing_harmonization.ipynb")
        self.qc_notebook = os.path.join(self.Notebooks_HOME, "RNAseq_QC_datasetDB_Rapids.ipynb")
        self.barcode_notebook = os.path.join(self.Notebooks_HOME, "BarcodeCountQC_datasetDB.ipynb")
        if len(self.Embeddings_param["subsample_dic"])>0:
            self.embedding_notebook = os.path.join(self.Notebooks_HOME, "Embeddings_subsample.ipynb")
        else:
            self.embedding_notebook = os.path.join(self.Notebooks_HOME, "Embeddings.ipynb")
        self.classifier_notebook = os.path.join(self.Notebooks_HOME, "pySingleCellNetCalibrated.ipynb")
        self.fisher_notebook = os.path.join(self.Notebooks_HOME, "Fisher_test.ipynb")
        self.edist_notebook = os.path.join(self.Notebooks_HOME, "EnergyDistance.ipynb")
        self.edist_notebook_P = os.path.join(self.Notebooks_HOME, "EnergyDistance_perPerturbations.ipynb")
        self.multivariate_notebook = os.path.join(self.Notebooks_HOME, "Multivariate_updated.ipynb")
        self.mixscape_notebook = os.path.join(self.Notebooks_HOME, "Mixscape.ipynb")
        self.hotspot_notebook = os.path.join(self.Notebooks_HOME, "HotSpot_template.ipynb")
        
        # Output notebooks
        self.crispr_hashing_harmonization_out = f'{self.Notebooks_Savepath}/{self.DatasetID}_dataset_CrisprHashing_harmonization.ipynb'
        self.qc_notebook_out = f'{self.Notebooks_Savepath}/{self.DatasetID}_gex_qc.ipynb'
        self.hashing_notebook_out = f'{self.Notebooks_Savepath}/{self.DatasetID}_hashing_qc.ipynb'
        self.crispr_notebook_out = f'{self.Notebooks_Savepath}/{self.DatasetID}_crispr_qc.ipynb'
        if len(self.Embeddings_param["subsample_dic"])>0:
            self.embedding_notebook_out = f'{self.Notebooks_Savepath}/{self.DatasetID}_embeddings_{list(self.Embeddings_param["subsample_dic"].items())[0][0]}_{list(self.Embeddings_param["subsample_dic"].items())[0][1]}_{list(self.Embeddings_param["subsample_dic"].items())[1][0]}_{list(self.Embeddings_param["subsample_dic"].items())[1][1]}.ipynb'
        else:
            self.embedding_notebook_out = f'{self.Notebooks_Savepath}/{self.DatasetID}_embeddings.ipynb'
        self.classifier_notebook_out = f'{self.Notebooks_Savepath}/{self.DatasetID}_pySCN.ipynb'
        self.fisher_notebook_out = f'{self.Notebooks_Savepath}/{self.DatasetID}_FisherTest.ipynb'
        self.edist_notebook_out = f'{self.Notebooks_Savepath}/{self.DatasetID}_EnergyDistance_{self.ED_param["n"]}NTCs_{self.ED_param["embeddings_key"].split("_",1)[-1]}.ipynb'
        self.edist_notebook_out_P = f'{self.Notebooks_Savepath}/{self.DatasetID}_EnergyDistance_perPerturbations_{self.ED_param["n"]}NTCs_{self.ED_param["embeddings_key"].split("_",1)[-1]}.ipynb'
        self.multivariate_notebook_out = f'{self.Notebooks_Savepath}/{self.DatasetID}_Multivariate.ipynb'
        self.mixscape_notebook_out = f'{self.Notebooks_Savepath}/{self.DatasetID}_Mixscape.ipynb'
        self.hotspot_notebook_out = f'{self.Notebooks_Savepath}/{self.DatasetID}_hotspot.ipynb'
    
    
        # Dataset setting ouput to be stored on the disk
        self.Dataset_HOME = f'{self.Dataset_HOME}/{self.DatasetID}'
        try:
            os.mkdir(self.Dataset_HOME)
        except OSError as error:
            print(error)
        
        self.Dataset_HOME_ft = f'{self.Dataset_HOME}/Fisher-test'
        try:
            os.mkdir(self.Dataset_HOME_ft)
        except OSError as error:
            print(error)
            
        self.Dataset_HOME_classifier = f'{self.Dataset_HOME}/Classifier'
        try:
            os.mkdir(self.Dataset_HOME_classifier)
        except OSError as error:
            print(error)
        
        self.Dataset_HOME_ed = f'{self.Dataset_HOME}/Energy Distance'
        try:
            os.mkdir(self.Dataset_HOME_ed)
        except OSError as error:
            print(error)
        
        self.Dataset_HOME_mv = f'{self.Dataset_HOME}/Multivariate'
        try:
            os.mkdir(self.Dataset_HOME_mv)
        except OSError as error:
            print(error)
        
        # Setting embeddings in and out parameters
        if len(self.Embeddings_param["subsample_dic"])>0:
            embeddings_in = f'{self.Dataset_HOME}/qc_embedding.h5ad'
            embeddings_out = f'{self.Dataset_HOME}/qc_embedding_{list(self.Embeddings_param["subsample_dic"].items())[0][0]}_{list(self.Embeddings_param["subsample_dic"].items())[0][1]}_{list(self.Embeddings_param["subsample_dic"].items())[1][0]}_{list(self.Embeddings_param["subsample_dic"].items())[1][1]}.h5ad'
        else:
            embeddings_in = f'{self.Dataset_HOME}/raw_qc.h5ad'
            embeddings_out = f'{self.Dataset_HOME}/qc_embedding.h5ad'
        
        # Setting metadata parameters
        ds = pydsdb.get_dataset(self.DatasetID, version=1, dev=self.DEV)
        
        s = ds.metadata['.internal']["metadata"]
        json_acceptable_string = s.replace("'", "\"")
        #d = json.loads(json_acceptable_string)
        d = json.loads(s)
        self.your_experiment = d["dataset"]["experiments"][0]["name"]
        self.title = d["dataset"]["title"]
        self.description = d["dataset"]["description"]
        self.name_space = [{"id": "GRCh38", "type": "genome"}]
        self.tech_name = "scRNA-seq"
        # See the alternative experiments
        adata, adatas = ds.experiments[self.your_experiment].toAnnData(alts=True)
        if adatas:
            self.alternative_experiments = list(adatas.keys())
        else:
            self.alternative_experiments = None
        
        # Get the mapping rate info and store it to use later
        try:
            self.Mapping_rate = adata.uns['adt_summary']
        except KeyError:
            self.Mapping_rate = None
        
        # Setting Notebook Parameters
        self.Crisprhashing_harmonization_parameters = dict(
                DatasetID = self.DatasetID,
                #auth = self.auth,
                DEV = self.DEV,
            test_DatasetID = self.test_DatasetID,
            DS_test = self.DS_test,
           title = self.title,
        description = self.description,
        name_space = self.name_space,
        organism=self.gene_param["organism"],
        sources = self.sources,
        tech_name = self.tech_name,
        author = self.author,
        your_experiment = self.your_experiment,
        )
        
        self.QC_transcriptome_parameters = dict(
                DatasetID=self.DatasetID,
                DEV = self.DEV,
                test_DatasetID = self.test_DatasetID,
                DS_test = self.DS_test,
                title = self.title,
                description = self.description,
                name_space = self.name_space,
                sources = self.sources,
                tech_name = self.tech_name,
                author = self.author,
                alt_experiments=self.alternative_experiments, 
                organism=self.gene_param["organism"],
                groupby=self.gene_param["groupby"],
                min_cells=self.gene_param["min_cells"],
                min_genes=self.gene_param["min_genes"],
                max_genes=self.gene_param["max_genes"],
                pct_mt=self.gene_param["pct_mt"],
                total_counts=self.gene_param["total_counts"],
                n_top_genes=self.gene_param["n_top_genes"],
                ribo=self.gene_param["ribo"],
                doublets=self.gene_param["doublets"],
                 Mapping_rate = self.Mapping_rate,
                raw_qc=f'{self.Dataset_HOME}/raw_qc.h5ad',
        raw_embedding =f'{self.Dataset_HOME}/qc_embedding.h5ad')
        
        self.hashing_parameters = dict(
                DatasetID=self.DatasetID,
                DEV = self.DEV,
                test_DatasetID = self.test_DatasetID,
                DS_test = self.DS_test,
                title = self.title,
                description = self.description,
                name_space = self.name_space,
                sources = self.sources,
                tech_name = self.tech_name,
                author = self.author,
                alt_experiments=self.alternative_experiments, 
                organism=self.gene_param["organism"],
                experiment = 'hashing',
                topBarcodesToPlot = self.hashing_param["topBarcodesToPlot"],
                bottomBarcodesToPlot = self.hashing_param["bottomBarcodesToPlot"],
                fix_barcodes = self.hashing_param["fix_barcodes"],
                valid_assignments = self.hashing_param["valid_assignments"])
        
        self.crispr_parameters = dict(
                DatasetID=self.DatasetID,
                DEV = self.DEV,
                test_DatasetID = self.test_DatasetID,
                DS_test = self.DS_test,
                title = self.title,
                description = self.description,
                name_space = self.name_space,
                sources = self.sources,
                tech_name = self.tech_name,
                author = self.author,
                alt_experiments=self.alternative_experiments, 
                organism=self.gene_param["organism"],
                experiment = 'crispr',
                topBarcodesToPlot = self.crispr_param["topBarcodesToPlot"],
                bottomBarcodesToPlot = self.crispr_param["bottomBarcodesToPlot"],
                fix_barcodes = self.crispr_param["fix_barcodes"],
                valid_assignments = self.crispr_param["valid_assignments"])
        
        
        self.embeddings_parameters = dict(
                DatasetID=self.DatasetID,
                DEV = self.DEV,
                test_DatasetID = self.test_DatasetID,
                DS_test = self.DS_test,
                title = self.title,
                description = self.description,
                name_space = self.name_space,
                sources = self.sources,
                tech_name = self.tech_name,
                author = self.author,
                organism=self.gene_param["organism"],
                groupby=self.gene_param["groupby"],
                min_cells=self.gene_param["min_cells"],
                min_genes=self.gene_param["min_genes"],
                max_genes=self.gene_param["max_genes"],
                pct_mt=self.gene_param["pct_mt"],
                total_counts=self.gene_param["total_counts"],
                n_top_genes=self.gene_param["n_top_genes"],
                ribo=self.gene_param["ribo"],
                doublets=self.gene_param["doublets"],
            subsample_dic = self.Embeddings_param["subsample_dic"],
            subsample_type = self.Embeddings_param["subsample_type"],
                raw_qc=embeddings_in,
        raw_embedding =embeddings_out)
        
        self.classifier_parameters = dict(
                DatasetID=self.DatasetID,
                DEV = self.DEV,
                test_DatasetID = self.test_DatasetID,
                DS_test = self.DS_test,
                title = self.title,
                description = self.description,
                name_space = self.name_space,
                sources = self.sources,
                tech_name = self.tech_name,
                author = self.author,
                organism=self.gene_param["organism"],
                training_set = self.classifier_param["training_set"],
                training_full = self.classifier_param["training_full"],
                dLevel = self.classifier_param["dLevel"],
                study_annotation = self.classifier_param["study_annotation"],
        min_genes = self.classifier_param["min_genes"],
        probability_cutoff = self.classifier_param["probability_cutoff"],
        margin_cutoff = self.classifier_param["margin_cutoff"],
            use_DSDB = self.classifier_param["use_DSDB"],
            raw_qc=f'{self.Dataset_HOME}/raw_qc.h5ad',
            raw_embedding =f'{self.Dataset_HOME}/qc_embedding.h5ad',
        annotated_outfile = f'{self.Dataset_HOME_classifier}/raw.qc.pySCN.h5ad')
        
        self.fisher_parameters = dict(
                DatasetID=self.DatasetID,
                DEV = self.DEV,
                test_DatasetID = self.test_DatasetID,
                DS_test = self.DS_test,
                title = self.title,
                description = self.description,
                name_space = self.name_space,
                sources = self.sources,
                tech_name = self.tech_name,
                author = self.author,
                organism=self.gene_param["organism"],
                genes = self.Fisher_param["genes"],
                controls = self.Fisher_param["controls"],
                x_category = self.Fisher_param["x_category"],
                y_category = self.Fisher_param["y_category"],
        remove_gRNA_doublets = self.Fisher_param["remove_gRNA_doublets"],
        permutations = self.Fisher_param["permutations"],
            final_raw_annotated_file = f'{self.Dataset_HOME_classifier}/raw.qc.pySCN.h5ad',
        out_file = f'{self.Dataset_HOME_ft}/fisher_general.csv')
        
        self.ED_parameters = dict(
                DatasetID=self.DatasetID,
                DEV = self.DEV,
                test_DatasetID = self.test_DatasetID,
                DS_test = self.DS_test,
                title = self.title,
                description = self.description,
                name_space = self.name_space,
                sources = self.sources,
                tech_name = self.tech_name,
                author = self.author,
                organism=self.gene_param["organism"],
                pVAL = self.ED_param["pVAL"],
            n = self.ED_param["n"],
            embeddings_key = self.ED_param["embeddings_key"],
            ED_file = f'{self.Dataset_HOME}/qc_embedding.h5ad',
            out_file = f'{self.Dataset_HOME_ed}/e_dist_{self.ED_param["n"]}NTCs_{self.ED_param["embeddings_key"].split("_",1)[-1]}.h5ad',
            out_csv = f'{self.Dataset_HOME_ed}/e_dist_{self.ED_param["n"]}NTCs_{self.ED_param["embeddings_key"].split("_",1)[-1]}.csv'
        )
        
        self.ED_parameters_P = dict(
                DatasetID=self.DatasetID,
                DEV = self.DEV,
                test_DatasetID = self.test_DatasetID,
                DS_test = self.DS_test,
                title = self.title,
                description = self.description,
                name_space = self.name_space,
                sources = self.sources,
                tech_name = self.tech_name,
                author = self.author,
                organism=self.gene_param["organism"],
                pVAL = self.ED_param["pVAL"],
            n = self.ED_param["n"],
            embeddings_key = self.ED_param["embeddings_key"],
            ED_file = f'{self.Dataset_HOME}/qc_embedding.h5ad',
            out_file = f'{self.Dataset_HOME_ed}/e_dist_perPerturbations_{self.ED_param["n"]}NTCs_{self.ED_param["embeddings_key"].split("_",1)[-1]}.h5ad',
            out_csv = f'{self.Dataset_HOME_ed}/e_dist_perPerturbations_{self.ED_param["n"]}NTCs_{self.ED_param["embeddings_key"].split("_",1)[-1]}.csv'
        )
        
        
        self.MV_parameters = dict(
                DatasetID=self.DatasetID,
                DEV = self.DEV,
                test_DatasetID = self.test_DatasetID,
                DS_test = self.DS_test,
                title = self.title,
                description = self.description,
                name_space = self.name_space,
                sources = self.sources,
                tech_name = self.tech_name,
                author = self.author,
                organism=self.gene_param["organism"],
            MV_file = f'{self.Dataset_HOME}/qc_embedding.h5ad',
            out_file = f'{self.Dataset_HOME_mv}/Multivariate_complete.h5ad'
        )
    
        self.MS_parameters = dict(
                DatasetID=self.DatasetID,
                DEV = self.DEV,
                test_DatasetID = self.test_DatasetID,
                DS_test = self.DS_test,
                title = self.title,
                description = self.description,
                name_space = self.name_space,
                sources = self.sources,
                tech_name = self.tech_name,
                author = self.author,
                organism=self.gene_param["organism"],
            MS_file = f'{self.Dataset_HOME}/qc_embedding.h5ad',
            out_file = f'{self.Dataset_HOME}/Mixscape_complete.h5ad'
        )
        
        self.HS_parameters = dict(
                DatasetID=self.DatasetID,
                DEV = self.DEV,
                test_DatasetID = self.test_DatasetID,
                DS_test = self.DS_test,
                title = self.title,
                description = self.description,
                name_space = self.name_space,
                sources = self.sources,
                tech_name = self.tech_name,
                author = self.author,
                organism=self.gene_param["organism"]
        )
    
    
    
    def Dataset_Harmonization (self):
        if len(self.alternative_experiments)==2:
            pm.execute_notebook(
            self.crispr_hashing_harmonization,
            self.crispr_hashing_harmonization_out,
            parameters = self.Crisprhashing_harmonization_parameters,
            kernel_name='gpuy310'
            )
            if self.Reports_format=="html":
                command = f'jupyter nbconvert {self.crispr_hashing_harmonization_out} --to html'
                subprocess.run(command, shell=True)
                os.remove(self.crispr_hashing_harmonization_out)
            elif self.Reports_format=="pdf":
                command = f'jupyter nbconvert {self.crispr_hashing_harmonization_out} --to pdf'
                subprocess.run(command, shell=True)
                os.remove(self.crispr_hashing_harmonization_out)
            
    def QC_Transcriptome (self):
        pm.execute_notebook(
            self.qc_notebook,
            self.qc_notebook_out,
            parameters = self.QC_transcriptome_parameters,
            kernel_name='gpuy310'
        )
        if self.Reports_format=="html":
            command = f'jupyter nbconvert {self.qc_notebook_out} --to html'
            subprocess.run(command, shell=True)
            os.remove(self.qc_notebook_out)
        elif self.Reports_format=="pdf":
            command = f'jupyter nbconvert {self.qc_notebook_out} --to pdf'
            subprocess.run(command, shell=True)
            os.remove(self.qc_notebook_out)
            
    def QC_Hashing (self):
        pm.execute_notebook(
            self.barcode_notebook,
            self.hashing_notebook_out,
            parameters = self.hashing_parameters,
            kernel_name='gpuy310'
        )
        if self.Reports_format=="html":
            command = f'jupyter nbconvert {self.hashing_notebook_out} --to html'
            subprocess.run(command, shell=True)
            os.remove(self.hashing_notebook_out)
        elif self.Reports_format=="pdf":
            command = f'jupyter nbconvert {self.hashing_notebook_out} --to pdf'
            subprocess.run(command, shell=True)
            os.remove(self.hashing_notebook_out)

    def QC_Crispr (self):
        pm.execute_notebook(
            self.barcode_notebook,
            self.crispr_notebook_out,
            parameters = self.crispr_parameters,
            kernel_name='gpuy310'
        )
        if self.Reports_format=="html":
            command = f'jupyter nbconvert {self.crispr_notebook_out} --to html'
            subprocess.run(command, shell=True)
            os.remove(self.crispr_notebook_out)
        elif self.Reports_format=="pdf":
            command = f'jupyter nbconvert {self.crispr_notebook_out} --to pdf'
            subprocess.run(command, shell=True)
            os.remove(self.crispr_notebook_out)
                        
    def embeddings (self):
        pm.execute_notebook(
            self.embedding_notebook,
            self.embedding_notebook_out,
            parameters = self.embeddings_parameters,
            kernel_name='gpuy310'
        )
        if self.Reports_format=="html":
            command = f'jupyter nbconvert {self.embedding_notebook_out} --to html'
            subprocess.run(command, shell=True)
            os.remove(self.embedding_notebook_out)
        elif self.Reports_format=="pdf":
            command = f'jupyter nbconvert {self.embedding_notebook_out} --to pdf'
            subprocess.run(command, shell=True)
            os.remove(self.embedding_notebook_out)
                        
    
    
    def classifier (self):
        pm.execute_notebook(
            self.classifier_notebook,
            self.classifier_notebook_out,
            parameters = self.classifier_parameters,
            kernel_name='gpuy310'
        )
        if self.Reports_format=="html":
            command = f'jupyter nbconvert {self.classifier_notebook_out} --to html'
            subprocess.run(command, shell=True)
            os.remove(self.classifier_notebook_out)
        elif self.Reports_format=="pdf":
            command = f'jupyter nbconvert {self.classifier_notebook_out} --to pdf'
            subprocess.run(command, shell=True)
            os.remove(self.classifier_notebook_out)
                       
    def fisher_test (self):
        pm.execute_notebook(
            self.fisher_notebook,
            self.fisher_notebook_out,
            parameters = self.fisher_parameters,
            kernel_name='gpuy310'
        )
        if self.Reports_format=="html":
            command = f'jupyter nbconvert {self.fisher_notebook_out} --to html'
            subprocess.run(command, shell=True)
            os.remove(self.fisher_notebook_out)
        elif self.Reports_format=="pdf":
            command = f'jupyter nbconvert {self.fisher_notebook_out} --to pdf'
            subprocess.run(command, shell=True)
            os.remove(self.fisher_notebook_out)
                       
    def e_dist (self, perPerturbation=False):
        if perPerturbation==False:
            pm.execute_notebook(
            self.edist_notebook,
            self.edist_notebook_out,
            parameters = self.ED_parameters,
            kernel_name='gpuy310'
            )
            if self.Reports_format=="html":
                command = f'jupyter nbconvert {self.edist_notebook_out} --to html'
                subprocess.run(command, shell=True)
                os.remove(self.edist_notebook_out)
            elif self.Reports_format=="pdf":
                command = f'jupyter nbconvert {self.edist_notebook_out} --to pdf'
                subprocess.run(command, shell=True)
                os.remove(self.edist_notebook_out)
        elif perPerturbation==True:
            pm.execute_notebook(
            self.edist_notebook_P,
            self.edist_notebook_out_P,
            parameters = self.ED_parameters_P,
            kernel_name='gpuy310'
            )
            if self.Reports_format=="html":
                command = f'jupyter nbconvert {self.edist_notebook_out_P} --to html'
                subprocess.run(command, shell=True)
                os.remove(self.edist_notebook_out_P)
            elif self.Reports_format=="pdf":
                command = f'jupyter nbconvert {self.edist_notebook_out_P} --to pdf'
                subprocess.run(command, shell=True)
                os.remove(self.edist_notebook_out_P)
        
            
    def multivariate (self):
        pm.execute_notebook(
            self.multivariate_notebook,
            self.multivariate_notebook_out,
            parameters = self.MV_parameters,
            kernel_name='gpuy310'
        )
        if self.Reports_format=="html":
            command = f'jupyter nbconvert {self.multivariate_notebook_out} --to html'
            subprocess.run(command, shell=True)
            os.remove(self.multivariate_notebook_out)
        elif self.Reports_format=="pdf":
            command = f'jupyter nbconvert {self.multivariate_notebook_out} --to pdf'
            subprocess.run(command, shell=True)
            os.remove(self.multivariate_notebook_out)
            
    
    def mixscape (self):
        pm.execute_notebook(
            self.mixscape_notebook,
            self.mixscape_notebook_out,
            parameters = self.MS_parameters,
            kernel_name='gpydev39'
        )
        if self.Reports_format=="html":
            command = f'jupyter nbconvert {self.mixscape_notebook_out} --to html'
            subprocess.run(command, shell=True)
            os.remove(self.mixscape_notebook_out)
        elif self.Reports_format=="pdf":
            command = f'jupyter nbconvert {self.mixscape_notebook_out} --to pdf'
            subprocess.run(command, shell=True)
            os.remove(self.mixscape_notebook_out)
            
            
            
    def Hotspot (self):
        pm.execute_notebook(
            self.hotspot_notebook,
            self.hotspot_notebook_out,
            parameters = self.HS_parameters,
            kernel_name='gpuy310'
        )
        if self.Reports_format=="html":
            command = f'jupyter nbconvert {self.hotspot_notebook_out} --to html'
            subprocess.run(command, shell=True)
            os.remove(self.hotspot_notebook_out)
        elif self.Reports_format=="pdf":
            command = f'jupyter nbconvert {self.hotspot_notebook_out} --to pdf'
            subprocess.run(command, shell=True)
            os.remove(self.hotspot_notebook_out)






