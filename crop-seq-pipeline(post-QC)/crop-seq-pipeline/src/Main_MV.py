#! /apps/user/gpy/envs/dev/GPUy310pascal/bin/python

#SBATCH -p himem
#SBATCH --mem=256G
#SBATCH --qos=long
#SBATCH -c 8

from Templates.tools.features import *
import os

auth = None

Notebooks_HOME = "./Templates/"   # the directory that contains input notebooks
A = os.getcwd()
Notebooks_Savepath = os.path.abspath(os.path.join(A ,"../Reports"))  # the directory that contains exported reports
Reports_format = "notebook" # the output report format can be one of "notebook", "html" and "pdf" 
Dataset_HOME = "/gstore/project/crc_recursion_gw/DLD1_Sublib4/" # the directory that contains all AnnData files stored on the disk

DatasetID = 'DS000016763'  # the dataset Id based on what you received when you submitted your dataset
DEV = False    # A boolean based on whether your dataset saved to the dev server

DS_test = False  #A boolean based on whether your want to create a test dataset to run the pipeline on 

# any identifiers you want to attach to the experiment, usually GEO ids, publication etc.
sources = [{"id": "Siavash-1234", "name": "Geo-ID"}]
#Authors names
author = "SG"

gene_param = {}
hashing_param = {}
crispr_param = {}
classifier_param = {}
Fisher_param = {}
ED_param ={}
MV_param={}

organism = 'human' # human or mouse
groupby = 'Sample' # For tabulation
min_cells = 3
min_genes = 500
max_genes = None
pct_mt = None
total_counts = None
n_top_genes = 2000
ribo = True
doublets = True # This refers to hashing. Only for hashing when you have demux_type in obs
# create a dictionary to keep the parameters organized
gene_param = dict({'organism':organism, 'groupby':groupby, 'min_cells':min_cells, 'min_genes':min_genes,
                  'max_genes':max_genes, 'pct_mt':pct_mt, 'total_counts':total_counts, 'n_top_genes':n_top_genes,
                  'ribo':ribo, 'doublets':doublets})


#experiment = 'hashing'
topBarcodesToPlot = 5
bottomBarcodesToPlot = 5
fix_barcodes=False
valid_assignments=None
# create a dictionary to keep the parameters organized
hashing_param = dict({'topBarcodesToPlot':topBarcodesToPlot, 'bottomBarcodesToPlot':bottomBarcodesToPlot, 
                                       'fix_barcodes':fix_barcodes, 'valid_assignments':valid_assignments})


#experiment = 'crispr'
topBarcodesToPlot = 5
bottomBarcodesToPlot = 5
fix_barcodes=False
valid_assignments=None
# create a dictionary to keep the parameters organized
crispr_param = dict({'topBarcodesToPlot':topBarcodesToPlot, 'bottomBarcodesToPlot':bottomBarcodesToPlot, 
                                       'fix_barcodes':fix_barcodes, 'valid_assignments':valid_assignments})

### Define Norm_control Parameters
use_DSDB = False   # A boolean based on whether use dataset from DSDB or sstored on the disk
# create a dictionary to keep the parameters organized
Norm_param = dict({'use_DSDB':use_DSDB})

#subsample_dic = {'CTNNB1':200, 'TCF7L2':200, 'CD81':200, 'MYC':200} # a dictionary with genes of interest as keys and subsample size
subsample_dic = {}
subsample_type='frac' # 'n': number of random sample of items  or 'frac': fraction of items to return
# create a dictionary to keep the parameters organized
Embeddings_param = dict({'subsample_dic':subsample_dic, 'subsample_type':subsample_type})


training_set = 'NGS4390'
training_full = f'/gstore/project/singlecell_sm/crc_base_profiling/matrices/{training_set}.raw.qc.label.h5ad'
dLevel = 'leiden_label_cell'
study_annotation = 'DemuxAssignment_hashing'
min_genes=500
probability_cutoff = 0.7
margin_cutoff = 1.5
use_DSDB = False   # A boolean based on whether use dataset from DSDB or sstored on the disk
# create a dictionary to keep the parameters organized
classifier_param = dict({'training_set':training_set, 'training_full':training_full, 
                                       'dLevel':dLevel, 'study_annotation':study_annotation,
                         'min_genes':min_genes, 'probability_cutoff':probability_cutoff, 'margin_cutoff':margin_cutoff, 'use_DSDB':use_DSDB})


# Parameters
genes = [           # Selected genes for visualization       
    "MYC",
    "CTNNB1",
    "CCND1",
    "AXIN2",
    "SOX9",
    "LGR5",
    "TCF4",
    "ASCL2",
    "CD81",
    "GSK3B",
    "APC",
]
controls = ["randg", "sgChr4", "NTC"] # This can have two or more classes, but if it has two, then remove to end up with two
# category = "SCN_class"
x_category = "Untreatedclus_SW1417"
y_category = "DOXclus_SW1417"
remove_gRNA_doublets = True  #True or False unncecessary
permutations = 10000  #Number of permutation
# create a dictionary to keep the parameters organized
Fisher_param = dict({'genes':genes, 'controls':controls, 'x_category':x_category,
                                       'y_category':y_category, 'remove_gRNA_doublets':remove_gRNA_doublets,
                         'permutations':permutations})

pVAL = 0.05 # the p-value cut-off for the significance
n = 500 # subsample number of NTCs, if you don't want subsampleing of NTCs set n=None 
embeddings_key = "X_pca" # "X_pca", "X_scVI", "X_pca_centered", "X_scVI_centered"
# create a dictionary to keep the parameters organized
ED_param = dict({'pVAL':pVAL, 'n':n, 'embeddings_key' : embeddings_key, 'i':1})


for key in ['X_pca', 'X_pca_sphered', 'X_pca_sphered_Norm', 'X_scVI']:
    embeddings_key = key # "X_pca" or "X_scVI"
    transcriptoprint = True #if True only for genes have phenotype if False for all genes
    # create a dictionary to keep the parameters organized
    MV_param = dict({'embeddings_key' : embeddings_key, 'transcriptoprint':transcriptoprint})


    DF = Features(DatasetID, auth, DEV = DEV,
            DS_test = DS_test,
            sources = sources, 
            author = author,
            Notebooks_HOME = Notebooks_HOME,
            Notebooks_Savepath = Notebooks_Savepath,
            Reports_format = Reports_format,
            Dataset_HOME = Dataset_HOME,
            gene_param = gene_param,
            hashing_param = hashing_param,
            crispr_param = crispr_param,
            Norm_param=Norm_param,
            Embeddings_param=Embeddings_param,
            classifier_param = classifier_param,
            Fisher_param = Fisher_param,
            ED_param = ED_param,
            MV_param=MV_param)

    #DF.embeddings()
    #DF.classifier()
    #DF.fisher_test()
    #DF.e_dist()
    #DF.e_dist(perPerturbation=True)
    DF.multivariate()
    #DF.mixscape()
    #DF.Hotspot()