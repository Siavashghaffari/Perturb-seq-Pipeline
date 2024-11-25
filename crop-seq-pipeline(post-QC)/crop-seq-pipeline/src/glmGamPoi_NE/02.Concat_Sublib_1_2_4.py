######## Load all libraries 

##########################################################################################################
########################### Part 1: Load Libraries #######################################################
##########################################################################################################

#single cell libraries
import scanpy as sc
import anndata as ad


##########################################################################################################
########################### Part 2: Load Data ############################################################
##########################################################################################################



sublib1_no_pos_controls_out_final = sc.read_h5ad("/gstore/project/elenakharitonova_melocars/downsampling_perturbseq/Pipeline_Development/data/Sublib1_DS000017114_remove_pos_cont_counts_obs_var.h5ad")
print("loaded sublib1")

# Load sublibrary 2 and 4 no positive controls
sublib2_4 = sc.read_h5ad("/gstore/project/elenakharitonova_melocars/downsampling_perturbseq/Pipeline_Development/data/Sublib2_Sublib4_DS000016652_remove_pos_cont_counts_obs_var.h5ad")
print("loaded sublib2 and 4")


##########################################################################################################
########################### Part 3: Concatanate ##########################################################
##########################################################################################################


## Concat sublib 1, 2 and 4
adata = ad.concat([sublib1_no_pos_controls_out_final, sublib2_4])
print("concatenated data")



##########################################################################################################
########################### Part 4: Save #################################################################
##########################################################################################################

## Save to h5ad file
adata.write_h5ad("/gstore/project/elenakharitonova_melocars/downsampling_perturbseq/Pipeline_Development/data/Sublib1_Sublib2_Sublib4_DS000017114_DS000016652_remove_pos_cont_counts_obs_var.h5ad")