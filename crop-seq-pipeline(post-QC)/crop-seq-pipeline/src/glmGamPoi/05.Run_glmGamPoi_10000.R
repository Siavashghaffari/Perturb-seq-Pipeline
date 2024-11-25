######## Load all libraries 


## Load libraries
library(glmGamPoi)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(parallel)
library(argparse)

##########################################################################################################
#############################        Part 1: Load Arguments        #######################################
##########################################################################################################

## Load arguments
parser <- ArgumentParser()
parser$add_argument('-i', '--input_directory', help='path to folders where by target data is saved (does not include the by target extension)', required = TRUE)
parser$add_argument('-g', '--genes_file', help='path to file naming genes (features) of target and NTCs', required = TRUE)
parser$add_argument('-t','--target', required =TRUE,
                    help='name of target to compare to')
parser$add_argument('-c','--control', 
                    help='name of control to compare target to. Default is NTC', default = "NTC")
parser$add_argument('-o','--output_folder', 
                    help='name of folder where to store results', required = TRUE)

## Print loaded arguments
args <- parser$parse_args()
print(args)



##########################################################################################################
########################        Part 2: Load Data by Target        #######################################
##########################################################################################################

## Base Directory
base_dir = paste(args$input_directory)

## File Names - Control vs Target
control_file = paste(base_dir, "by_target/", "count_mat",args$control,".mm", sep = "")
target_file = paste(base_dir, "by_target/", "count_mat",args$target,".mm", sep = "")
control_cellnames_file = paste(base_dir, "by_target/", "count_mat",args$control,"_rownames.csv", sep = "")
target_cellnames_file = paste(base_dir, "by_target/", "count_mat",args$target,"_rownames.csv", sep = "")
genes_file = args$genes_file
sample_metadata_file = paste(base_dir, "coldata_for_glmgampois.csv", sep = "")

## Load Control and Target Data
mat.control =  t(readMM(control_file))
mat.target =  t(readMM(target_file))

## Load information about samples aka cell names
cell_names.control = fread(control_cellnames_file)
cell_names.target = fread(target_cellnames_file)

## Information about genes
gene_names = fread(genes_file, header = T)

## Load bdev Genes
sublib1_bdev =fread("/gstore/scratch/u/ghaffars/glmGamPoi/bdev/nopos_controls_sublib1_top3000_bdev.csv")
sublib2_bdev =fread("/gstore/scratch/u/ghaffars/glmGamPoi/bdev/nopos_controls_sublib2_top3000_bdev.csv")
sublib3_bdev =fread("/gstore/scratch/u/ghaffars/glmGamPoi/bdev/nopos_controls_sublib3_top3000_bdev.csv")
sublib4_bdev =fread("/gstore/scratch/u/ghaffars/glmGamPoi/bdev/nopos_controls_sublib4_top3000_bdev.csv")
hvg_union = unique(c(sublib1_bdev$Symbol, sublib2_bdev$Symbol, sublib3_bdev$Symbol, sublib4_bdev$Symbol))

## Cell Sample Meta Data
dat.meta = fread(sample_metadata_file)

## Data Sample Names
control_file = paste(base_dir, "by_target/", "count_mat",args$control,".mm")
target_file = paste(base_dir, "by_target/", "count_mat",args$target,".mm")

## Add rownames and colnames
rownames(mat.control) = rownames(mat.target) = gene_names$Symbol
colnames(mat.control) = cell_names.control$samples
colnames(mat.target) = cell_names.target$samples

## Print info about data
print('Control mat size')
print(dim(mat.control))
print('Target mat size')
print(dim(mat.target))

## Keep only highly variable genes
mat.control = mat.control[rownames(mat.control) %in% hvg_union,]
mat.target = mat.target[rownames(mat.target) %in% hvg_union,]

## Print info about data
print("Keep only HVG")
print('Control mat size')
print(dim(mat.control))
print('Target mat size')
print(dim(mat.target))

## Combine Target and Contol matrix
mat.combined <- cbind(mat.target, mat.control)
print("done loading data")

##########################################################################################################
#############################     Part 3: Filter Based on Meta Data       #################################
##########################################################################################################

## Cells to Keep
cells.keep.fromcountmat <- colnames(mat.combined)

## Load meta data and keep only those cells relevatnt to control and target
dat.meta <- dat.meta %>%
  dplyr::filter(cell %in% cells.keep.fromcountmat) %>% as.data.frame()

# Remove gems without any target cells and add rownames
dat.meta <- dat.meta %>%
  group_by(batchid) %>%
  mutate(gem.has.target = length(label[label == 'zAll']) > 0) %>%
  dplyr::filter(gem.has.target) %>%
  as.data.frame()
rownames(dat.meta) <- dat.meta$cell

## Size of Meta Data
print('Meta quality size')
print(dim(dat.meta))

## Examine several rows
print('Peek at meta')
print(head(dat.meta))

## Cells to keep (aka those that are in a gem with at least one target sample)
cells.keep.gemfilt <- dat.meta$cell 

## Print info about size of matrix
print('dim before filtering out gems with no target cells')
print(dim(mat.combined))
print('dim after filtering out gems with no target cells')
mat.combined <- mat.combined[, cells.keep.gemfilt]
print(dim(mat.combined))

print("done with filtering data")

##########################################################################################################
######################.         Part 4: Prepare for Running glmGamPoi             ########################
##########################################################################################################

## Re order cells of mat combined to match meta
cells.order <- rownames(dat.meta)
print('Order cells by meta')
print(head(cells.order))
mat.combined <- as.matrix(mat.combined)[, cells.order]

## Stop if the row names of dat meta do not match column names of mat combined
assertthat::assert_that(identical(rownames(dat.meta), colnames(mat.combined)))

## Define Fit Gam Poi Function
FitGlmGamPoi <- function(countmat, meta, jform, jparam = 'label'){
  fit.ggp <- glmGamPoi::glm_gp(data = countmat, design = jform, col_data = meta, on_disk = FALSE)
  pnames <- colnames(fit.ggp$Beta)
  jparam.out <- grep(paste0("^", jparam), pnames)
  de.ggp <- test_de(fit.ggp, contrast = labelzAll)
  identity_design_matrix <-  diag(nrow = ncol(fit.ggp$Beta))
  pred <- predict(fit.ggp, se.fit = TRUE, newdata = identity_design_matrix)
  colnames(pred$fit) <- pnames
  colnames(pred$se.fit) <- pnames
  dat.out <- data.frame(gene = rownames(pred$fit), beta = pred$fit[, jparam.out], beta.se = pred$se.fit[, jparam.out], stringsAsFactors = FALSE)
  out.lst <- list(de.ggp = de.ggp, beta = fit.ggp$Beta, dat.out = dat.out, meta = meta)
  return(out.lst)
}

## Form of equation to fit, in this case we correct for batchid/gem
jjform <-  as.formula(~ label + batchid)
print(paste('model:', as.character(jjform), collapse = " "))

print("done with set up")

##########################################################################################################
###########################    Part 5: Run glmGamPoi and Save Results       ##############################
##########################################################################################################

## Run glmGamPoi
system.time(
  out.lst <- FitGlmGamPoi(countmat = mat.combined, meta = dat.meta, jform = jjform, jparam = 'label')
)

## Create output directory if not exist
glm_dir = paste(args$output_folder, "unweightedglmGamPoi_gembatch_results/",sep = "")
if(!file.exists(glm_dir)){ 
  dir.create(glm_dir)
}
      
## Save Results            
saveRDS(out.lst, file = paste(glm_dir,"glm_unweighted_target_",args$target,"_top10000_sublib124union_genes.rds", sep = ""))
print("all done")

