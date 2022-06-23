#GBLUP implemented with either one or two kinship matrices
#This script is supposed to be run via command line, e.g.:
#> Rscript multikinship_regression.R <config>
#Where <config> is a configuration file describing the 
#regression to be tested.
#Please note that you need the input data in the proper folder
#structure to be able to run this script.

# KEY LIBRARIES -----------------------------------------------------------
library(GROAN)

# INPUT CONFIGURATION MANAGEMENT ------------------------------------------
args = commandArgs(trailingOnly=TRUE)
if (length(args) == 1){
  #loading the parameters
  source(args[1])
}else{
  stop('a config file is needed')
}

# REGRESSION --------------------------------------------------------------
for (i in 1:nrow(config)){
  #current configuration
  curr = config[i,]
  print(curr)

  #reading phenotypes
  phenos_infile = file.path(curr$base_folder, 'phenotypes', curr$phenotype_filename)
  stopifnot(file.exists(phenos_infile))
  phenos = read.csv(phenos_infile, stringsAsFactors = FALSE)

  #reading kinship matrices
  kinships = list()
  for (k in strsplit(curr$kinships, split = ' ')[[1]]){
    kinship_infile = paste(sep='', 
                           'kinship_', k,
                           '_minMAF', curr$MAF_min, 
                           '_maxMAF', curr$MAF_max, '.csv.gz')
    kinship_infile = file.path(curr$base_folder, kinship_infile)
    stopifnot(file.exists(kinship_infile))
    kinships[[k]] = read.csv(kinship_infile, stringsAsFactors = FALSE, row.names = 1)
  }

  #first kinship matrix will be used as covariance, from the second one going
  #on (if present) they will be used for extracovariates
  covariance_curr = kinships[[1]]
  extraCovariates_curr = NULL
  if (length(kinships) > 1){
    for (k in 2:length(kinships)){
      if (is.null(extraCovariates_curr)){
        extraCovariates_curr = kinships[[k]]
      }else{
        extraCovariates_curr = cbind(extraCovariates_curr, kinships[[k]])
      }
    }
    #this prevents problems down the line: all covariates gets a unique name
    colnames(extraCovariates_curr) = paste(sep='', 'k', 1:ncol(extraCovariates_curr))
  }

  #taking note of the available kinships
  kinship_names = paste(collapse='___', names(kinships))

  #room for results
  outfolder = file.path(curr$base_folder, 'regressions', kinship_names)
  dir.create(outfolder, showWarnings = FALSE, recursive = TRUE)

  #building GROAN structures
  ds_name = paste(sep = '',
                  'dataset=', basename(curr$base_folder),
                  ' phenotype=', curr$phenotype,
                  ' minMAF=', curr$MAF_min,
                  ' maxMAF=', curr$MAF_max,
                  ' kinship=', kinship_names,
                  ' regressor=', curr$regressor
  )
  ds = createNoisyDataset(name = ds_name,
                          covariance = covariance_curr, extraCovariates = extraCovariates_curr,
                          phenotypes = scale(phenos[,curr$phenotype]))
  wb = createWorkbench(folds = curr$folds, reps = curr$reps, outfolder = outfolder,
                       regressor = phenoregressor.BGLR.multikinships , regressor.name = 'RKHS_multikin',
                       nIter = curr$nIter, burnIn = curr$burnIn)

  #running regression
  GROAN.run(ds, wb)
}

