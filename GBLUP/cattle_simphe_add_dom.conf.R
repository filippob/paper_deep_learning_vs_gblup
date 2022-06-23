#This config file generates the configuration for testing GBLUP on additive
#matrix on all the simphe phenotypes

#We will be filtering on:
QTN = 1000

#these are the common parameters, valid for all configurations
baseconfig = data.frame(
  base_folder = '~/research/deep_learning_for_breeding/datasets/cattle/',
  kinships = paste('additive', 'dominance'),
  MAF_min = 0.05,
  MAF_max = 0.5,
  regressor = 'RKHS',
  folds = 5,
  reps = 10,
  phenotype = NA, #this needs to vary
  phenotype_filename = 'phenotypes_iteration3.csv',
  nIter = 10000,
  burnIn = 500
)

#we are now reading the list of available phenotypes from the proper pheno file
phenos = read.csv(file.path(baseconfig$base_folder, 'phenotypes', baseconfig$phenotype_filename))
phenos = colnames(phenos)[grepl(pattern = '^simphe_', x = colnames(phenos))]

#keeping only phenotypes related to the current QTN
tmp = paste(sep='', '_QTN', QTN, '_')
phenos = phenos[grepl(phenos, pattern = tmp, fixed = TRUE)]

#keeping only phenotypes with no epistasis
phenos = phenos[grepl(phenos, pattern = '_AA0_AD0_DA0_DD0_', fixed = TRUE)]

#creating the config object, used in the baseline_regression.R script
config = NULL
for (p in phenos){
  tmp = baseconfig
  tmp$phenotype = p
  config = rbind(config, tmp)
}

#cleanup, so that no variables are declared except from config
rm(tmp, baseconfig, p, phenos, QTN)
