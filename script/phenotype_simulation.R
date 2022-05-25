#install.packages('SimPhe')
library(SimPhe)

#steps to success
# - decide the number of markers for the main effects (additive, dominant)
# - those same markers will be used, in pairs, for epistasis
# - decide trait heritability (influences added final noise)
# - decide overall mean (not strictly enforced by simphe)
# - decide a plan for effect sizes

# INPUT PARAMETERS --------------------------------------------------------

#dataset folder, change to your own path
root = '~/research/paper_deep_learning_vs_gblup'

#markers file, expected in 0/1/2, markers x samples
SNP_file = file.path(root, 'data', 'cattle_SNPs', 'SNPs.csv.gz')

#output file
pheno_file = file.path(root, 'simulated_phenotypes', 'phenotypes_iteration3.csv')

#overall mean
pheno_mean = 0

#overall heritability
pheno_h2 = 0.7

#coefficient of variation of QTN
QTN_cv = 0.1

#number of QTN in the main effect (must be even, because epistasis
#works in pairs of SNPs)
QTN_num = 1000

#mean absolute effect of QTN, one element per simulated scenario
#actual final effect has 50% chance of being negative,
#i.e. being multiplied by -1
QTN_additive_mean            = c(100,  75,  50, 25,    0,  33,  33,  33,   33)
QTN_dominance_mean           = c(  0,  25,  50, 75,  100,  33,  33,  33,   33)
QTN_additive_additive_mean   = c(  0,   0,   0,  0,    0,  34,   0,   0,   11.3)
QTN_additive_dominance_mean  = c(  0,   0,   0,  0,    0,   0,  34,   0,   11.3)
QTN_dominance_additive_mean  = c(  0,   0,   0,  0,    0,   0,   0,   0,   0)
QTN_dominance_dominance_mean = c(  0,   0,   0,  0,    0,   0,   0,  34,   11.3)

#tmp folder: for SimPhe interface I'll need to write there genos and parameters
#files, and then read everything again...
tmp_folder = tempdir(check = TRUE)
SNP_file_simphe = file.path(tmp_folder, 'genotypes.txt')
param_file_simphe = file.path(tmp_folder, 'params.txt')

# SUPPORT FUNCTIONS -------------------------------------------------------
pick_main_effects = function(SNP_names, QTN, A_mean, D_mean, cv){
  #sampling the markers with an effect
  SNP_selected = sample(x = SNP_names, size = QTN, replace = FALSE)

  #generating effects: additive and dominant
  A_eff = rnorm_posneg_effect(n = QTN, mean = A_mean, cv = cv)
  D_eff = rnorm_posneg_effect(n = QTN, mean = D_mean, cv = cv)
  
  #putting all together
  return(data.frame(
    SNP = SNP_selected,
    additive = A_eff,
    dominance = D_eff
  ))
}

#picks from the normal distribution with specification of mean and coefficient
#of variation, plus results have 50% chance to have the sign flipped
rnorm_posneg_effect = function(n, mean, cv){
  res = rnorm(n = n, mean = mean, sd = mean * cv)
  sel = runif(n = n) > 0.5
  res[sel] = res[sel] * -1
  return(res)
}

#SNP with epistasis effect are the same as those with the main effect,
#picked in pairs
pick_epistasis_effects = function(SNP_names, AA_mean, AD_mean, DA_mean, DD_mean, cv){
  QTN_epi = length(SNP_names) / 2
  SNPA_selected = SNP_names[1:QTN_epi]
  SNPB_selected = SNP_names[(QTN_epi + 1):(2 * QTN_epi)]
  
  #generating effects: combinations of additive and dominant
  AA_eff = rnorm_posneg_effect(n = QTN_epi, mean = AA_mean, cv = cv)
  AD_eff = rnorm_posneg_effect(n = QTN_epi, mean = AD_mean, cv = cv)
  DA_eff = rnorm_posneg_effect(n = QTN_epi, mean = DA_mean, cv = cv)
  DD_eff = rnorm_posneg_effect(n = QTN_epi, mean = DD_mean, cv = cv)
  
  #putting all together
  return(data.frame(
    SNPA = SNPA_selected,
    SNPB = SNPB_selected,
    additive_additive = AA_eff,
    additive_dominance = AD_eff,
    dominance_additive = DA_eff,
    dominance_dominance = DD_eff
  ))
}

# GENOTYPES PREP ----------------------------------------------------------
#loading genotypes, bringing them to the SimPhe format
writeLines('Preparing SNP data...')
SNP = read.csv(SNP_file, row.names = 1, stringsAsFactors = FALSE)
SNP = data.frame(t(SNP))
writeLines('DONE')

# PHENOTYPE SIMULATION ----------------------------------------------------
#epoch for this run (so that we can have multiple repetitions of the
#current config and still have unique phenotype names)
epoch = as.integer(Sys.time())

#number of main QTN must be even
if (QTN_num %% 2 == 1){
  stop('Number of QTN must be even')
}

#if there are already phenotypes, we are going to add to the table
pheno_values = NULL
if(file.exists(pheno_file)){
  pheno_values = read.csv(pheno_file, stringsAsFactors = FALSE)
}

#for each requested scenario
for (i in 1:length(QTN_additive_mean)){
  #building a string that will be used to name this phenotype and for interface purposes
  phenoname = paste(sep='', 'simphe',
                    '_mean', pheno_mean,
                    '_hSquare', pheno_h2,
                    '_cv', QTN_cv,
                    '_QTN', QTN_num,
                    '_A', QTN_additive_mean[i],
                    '_D', QTN_dominance_mean[i],
                    '_AA', QTN_additive_additive_mean[i],
                    '_AD', QTN_additive_dominance_mean[i],
                    '_DA', QTN_dominance_additive_mean[i],
                    '_DD', QTN_dominance_dominance_mean[i],
                    '_epoch', epoch
  )
  writeLines(paste(sep='', ' - ', phenoname))
  
  #-------- building the param file - mean part
  fp = file(param_file_simphe, open = 'w')
  param = c(
    #----mean
    '[P1mean]',
    'mean',
    paste(pheno_mean),
    ''
  )
  writeLines(param, con = fp)
  close(fp)
  
  #-------- building the param file - main part
  fp = file(param_file_simphe, open = 'a')
  param = c(
    '[P1main]',
    'SNP additive dominance'
  )
  writeLines(param, con = fp)
  close(fp)

  main_effects = pick_main_effects(SNP_names = colnames(SNP), 
                    QTN = QTN_num, 
                    A_mean = QTN_additive_mean[i],
                    D_mean = QTN_dominance_mean[i],
                    cv = QTN_cv
                    )
  write.table(main_effects, file = param_file_simphe, sep = '\t', row.names = FALSE, quote = FALSE, append = TRUE, col.names = FALSE)
  
  #-------- building the param file - epistasis part
  fp = file(param_file_simphe, open = 'a')
  param = c(
    '',
    '',
    '[P1epistasis]',
    'SNPA SNPB additive_additive additive_dominance dominance_additive dominance_dominance'
  )
  writeLines(param, con = fp)
  close(fp)
  epi_effects = pick_epistasis_effects(SNP_names = main_effects$SNP, 
                                   AA_mean = QTN_additive_additive_mean[i],
                                   AD_mean = QTN_additive_dominance_mean[i],
                                   DA_mean = QTN_dominance_additive_mean[i],
                                   DD_mean = QTN_dominance_dominance_mean[i],
                                   cv = QTN_cv
  )
  write.table(epi_effects, file = param_file_simphe, sep = '\t', row.names = FALSE, quote = FALSE, append = TRUE, col.names = FALSE)
  
  #-------- building the param file - heritability part
  fp = file(param_file_simphe, open = 'a')
  param = c(
    '',
    '[P1heritability]',
    'heritability',
    pheno_h2
  )
  writeLines(param, con = fp)
  close(fp)
  
  #------- invoking simphe
  #here we need a safeguard: sometimes sim.phe fails with warning:
  #In rnorm(nrow(geno), sd = sqrt(exp.noise.var)) : NAs produced
  #and returns NAs
  maximum_cycles = 3
  while(TRUE){
    phenotypes = sim.phe(
      sim.pars = param_file_simphe, 
      fgeno = SNP, 
      ftype = "ind.head", 
      seed = floor(runif(1) * 100000), #this is so each invocation is different 
      fwrite = FALSE, 
      pattern = "[[:alpha:]]+",
      genetic.model	 = 'epistasis')
    
    #if data is good let's break the cycle
    if (sum(is.na(phenotypes[,1])) == 0){
      break
    }
    
    #but we are not doing this more than 100 times
    maximum_cycles = maximum_cycles - 1
    if (maximum_cycles == 0){
      stop('sim.phe keeps returning NAs, stopping')
    }
    
    #if we get here we have to cycle again and try another round of generation
    writeLines('sim.phe returned NAs, trying again')
  }
  
  #putting phenotype with already present ones
  newcolnames = c(colnames(pheno_values), phenoname)
  if (is.null(pheno_values)){
    pheno_values = phenotypes
  }else{
    pheno_values = cbind(pheno_values, phenotypes)  
  }
  colnames(pheno_values) = newcolnames
  
  stop('here')
  
  #saving, overwriting at each iteration
  write.csv(pheno_values, pheno_file, row.names = FALSE)
}