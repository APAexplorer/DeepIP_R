## ----include = FALSE-------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(
  fig.width = 6,
  fig.height = 5.5,
  collapse = TRUE,
  warning = FALSE,
  comment = "#>"
)


## ----libs, warning=FALSE, message=FALSE------------------------------------------------------------------------------------------------------
library(DeepIP, warn.conflicts = FALSE, quietly=TRUE)

## genome assembly
library(BSgenome.Athaliana.TAIR.TAIR9, 
        warn.conflicts = FALSE, quietly=TRUE)
bsgenome=BSgenome.Athaliana.TAIR.TAIR9

## genome annotation
library(TxDb.Athaliana.BioMart.plantsmart51, 
        warn.conflicts = FALSE, quietly=TRUE)
txdb=TxDb.Athaliana.BioMart.plantsmart51


## ----data_input------------------------------------------------------------------------------------------------------------------------------
## high quality PAs (HQ) -- for getting real A-rich PAs
hqPAfiles=system.file("demo_data_ath", "DRS_Ath.bed", package = "DeepIP")

## other potential PAs (LQ)
## -- together with HQ PAs and TXDB for getting free-PA regions for IP data.
allPAfiles=system.file("demo_data_ath", "PlantAPAdb_Ath.bed", 
                       package = "DeepIP")

## check chr names in bsgenome
## it is better to make chr names in PA files, TXDB, and bsgenome consistent
## although DeepIP will check and make the chr consistency.
seqnames(bsgenome) 

## used chrs
chrs=paste0('Chr', c(1:5))

## known PA signals for hierarchical screening sequences
## For Arabidopsis, only AATAAA is the most dominant polyA signal, 
## accouting for ~10% PAs.
## If the polyA signal is unknown for the species, simply set grams=NULL.
grams='AATAAA'
gramsPriority=NULL

## for animal species, we can use:
# grams=c('AATAAA','ATTAAA','TATAAA','AGTAAA','AATACA','CATAAA','AATATA',
#           'GATAAA','AATGAA','AATAAT','AAGAAA','ACTAAA','AATAGA','ATTACA',
#           'AACAAA','ATTATA','AACAAG','AATAAG')
# gramsPriority=c(1, 2, rep(3, length(gramsMM)-2))


## output directory
## all files will be generated in this folder
outputDir='D:/DeepIP_data_ath/'

if (!file.exists(outputDir)) dir.create(outputDir)


## ----getArichPAseqsIn3UTR--------------------------------------------------------------------------------------------------------------------

getArichPAseqsIn3UTR(paBedFiles=hqPAfiles, txdb=txdb, extUTRLen=5000,
                     bsgenome=bsgenome, chrs=chrs, grams=grams,
                     gramsPriority=gramsPriority,
                     outputSeqPre=paste0(outputDir, 
                                         "realPA_Arich_seq/realPA_Arich_seq"))


## ----getFreePARangesIn3UTR-------------------------------------------------------------------------------------------------------------------
## The output file utrsNoPA.RDS contains regions without any potential PAs
utrsNoPA=getFreePARangesIn3UTR(paBedFiles=c(hqPAfiles, allPAfiles), 
                         txdb=txdb, extUTRLen=5000, extPA=200, 
                         outputRds='utrsNoPA.RDS', chrs=chrs, 
                         plotFA=TRUE, N=5000, bsgenome=bsgenome,
                         outputDir=outputDir)


## ----getArichIPseqs--------------------------------------------------------------------------------------------------------------------------
## It may take about 0.5H for 90,000 sequences.
## Here for demonstration, we set `N=3000` to speed up, 
## which randomly selects 3000 A-rich sequences for processing.
## In practive, we can set `N=NULL` to process all sequences.
getArichIPseqs(regionsNoPARDS=paste0(outputDir,'utrsNoPA.RDS'), 
               bsgenome=bsgenome, 
               grams=grams,
               N=3000, 
               outputSeqPre=paste0(outputDir, "IP_Arich_seq/IP_Arich_seq"))


## ----statCntFas------------------------------------------------------------------------------------------------------------------------------
## count the number of real and IP sequences
fas=statCntFas(path=paste0(outputDir, 'realPA_Arich_seq/'), 
               filePre='realPA_Arich_seq.')  

fas=statCntFas(path=paste0(outputDir, 'IP_Arich_seq'), 
               filePre='IP_Arich_seq.')


## ----split2TrainTest-------------------------------------------------------------------------------------------------------------------------

dir.create(paste0(outputDir,'modelDataSplits_ath'))

## true train and test files
split2TrainTest(path=paste0(outputDir,'realPA_Arich_seq/'), 
                filePre='realPA_Arich_seq.',
                dir1=paste0(outputDir,
                            'modelDataSplits_ath/realPA_Arich_seq_train_per70'),
                dir2=paste0(outputDir,
                            'modelDataSplits_ath/realPA_Arich_seq_test_per30'), 
                per1=0.7,
                label=":1")

## false train and test files
split2TrainTest(path=paste0(outputDir,'IP_Arich_seq'), 
                filePre='IP_Arich_seq.',
                dir1=paste0(outputDir,
                            'modelDataSplits_ath/IP_Arich_seq_train_per70'),
                dir2=paste0(outputDir,
                            'modelDataSplits_ath/IP_Arich_seq_test_per30'), 
                per1=0.7,
                label=":0")


## ----getNseqs_train--------------------------------------------------------------------------------------------------------------------------
## get the polyA signal distributions for the subsequent random sampling
REALPA_PERC=statCntFas(path=paste0(outputDir, 'realPA_Arich_seq/'), 
                       filePre='realPA_Arich_seq.') 

seqDir=paste0(outputDir,'modelDataSplits_ath/')

## randomly sample 10,000 sequences as real (positive) training data
files=getNseqs(seqDir=paste0(seqDir, 'realPA_Arich_seq_train_per70'), 
               N=10000,
               nsplits=1,
               outputPre=paste0(seqDir, 
                                "realPA_Arich_seq_train_per70_10000s_1splits"),
               perc=REALPA_PERC) 


## plot the real training fa file for validation
movAPA::plotATCGforFAfile(faFiles=files, ofreq=FALSE,
                          opdf=T, refPos=101)

## randomly sample 10,000 sequences as false (IP, negative) training data
## Here for demonstration, we set `N=1000`.
files=getNseqs(seqDir=paste0(seqDir, 'IP_Arich_seq_train_per70'), 
               N=1000, 
               nsplits=1,
               outputPre=paste0(seqDir, 
                                'IP_Arich_seq_train_per70_10000s_1splits'),
               perc=REALPA_PERC) 

## plot the IP training fa file for validation
movAPA::plotATCGforFAfile(faFiles=files, ofreq=FALSE,
                          opdf=T, refPos=101)


## ----getNseqs_test---------------------------------------------------------------------------------------------------------------------------

## randomly sample 5,000 sequences as real (positive) test data
files=getNseqs(seqDir=paste0(seqDir, 'realPA_Arich_seq_test_per30'), 
               N=5000, 
               nsplits=1,
               outputPre=paste0(seqDir, 
                                'realPA_Arich_seq_test_per30_5000s_1splits'),
               perc=REALPA_PERC) 


movAPA::plotATCGforFAfile(faFiles=files, ofreq=FALSE, opdf=T, refPos=101)

## randomly sample 5,000 sequences as IP (negative) test data
## Here for demonstration, we set `N=500`.
files=getNseqs(seqDir=paste0(seqDir, 'IP_Arich_seq_test_per30'), 
               N=500, 
               nsplits=1,
               outputPre=paste0(seqDir, 
                                'IP_Arich_seq_test_per30_5000s_1splits'),
               perc=REALPA_PERC) 

movAPA::plotATCGforFAfile(faFiles=files, ofreq=FALSE, opdf=T, refPos=101)



## ----combineFaFiles--------------------------------------------------------------------------------------------------------------------------
combineFaFiles_fast(paste0(seqDir, 
                    c('realPA_Arich_seq_train_per70_10000s_1splits.split.1.fa',
                      'IP_Arich_seq_train_per70_10000s_1splits.split.1.fa')),
                    ofile = paste0(seqDir, 'train.10000T.10000F.fa'),
                    verbose = TRUE)


combineFaFiles_fast(paste0(seqDir, 
                    c('realPA_Arich_seq_test_per30_5000s_1splits.split.1.fa',
                      'IP_Arich_seq_test_per30_5000s_1splits.split.1.fa')),
                    ofile = paste0(seqDir, 'test.5000T.5000F.fa'),
                    verbose = TRUE)



## ----conda, eval=FALSE-----------------------------------------------------------------------------------------------------------------------
# 
# # create the conda env in the command window
# conda create -n DeepIP python=3.7
# 
# # test the env
# conda env list
# conda activate DeepIP
# 
# # install the following modules
# pip install Keras
# pip install tensorflow
# pip install pandas
# pip install sklearn
# 
# # or use conda install, like
# conda install -c anaconda protobuf


## ----train, eval=FALSE-----------------------------------------------------------------------------------------------------------------------
# setwd(seqDir)
# trainDeepIP(condaEnv="A_CONDA_ENV",
#             inTrainSeq='train.10000T.10000F.fa',
#             outTrainedModel='train.10000T.10000F.epoch100.hdf5',
#             epoch=100)
# 


## ----test, eval=FALSE------------------------------------------------------------------------------------------------------------------------
# testDeepIP(condaEnv="A_CONDA_ENV",
#            inTestSeq='test.5000T.5000F.fa',
#            inTrainedModel='train.10000T.10000F.epoch100.hdf5',
#            outTestCsv='train.10000T.10000F.epoch100_ON_test.5000T.5000F.csv',
#            seqLabel='')
# 
# # calculate ROC/AUC/F1/... metrics
# statDeepRes('train.10000T.10000F.epoch100_ON_test.5000T.5000F.csv',
#             ofile='train.10000T.10000F.epoch100_ON_test.5000T.5000F.stat.csv')
# 
# # Plot single nucleotide profiles for TP/FP/TN/FN sequences
# plotFaDeepRes('train.10000T.10000F.epoch100_ON_test.5000T.5000F.csv',
#               fafile='test.5000T.5000F.fa')

