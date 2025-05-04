#' trainDeepIP trains DeepIP model for given input sequences
#' @param condaEnv the python conda environment
#' @param inTrainSeq fa file containing positive (/1, _1, or :1) and negative (end with /0, _0, or :0) sequences
#' @param outTrainedModel output hdf5 file
#' @param epoch iteration number, default is 100,
#' @param seqLabel define the label of the sequences in `inTrainSeq`, could be '' or 0/1, 1/0, 0:1, 1:0, 01, 10
#' @return outTrainedModel file name.
#' @examples
#' \dontrun{
#' trainDeepIP("C:/Users/LENOVO/anaconda3/envs/DeepPASTA/",
#'             'train.10000T.10000F/train.10000T.10000F.1.fa',
#'             'train.10000T.10000F/train.10000T.10000F.1.epoch100.hdf5',
#'             100)
#' }
#' @export
trainDeepIP<-function(condaEnv, inTrainSeq, outTrainedModel, epoch=100, seqLabel='') {
  DeepIP_TRAIN <- system.file("py", "DeepIP_train.py", package = "DeepIP")
  reticulate::use_condaenv(condaEnv)

  ## train: train.10000T.10000F.1.fa (~90s/epoch)
  cmd=sprintf("trainSeq='%s'; trainedModel='%s'; epoch=%d; seqLabel='%s'",
              inTrainSeq,
              outTrainedModel,
              epoch,
              seqLabel)
  reticulate::py_run_string(cmd)
  reticulate::py_run_file(DeepIP_TRAIN)
  return(outTrainedModel)
}


#' testDeepIP test given input sequences using given trained model
#' @param condaEnv the python conda environment
#' @param inTestSeq fa file for test, could have labels like (/1, _1, or :1) or (/0, _0, or :0) in the sequence title or not.
#' @param inTrainedModel trained hdf5 file
#' @param outTestCsv output csv file with five columns: title,score,true_label,predict_label,res=TP/FP/FN/TN
#' @param seqLabel define the label of the sequences in `inTestSeq`, could be '' or '0/1','0:1','01','1/0','1:0','10','00','0:0','0/0','11','1/1','1:1','1','0'.
#' @return outTestCsv file name.
#' @examples
#' \dontrun{
#' testDeepIP("C:/Users/LENOVO/anaconda3/envs/DeepPASTA/",
#'             'test.all.fa',
#'             'train.10000T.10000F/train.10000T.10000F.1.epoch100.hdf5',
#'             'train.10000T.10000F.1.epoch100_ON_test.all.csv'
#'             '')
#' # metrics
#' statDeepRes('train.10000T.10000F.1.epoch100_ON_test.all.csv',
#'             ofile='train.10000T.10000F.1.epoch100_ON_test.all.stat.csv')
#'
#' # TP/FP fa profile
#' plotFaDeepRes('train.10000T.10000F.1.epoch100_ON_test.all.csv',
#'               fafile='test.all.fa')

#' # grams metrics
#' stats=statDeepResByGrams('train.10000T.10000F.1.epoch100_ON_test.all.csv',
#'          ofile='train.10000T.10000F.1.epoch100_ON_test.all.statByGrams.csv',
#'                          plot=TRUE)
#' }
#' @export
testDeepIP<-function(condaEnv, inTestSeq, inTrainedModel, outTestCsv, seqLabel='') {
  DeepIP_TEST <- system.file("py", "DeepIP_test.py", package = "DeepIP")
  reticulate::use_condaenv(condaEnv)

  cmd=sprintf("testSeq='%s'; trainedModel='%s'; outputFile='%s'; seqLabel='%s'",
              inTestSeq,
              inTrainedModel,
              outTestCsv,
              seqLabel)
  reticulate::py_run_string(cmd)
  reticulate::py_run_file(DeepIP_TEST)
  return(outTestCsv)
}
