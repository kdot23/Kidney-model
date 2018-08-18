# Kidney-model
------------------------------

## This is a repository housing our experimental code for doing online matching with Kidney Exchange.
## Files are organized as follows:
---------------------------------

KidneyDataGen - Generates a sample population to experiment on
  * -T number of compatible pairs
  * -K number of incompatible pairs
  * -o name of output file

generateBeta - Generates beta values as training data for sample population
  * --inputFiles list of inputfiles to run on
  * -o name of output file
  * --quality should be present if optimizing for QUALITY instead of COUNT
  * --graph_state should be present if betas predicted by LP on incompatible pool should be included as training data

greedyMatching - performs online matching using the greedy algorithm
  * --inputFiles list of files to perform matching on
  * -o name of output file
  * --quality should be present if optimizing for QUALITY instead of COUNT
  * --agents name of output file containing agent information
  
onlineMatching - performs online matching using ODASSE with various methods for estimating beta
  * --trainFiles list of training files to use
  * --testFiles list of files to perform matching on
  * -o name of ouptut file
  * --agents name of file to output agent information
  * -v list of variables to use as training data
  * --quality should be present if optimization should be done for QUALITY instead of COUNT
  * --forestRegression should be present if forestRegression is to be used instead of linear
  * --lpEstimator should be present if betas from LP on incompatible pool should be used as predicted betas
  * --lpRepeat like lpEstimator but reruns LP after each match
  * --graph_state flag should be present if betas predicted by LP on incompatible pool should be used as training data

oracleMatching - performs online matching assuming we know the future
  * --inputFiles list of inputfiles to run on
  * --quality flag should be present if optimization should be run for QUALITY instead of COUNT
  * -o name of output file
  * --agents name of file agent info should be outputted to

2Cycle and 3Cycle matching are in their corresponding directories
### Note: Make sure you add the --quality flag to GenerateBeta, GreedyMatching and OnlineMatching as optimizing for COUNT in the online setting does very poorly
