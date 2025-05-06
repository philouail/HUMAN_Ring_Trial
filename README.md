# HUMAN_Ring_Trial

This repo is to group the analysis of the Ring trial study which is part of the
HUMAN doctoral network main studies.


The study is defined as such: 

Within HUMAN, we would like to understand better where the differences between
different LC-MS setups used for metabolomics come from. Therefore, different
setups are compared using different methods. First, each participating
laboratory uses a specific method in “daily” laboratory business. Second, a
standardized method using the same column and gradient is used. Therefore, a
total of 83 mixtures will be analyzed twice: by a “HUMAN reference method” and
by a lab specific method.  In the first step for comparison of the methods and
laboratories, mixtures obtained from the MetaSci metabolite standard library are
measured. Since the ground truth (number and identity of metabolite) is known,
we can compare the different chromatographic methods within a laboratory and
compare the different LC-MS setups from different laboratories using the
standardized method.

The first step of this would be to automatise the library building as much as possible. 

each lab will have their own preprocessing files which will output a certain number of 
object used in the library building. 

The library file however will be common to all. and will generate 2 csv files as well as a number of plot to help in library building
and making informed choice in case of ambiguity. 