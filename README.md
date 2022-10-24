# Read me
#work
### Running isres_plus
isres_plus takes three arguments to run.  The first argument is the function handle to be optimized. The second argument contains the lower and upper limits of the unknown model parameters. The function takes a third argument called options, which is a structure containing three sub-structures called “evo”, “plus”, and “model”. The “evo” structure contains algorithm parameters that define the evolutionary part of the algorithm while the “plus” structure contains controls for the linstep and Newton step part of the algorithm.  The “model” substructure can be used to send in model-specific parameters to the optimizing function. isres_plus.m contains default values of the “evo” and “plus” parts of the options structure in case none of the values are supplied. 

### Other details. 
Use the sc_isres_plus.m script file to run the dl/cact, gap gene circuit, and smad signaling models and as a guide to run any model of your own. isres_plus.m additionally implements an island evolution strategy for the evolutionary part of the alogrithm as per [Fomekong-Nanfack, 2007] (https://academic.oup.com/bioinformatics/article/23/24/3356/262640). However, isres_plus does not run it by default.  The default values of the evolutionary part of isres_plus are those recommended by [Runarsson & Yao, 2005](https://ieeexplore.ieee.org/document/1424197)). 
