# FIT_FecalContamination

These are the core functions needed to 

1) Run the code to obtain the results for the manuscript, "Characterizing differences in sources and contributions to fecal contamination of sediment and surface water with the microbial FIT framework"
(associated with the ID es-2022-00224x)

This can be run using the code, *RunPaperResults.m*

2) Run an example of the "Find" stage of FIT

This can be run using the code, *FindExample.m*

3) Run examples of the "Inform" stage of FIT with the Overland and River distance with Flow (ORF) model and the Ground hauling, Overland and River distance with Flow (GORF) spatial predictor models 

This can be run using the code, *InformExample.m*

The original source data can be accessed publically via the sites described in the supporting information of the manuscript (See section S6). 
For data that is not publically accessible, you may request it from the organizations described in the supporting information of the manuscript (See section S6).

Here, we only include the code for obtaining the results of the manuscript. However, to run this code, data needs to be downloaded. The data can be downloaded from https://figshare.com/projects/Characterizing_differences_in_sources_and_contributions_to_fecal_contamination_of_sediment_and_surface_water_with_the_microbial_FIT_framework_DATA/133065

# MATLAB Toolboxes that are needed to run the code: 
1) Optimization Toolbox
2) Global Optimization Toolbox
3) Bioinformatics Toolbox

Also, I will note here that these functions have been tested in MATLAB version 2020b and up. 

# Different versions of this code: 
## version01 included the following functions 
coord2dist.m
FindExample.m
FindSitesonNetwork.m
GORFtinypartv2.m
InformExample.m
KewauneeObjectiveFunction1.m
KewauneeObjFunLoad.m
KewauneePenaltyFunction1.m
KewauneePenaltyFunction2.m
KewRoadDistanceMatrix.m
PrecipitationObjectiveFunctionV2.m
ReliabilityScore.m
RunPaperResults.m
sectionedcoord2dist.m
sedc2.m
sedc3.m
sedc4.m
SEDCparamOptim.m 

## version02 
The file name being called to run the core functions was changed to "FITBacteroides.mat" for...
FindExample.m
InformExample.m
RunPaperResults.m 

The following code was missing and was added: 
CVlinmdl2.m
randinterval.m

It was noted that the following MATLAB Toolboxes were required to run the code: 
1) Optimization Toolbox
2) Global Optimization Toolbox

## version03 
The following code was missing and was added: 
KewauneeSEDCresults.m 

It was noted that the following MATLAB Toolboxes were required to run the code: 
1) Bioinformatics Toolbox 






