# polyandry_meta_analysis
Raw data and code used to produce figures and analyses for the manuscript titled "The adaptive significance of polyandry: a meta-analysis"

## Data folder
The polyandry_meta_full_dataset.csv file contains the full set of extracted data from eligible studies. We then split this dataset into a fecundity and a longevity data file using the cleaning_data.R file that can be found in the scripts folder. The resulting fecundity and longevity datasets are also included in the data folder. Finally, the extracted model results from the global fecundity model (no moderators) and global longevity model are included in the globel_model_results.csv file.

## Scripts folder
This folder contains the cleaning_data.R script which separates the full dataset into fecundity and longevity data. It also contains code that generated the global model results figure (Fig. 2 in the manuscript). 

All of the code used to produce the fecundity analyses and figures can be found in the fecund_analysis.R script while all of the code used to produce the longevity analyses and figures can be found in the longevity_analysis.R script. 

