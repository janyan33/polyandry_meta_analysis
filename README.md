# polyandry\_meta\_analysis

Raw data and code used to produce figures and analyses for the manuscript titled "The adaptive significance of polyandry: a meta-analysis"


## Data folder

The polyandry\_meta\_full\_dataset.csv file contains the full set of extracted data from eligible studies. We then split this dataset into a fecundity and a longevity data file using the cleaning\_data.R file that can be found in the scripts folder. The resulting fecundity and longevity datasets are also included in the data folder. Finally, the extracted model results from the global fecundity model (no moderators) and global longevity model are included in the globel\_model\_results.csv file.


Scripts folder
---

This folder contains the cleaning\_data.R script which separates the full dataset into fecundity and longevity data. It also contains code that generated the global model results figure (Fig. 2 in the manuscript).

All of the code used to produce the fecundity analyses and figures can be found in the fecund\_analysis.R script while all of the code used to produce the longevity analyses and figures can be found in the longevity\_analysis.R script.


Other files
---

I've also included the PRISMA flowchart in the main folder and an Excel file with all studies that were read in detail during the full text screen. The first tab includes all 422 which were screened for eligibility with red highlighted rows denoted excluded studies and green highlighted rows denoting included studies. The second tab in this file is a list of all the excluded studies with rationales for each exclusion and the third tab is a full list of only the included studies.

