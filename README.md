# Cognitive_Control_Developmental_Trajectory
This repository includes data and codes related to the manuscript: 
The Lifespan Trajectories of Brain Activities Related to Cognitive Control

Link: https://www.biorxiv.org/content/10.1101/2023.08.20.554018v3

Cite: Li, Z., Petersen, I. T., Wang, L., Radua, J., Yang, G., & Liu, X. (2023). The Lifespan Trajectories of Brain Activities Related to Cognitive Control. bioRxiv.

[Note]:
The codes are only for the analyses with Rstudio, but not for the analyses performed within the GingerALE and SDM-PSI softwares, which can be done manually with the GUIs.

[System requirements]:
R-4.2.3
Rstudio-2023.06.0

[Required libraries] (other version not tested, but should be generally fine):
metafor-4.0-0
stringr-1.5.1
performance-0.10.3
car-3.1-2
visreg-2.7.0
MuMIn-1.47.5
mgcv-1.8-42
MASS-7.3-58.2
readxl-1.4.3
fastDummies-1.7.3

[Instructions to run on data]:
Before running the codes, please set the working directory to the ./code.

[Expected output]:
a folder named "laterality" within "meta_data" when running the codes laterality+groupcompare*.R;
a folder named "plot" within "meta_data/analysis_Adult-ChildOld" when running other R codes;
The plot_agedist.m generates Figure S2.

[Expected run time]:
Typically within 1 min for all codes except the code\GAM+modelcomp+plot_meanROI_testmeanage.R, which could take up to one hour (depending on the computer performance) with 1000 iterations.

[Instruction for use]:
Before using these code, you need to run the analyses within SDM-PSI, and then extract data from ROIs, then put the extracted data into the directory data\extracted_data_contrastanalysis. Then you need to change the directories (basedir and scriptdir) within each code.

[Reproduction instruction]: 
Simply run these codes, you can reproduce the plots and statistics reported within our manuscript
