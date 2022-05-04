
### Structure

This repository is organized as follows:

#### measurement_data  
Contains 6 sofa files of HRTFs of a Kemar dummy head wearing each of the tested conditions.  

#### listening_test  
Contains
  - Results: raw results for colouration and localization tests  
  - Analysis: scripts for the analysis of the listening test results  
  - Stimuli for colouration and localization tests  
  
#### model  
Contains:
  - Colouration: CLL-based model to predict colouration from the measurement_data (sofa files)  
  - Localization: Evaluation of the models may2011 and baumgartner2014 to predict the localization data.  
                  These models are available in the Auditory Modelling Toolbox (https://amtoolbox.org/).  

#### dependencies  
Contains some auxiliary functions needed to run the main scripts.  


### Usage  

If you want to run the models, first run startup.m to add the relevant 
folders to your matlab path  

For the coloration model go to model/coloration and simply run colorationModel.m  

### Dependencies  
Apart from the included files, the auditory modelling toolbox (https://amtoolbox.org/) is required 
for running the localization models 



