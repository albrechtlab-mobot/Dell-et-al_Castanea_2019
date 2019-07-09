# Dell-et-al_Castanea_2019
A repository that contains data and R scripts for running analyses performed in the paper "Germination traits in the threatened southeastern grassland endemic, Marshallia mohrii (Asteraceae)" written by Noah Dell, Quinn Long, and Matthew Albrecht and published in Castanea (2019). 

The files in this repository consists of 5 data sets and 2 R scripts. The files are described as follows:

  Marmoh_Analyses.r - An R script which performs time to event analysis and GLMs for the following 4 data sets as described in Dell 
  et al (2019). It uses the following data:
    
   Marmoh_Analysis 1_T2E.csv - A dataset containing cumulative germination totals for each replicate in each week of observation for
   Experiment 1. The spreadsheet contains the following columns:
    
      Species	- Species Name
      SpAbrev - Species Abbreviation
      Accession #	- Accession number for the seed collection from which seeds were sourced for the experiment
      Collection Location - Location at which seeds were collected	
      MaternalLine - If seeds were collected from plants individually, this is a unique identifier, otherwise it is labelled as 
        "NA - bulk"
      Treatment	- Column that indicates the cold stratification temperature (if any) and the incubation temperature for that replicate
      TreatNotes - Notes on how plants were watered and treated during the experiment
      Replicate# - Unique identifier for each dish within each treatment
      dLabel - Unique identifier for each dish among all treatments
      Seeds/Dish - Describes the number of seeds in each dish at the start of the experiment
      Substrate	- Lists the substrate on which seeds were placed
      PreTreatment - NA, notes on how seeds were treated before the experiment
      StartTemp	- Temperature at start of experiment
      StratifLength - Length of cold stratification treatment	
      StartDate	- Date on which seeds were placed on petri dish and treatment began
      StartPhotoperiod - Light regime at start of experiment
      Temp2	- Temperature after cold stratification (if applicable)
      Temp2Date	- Date on which seeds were moved from cold stratification to incubation (if applicable)
      Temp2Photo	- Light regime during incubation (if applicable)
      IncubationTemp	- Temperature during incubation
      Start	- First week of an observation period
      End	- Second week of an observation period
      Germ	- number of germinants during the current observation period
      Unvia	- number of unviable seeds during the current observation period
      tGrm	- number of germinants from the start of the experiment to the current observation period
      tUnv - number of unviable seeds from the start of the experiment to the current observation period
      dChk - Seeds remaining on the dish at the current observation period. 
      
   Marhmoh_Analysis 2_cumulative.csv - A dataset containing cumulative germination totals for each replicate across all four weeks for
   Experiment 2. Columns are the same as in Experiment 1, except there are no "Start" and "End" columns, and two columns at the end
   (PropGerm and PropInv) denote the proportion of seeds that germinated and the proportion of seeds that became inviable, 
   respectively.
     
   Marmoh_Analsysis 3_cumulative.csv - A dataset containing cumulative germination totals for each replicate across all four weeks for
   Experiment 3. Columns are the same as in Experiment 2. 
     
   Marmoh_Analsysis 4_cumulative.csv - A dataset containing cumulative germination totals for each replicate across all four weeks for
   Experiment 4. Columns are the same as in Experiment 2. 
     
    
  Marmoh_SoilTemp.r - A script to calculate daily and monthly average, high, and low temperatures. It uses the following data:
    
   Marmoh_SoilTemp.csv - A dataset containing soil and airtemperatures collected every 90 minutes at the Ketona Dolomite Glades. 
   Soil temperature was recorded using two iButtons placed 1 m apart. The spreadsheet contains the following columns:
   
      Site - The site at which the temperatures were collected
      Position - The position of the iButton (upslope or downslope)
      Medium - States whether the sensor was recording air or soil temperature
      Date/Time - Date and time at which each sample was taken
      Unit - Unit of measurement (degrees celsius)
      Temp - Recorded temperature

To run either of these R scripts, you will need to save the appropriate files in a folder and adjust the file paths in the command to 
set the working directory at the beginning of each script.  
