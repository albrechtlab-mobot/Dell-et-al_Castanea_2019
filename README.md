# Dell-et-al_Castanea_2019
A repository that contains data and R scripts for running analyses performed in the paper "Germination traits in the threatened southeastern grassland endemic, Marshallia mohrii (Asteraceae)" written by Noah Dell, Quinn Long, and Matthew Albrecht and published in Castanea (2019). 

The files in this repository consists of 5 data sets and 2 R scripts. The files are described as follows:

  Marmoh_Analyses.r - An R script which performs time to event analysis and GLMs for the following 4 data sets as described in Dell et al (2019)
    
    Marmoh_Analysis 1_T2E.csv - A dataset containing cumulative germination totals for each replicate in each week of observation for Experiment 1.
      Species	- Species Name
      SpAbrev - Species Abbreviation
      Accession #	- Accession number for the seed collection from which seeds were sourced for the experiment
      Collection Location - Location at which seeds were collected	
      MaternalLine - If seeds were collected from plants individually, this is a unique identifier, otherwise it is labelled as "NA - bulk"
      Treatment	- Column that indicates the cold stratification temperature (if any) and the incubation temperature for that replicate
      TreatNotes - Notes on how plants were watered and treated during the experiment
      Replicate# - Unique identifier for each dish within each treatment
      dLabel - Unique identifier for each dish among all treatments
      Seeds/Dish - Describes the number of seeds in each dish at the start of the experiment
      Substrate	- Lists the substrate on which seeds were placed
      PreTreatment - NA, notes on how seeds were treated before the experiment
      StartTemp	- Stratification temperature. 5 degrees celsius, unless no cold stratification was used, then incubation temp is listed
      StratifLength	
      StartDate	
      StartPhotoperiod	
      Temp2	
      Temp2Date	
      Temp2Photo	
      IncubationTemp	
      Start	
      End	
      Germ	
      Unvia	
      tGrm	
      tUnv	
      dChk
