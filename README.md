# VA Scoring Package

Collection of files and scripts to score and grade new districting maps in VA

## Contents
- PGP grading system for Virginia metrics November 2021.pdf: Scoring system description for VA specific maps.

- Shapefile: shapefile of VA precincts with attached census and election data
  - VA_precs_all_data20.shp

- Aggregation Scripts: folder for all scripts that can attach a plan to the shapefile
  - VA_add_new_plan_cong.py: script for adding VA congressional plans
  - VA_add_new_plan_hod.py: script for adding VA House of Delegates plans
  - VA_add_new_plan_ss.py: script for adding VA state senate plans

- Plans: folder for all plans you want to score
  - archive: a folder for plans that have already been added

- Scratch: empty folder for scratch work in VA_add_new_plan scripts

- Shapefiles: folder for all input data files
  - VA_precs_cong_scoring20.shp: all input data for Congressional scoring
  - VA_precs_house_scoring20.shp: all input data for House of Delegates scoring
  - VA_precs_sen_scoring20.shp: all input data for State Senate scoring

- Scoring Scripts: folder for all scripts that score plans on a number of metrics
  - VA_congressional_scoring.py: scores Congressional plans
  - VA_state_senate_scoring.py: scores State Senate plans
  - VA_house_of_delegates_scoring.py: scores House of Delegates plans
- Benchmark Thresholds for Grading


## Workflow

When you have a plan, or set of plans, in the form of shapefiles that you would like to score/grade,
use the chamber-specific aggregation script and the provided shapefile to attach the plan.

Make sure to only plans for one chamber are in the /Plans folder at once.
Follow the instructions in VA_add_new_plan script and run through it cell by cell.
You will want to run through this once per chamber, and may have to do it every
time a new plan is out you want to score. The /archive folder is useful for keeping old maps.

When you finish the VA_add_new_plan.py script, it will output a file called
VA_precs_chamber_scoring20.shp where the chamber will be cong, senate, or house.
This is an input for scoring scripts. Once you have
the VA_precs_chamber_scoring20.shp file, open up the scoring script that corresponds to the chamber.
Fill in the information in the "# ----- CODE TO CHANGE FOR EACH RUN OF CODE" section.
The district_column_names refers to the column names of your district plans in the shapefile.
Make sure you have a folder on your local computer where the resulting json will be saved in the folder_to_store variable.
Once you have filled out this section, run the file and it
should output jsons for each district plan that you are scoring into your
folder_to_store. You can open up the jsons to see a dictionary of data that you have saved for each plan.

With the json data, you can now compare the metric results to the benchmark thresholds to determine the grades
that the plan would receive for competitiveness, partisan fairness, geography as well as the final grade.


## Metadata

- 'GEOID20': Voting district identifier
- 'NAMELSAD20': 2020 Census name and the translated legal/statistical
area description for voting district
- 'NAME20': 2020 Census voting district name
- 'VTDI20': 2020 Census voting district indicator
- 'VTDST20': 2020 Census voting district code
- 'COUNTYFP20': 2020 Census county FIPS code
- 'STATEFP20': 2020 Census state FIPS code
- 'LSAD20': 2020 Census legal/statistical area description code for
voting districts
- 'MTFCC20': MAF/TIGER feature class code (G5240)
- 'FUNCSTAT20': 2020 Census functional status
- 'ALAND20': 2020 Census land area
- 'AWATER20': 2020 Census water area
- 'INTPTLAT20': 2020 Census latitude of the internal point
- 'INTPTLON20': 2020 Census longitude of the internal point
- 'CD116': 116th congressional district FIPS code
- 'SLDU18': Current state legislative district upper chamber code
- 'SLDL18': Current state legislative district lower chamber code
- 'C12Dist': 2012 Enacted Congressional district number
- 'C16Dist': 2016 Enacted Congressional district number
- 'S11Dist': 2011 Enacted State Senate district number
- 'H11Dist': 2011 Enacted House of Delegates district number
- 'H19Dist': 2019 Enacted House of Delegates district number
- 'TOT': total population
- 'WH_A': white alone
- 'BL_A': black or african american alone
- 'NA_A': american indian and alaska native alone
- 'AS_A': asian alone
- 'PC_A': native hawaiian and other pacific islander alone
- 'SO_A': some other race alone
- 'TMO': population of two or more races
- 'TOT_H': total hispanic population
- 'TOT_NH': total nonhispanic population
- 'NHWH_A': nonhispanic white alone
- 'NHBL_A': nonhispanic black alone
- 'NHNA_A': nonhispanic american indian and alaska native alone
- 'NHAS_A': nonhispanic asian alone
- 'NHPC_A': nonhispanic native hawaiian and other pacific islander alone
- 'NHSO_A': nonhispanic some other race alone
- 'NHTMO': nonhispanic population of two or more races
- 'TOT_VAP': total voting age population
- 'WH_A_VAP': white alone voting age population
- 'BL_A_VAP': black or african american alone voting age population
- 'NA_A_VAP': american indian and alaska native alone voting age population
- 'AS_A_VAP': asian alone voting age population
- 'PC_A_VAP': native hawaiian and other pacific islander alone voting age population
- 'SO_A_VAP': some other race alone voting age population
- 'TMO_VAP': population of two or more races voting age population
- 'TOT_HVAP': total hispanic voting age population
- 'TOT_NHVAP': total nonhispanic voting age population
- 'NHWH_A_VAP': nonhispanic white alone voting age population
- 'NHBL_A_VAP': nonhispanic black or african american alone voting age population
- 'NHNA_A_VAP': nonhispanic american indian and alaska native alone voting age population
- 'NHAS_A_VAP': nonhispanic asian alone voting age population
- 'NHPC_A_VAP': nonhispanic native hawaiian and other pacific islander alonevoting age population
- 'NHSO_A_VAP': nonhispanic some other race alone voting age population
- 'NHTMO_VAP': nonhispanic population of two or more races voting age population
- 'G16DPRS': Number of votes for 2016 Democratic presidential candidate
- 'G16RPRS': Number of votes for 2016 Republican presidential candidate
- 'G16OPRS': Number of votes for 2016 other presidential candidates
- 'G17DGOV': Number of votes for 2017 Democratic gubernatorial candidate
- 'G17RGOV': Number of votes for 2017 Republican gubernatorial candidate
- 'G17OGOV': Number of votes for 2017 other gubernatorial candidates
- 'G18DSEN': Number of votes for 2018 Democratic senate candidate
- 'G18RSEN': Number of votes for 2018 Republican senate candidate
- 'G18OSEN': Number of votes for 2018 other senate candidates


## Sources

Precinct-level election data was downloaded from [MGGG](https://github.com/mggg-states/VA-shapefiles).  Population and demographic comes from the U.S. Census Bureau, and does not include prison adjusted population.
