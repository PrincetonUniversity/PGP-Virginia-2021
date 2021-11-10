
"""

Add new plan script for VA Congressional maps
11/8/21

Attaches plans in shapefile format to a precinct dataset

"""

#%%  #### IMPORTS

import pandas as pd
import geopandas as gpd
import maup
import os
import zipfile


#%% #### STATE SETUP

'''
Set the state abbreviation, FIPS code, epsg and chamber details for new maps 

'''
state_name = 'VA'
state_fips = '51'
state_epsg = 26917
chamber = 'Congressional'
dist_chamber = 11

'''
Change working directory to folder storing data for state
'''
os.chdir('Path/to/folder/')
cwd = os.getcwd()
cwd

#%% ### FUNCTIONS

def checkGeo(gdf):
    '''
    
    Parameters
    ----------
    gdf : geoDataFrame
        geoDataFrame to check for topology errors

    Returns
    -------
    None.

    '''
    # tests for valid geometries
    if gdf['geometry'].is_valid.all():
        print('no invalid geometries')
    else:
        print('warning: invalid geometries')
        
    # tests for empty geometries
    if gdf['geometry'].is_empty.any():
        print('warning: empty geometries')
    else:
        print('no empty geometries')
        
    # tests for missing geometries
    if gdf['geometry'].isna().any():
        print('warning: missing geometries')
    else:
        print('no missing geometries')
        
    # test for duplicated indices
    if gdf.index.duplicated().any():
        print('warning: duplicated index')
    else:
        print('no duplicated indices')
        
    ## test for duplicated geometries 
    ## uncomment if needed - time intensive
    # if gdf.geometry.duplicated().any():
    #     print('warning: duplicated geometry')
    # else:
    #     print('no duplicated geometries')

def assignPlans(pl_dict):  
    '''
    
    assigns plans given a dictionary pl_dict 
    run once after each set of plans is in a pl_dict

    '''
    
    for k,v in pl_dict.items():    
        print(k)   
        
        plan_shp = v
        dist_col = 'DISTRICT'
    
        # detect number of districts
        precs.dtypes
        num_dists = len(plan_shp)
        dist_type = "None"
        
        if(num_dists == dist_chamber):
            dist_type = chamber
        else:
            dist_type = 'Other Chamber'
            
        print('District type: ' + dist_type)
        print('Number of districts: ' + str(num_dists))
        
        # check the new plan
        plan_shp = plan_shp.to_crs(proj_crs)
        plan_shp['geometry'] = plan_shp.geometry.buffer(0)
        checkGeo(plan_shp)
        
        # check for district column
        if dist_col in plan_shp.columns:
            print('District assignment column PRESENT')
        else:
            print('WARNING: District assignment column MISSING')
         
        # check for unique district columns
        if plan_shp[dist_col].nunique() == num_dists:
            print('District assignment column contains unique values')
        else:
            print('WARNING: District assignment column contains duplicate values')
            
        # assign plan with spatial join 
        print("Starting assignment for plan " + k)
        assignment = maup.assign(precs, plan_shp)
        precs[k] = assignment
        
        if not os.path.exists('./Dashboard/'):
            os.makedirs('./Dashboard/')
       
        # group by new columns - all census + election columns
        agg_cols = prec_cols + [k] 
        precs_grouped = precs[agg_cols].groupby(k).sum()
        
        precs_grouped.to_csv('./Dashboard/{0}_precs_grouped.csv'.format(k))
        
#%%  ### IMPORT PRECINCT DATA

path_to_zip_file = '{0}/{1}_precs_all_data20.zip'.format(cwd, state_name) 
    
with zipfile.ZipFile(path_to_zip_file, 'r') as zip_ref:
    zip_ref.extractall('./{0}_precs_all_data20/'.format(state_name))
    
precs = gpd.read_file('./{0}_precs_all_data20/{0}_precs_all_data20.shp'.format(state_name)) 

#%%  ### INITIAL PRECINCT DATA CLEAN

## Create unique GEOID20 column if it is missing in precinct dataset
# precs.dtypes
# precs['GEOID20'] = precs['STATEFP20'] + precs['COUNTYFP20'] + precs['VTDST20']  

# Check for unique GEOID20 values
if (precs['GEOID20'].nunique() == len(precs)):
    print('GEOID20 column contains unique values')
else:
    print('WARNING: GEOID20 column contains duplicate values')
  
# Set unique GEOID20 column to index for assignment
precs.set_index('GEOID20', inplace=True)

    
# CRS
print(precs.crs)

proj_crs = precs.crs

# buffer(0) geometry in case of invalid geometries
precs['geometry'] = precs.geometry.buffer(0)

# Check for valid geometries
checkGeo(precs)


all_columns = list(precs.columns)

#%%  ### SET COLUMNS FOR AGGREGATION

'''
Preview the columns listed in all_columns to see which variables and indices
we need to aggregate to get district estimates. Change the line below to 
correspond to the indices of these variables in the precinct data.  This should
include all columns containing census or election data.
'''

prec_cols = list(precs.columns)[21:64]


#%%  ### ASSIGN PROPOSED PLANS

'''
Reads in all folders within a "Proposed" folder and scans for shapefiles
within those folders, then attaches them to the precincts and aggregates data.


Change the directory of the "Proposed" folder where you are storing shapefiles,
and change the variable assign_col to the column that contains the district 
number for plans. 
'''

# load shapefiles from folder
folders = os.listdir('./Plans/')
plan_dict = {}

# get all proposed plans for scoring
for pl in folders:  
    print('Folder: {0}'.format(pl))
    in_folder = os.listdir('./Plans/{0}'.format(pl))   
    for cur in in_folder:
        if(cur[-4:] == '.shp'):
            basename = cur[:-4]
            filepath = '{0}/Plans/{1}/{2}'.format(cwd, pl, cur)
            plan_dict[basename] = gpd.read_file(filepath)
            print('File: {0}'.format(cur))
            print('File path: {0}\n'.format(filepath))
            
    
plans = list(plan_dict.values())

# change this column name if it differs in the shapefile
assign_col = 'DISTRICTNO'

# create int and str district columns for assignment and geojson
for pl in plan_dict.values():
    print(pl[assign_col].nunique())
    pl['DISTRICT'] = pl[assign_col].map(lambda x:int(x))
    pl['DIST_NAME'] = pl['DISTRICT'].map(lambda x:str(x).zfill(3))
    if ((assign_col != "DISTRICT") and (assign_col != "DIST_NAME")):
        pl.drop(columns=[assign_col], inplace=True)
    pl.set_index('DIST_NAME',inplace=True) 

  
# save geojson with shapes and district numbers only
for k, v in plan_dict.items():
    print(k)

    plan_gj = v[['DISTRICT', 'geometry']]
    plan_gj = plan_gj.to_crs(epsg=4326)
    
    if not os.path.exists('./Plans/GeoJSON'):
       os.makedirs('./Plans/GeoJSON')
        
    plan_gj.to_file('./Plans/GeoJSON/{0}_{1}.geojson'.format(state_name, k), driver='GeoJSON')

# assign proposed plans
assignPlans(plan_dict)



#%% ### RENAME COLUMNS 

'''
Rename columns from the long file name to an appropriate and unique
column name with 10 characters or fewer to use as the plan assignment
in the precinct shapefile
'''

cols  = list(precs.columns)

# rename columns for new plans
precs.rename(columns={'OriginalColumnName':'NewColumnName'}, inplace=True)

#%% ### FIND UNASSIGNED PRECINCTS

'''

Determines if there are any issues with unassigned precincts
Unassigned precincts will get a null value for their assignment, which can 
throw off the population equivalence because there are now X+1 districts.

'''
cols  = list(precs.columns)

# set the index of the column(s) corresponding to the newly added plan
plancols = list(precs.columns)[68:]  


# find number of precincts unassigned in current plan
# find number of populated precincts unassigned in current plan

for pl in plancols:
    print('Plan: ' + str(pl))
    print('Number of plan districts: ' + str(precs[pl].nunique()))
    unassigned_precs = precs.loc[precs[pl].isnull()]
    unassigned_pop_precs = unassigned_precs.loc[unassigned_precs['TOT'] >0]
    print('Number of unassigned precincts: ' + str(len(unassigned_precs)))
    print('Number of unassigned and populated precincts: ' + str(len(unassigned_pop_precs))+ '\n')
    if (len(unassigned_precs) > 0):
        unassigned_precs.to_file('./Scratch/{0}_unassigned_precs.shp'.format(pl))
    if (len(unassigned_pop_precs) > 0):
        unassigned_pop_precs.to_file('./Scratch/{0}_unassigned_pop_precs.shp'.format(pl))



#%% ### SAVE UPDATED PRECINCTS FILE

precs.to_file('./Shapefiles/{0}_precs_cong_scoring20.shp'.format(state_name))  # CONGRESSIONAL


#%%



