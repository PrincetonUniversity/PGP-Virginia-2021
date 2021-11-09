# ---- Imports -------
import csv
import pickle
import os
from functools import partial
import json
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
import time
from itertools import groupby
from statistics import mean
import math, random
from scipy.stats import norm
import pandas as pd
plt.style.use('seaborn-whitegrid')
from enum import Enum
import gerrymetrics as gm
import sys
import argparse
from locality_splitting import metrics


class CountySplit(Enum):
    NOT_SPLIT = 0
    NEW_SPLIT = 1
    OLD_SPLIT = 2



# ----- CODE TO CHANGE FOR EACH RUN OF THE CODE

# THIS IS THE PATH TO THE SHAPEFILE WITH THE ATTACHED PLAN
# FOR EXAMPLE IT MIGHT LOOK LIKE:
# path_to_shapefile = "/Users/ari/Documents/PGP/Dashboard/Virginia/VA_precs_sen_scoring20/VA_precs_sen_scoring20.shp"
path_to_shapefile = 

# THIS IS THE NAMES OF THE DISTRICT ASSIGNMENT COLUMNS THAT YOU WANT TO SCORE
# IT HAS TO BE IN AN ARRAY EVEN IF IT'S JUST ONE COLUMN
# FOR EXAMPLE IT MIGHT LOOK LIKE:
# district_column_names = ["S21_1Dist"]
# or
# district_column_names = ["S21_1Dist", "S21_2Dist"]
district_column_names = 

# THIS IS THE PATH TO THE FOLDER THAT YOU WANT TO STORE THE RESULTING JSON
# FOR EXAMPLE IT MIGHT LOOK LIKE:
# folder_to_store = "/Users/ari/Documents/PGP/Dashboard/Virginia/"
folder_to_store = 
os.chdir(folder_to_store)

# ---- DON'T CHANGE ANYTHING BELOW THIS LINE

# --- Internal Functions & Definitions ---
state = "VA"
num_dist = 40
county_assignment_col = "COUNTYFP20"
cong = False
sen = True
house = False
election_names = [
    "PRES16",
    "GOV17",
    "USSEN18"
]
election_columns = [
    ["G16DPRS", "G16RPRS"],
    ["G17DGOV", "G17RGOV"],
    ["G18DSEN", "G18RSEN"]
]
num_elections = len(election_names)
pop_col = "TOT"
WVAP_col = "NHWH_A_VAP"
BVAP_col = "BL_A_VAP"
HVAP_col = "TOT_HVAP"
AVAP_col = "AS_A_VAP"
Native_col = "NA_A_VAP"
Pacific_col = "PC_A_VAP"
VAP_col = "TOT_VAP"
cong_seat_pop = 761169
cong_factor = 1.01
leg_factor = 1.05
comp_lower = 0.465
comp_higher = 0.535


# ---- Smallest enclosing circle -----
# Smallest enclosing circle - Library (Python)
# 
# Copyright (c) 2020 Project Nayuki
# https://www.nayuki.io/page/smallest-enclosing-circle
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public License
# along with this program (see COPYING.txt and COPYING.LESSER.txt).
# If not, see <http://www.gnu.org/licenses/>.
# 

import math, random


# Data conventions: A point is a pair of floats (x, y). A circle is a triple of floats (center x, center y, radius).

# Returns the smallest circle that encloses all the given points. Runs in expected O(n) time, randomized.
# Input: A sequence of pairs of floats or ints, e.g. [(0,5), (3.1,-2.7)].
# Output: A triple of floats representing a circle.
# Note: If 0 points are given, None is returned. If 1 point is given, a circle of radius 0 is returned.
# 
# Initially: No boundary points known
def make_circle(points):
	# Convert to float and randomize order
	shuffled = [(float(x), float(y)) for (x, y) in points]
	random.shuffle(shuffled)
	
	# Progressively add points to circle or recompute circle
	c = None
	for (i, p) in enumerate(shuffled):
		if c is None or not is_in_circle(c, p):
			c = _make_circle_one_point(shuffled[ : i + 1], p)
	return c


# One boundary point known
def _make_circle_one_point(points, p):
	c = (p[0], p[1], 0.0)
	for (i, q) in enumerate(points):
		if not is_in_circle(c, q):
			if c[2] == 0.0:
				c = make_diameter(p, q)
			else:
				c = _make_circle_two_points(points[ : i + 1], p, q)
	return c


# Two boundary points known
def _make_circle_two_points(points, p, q):
	circ = make_diameter(p, q)
	left  = None
	right = None
	px, py = p
	qx, qy = q
	
	# For each point not in the two-point circle
	for r in points:
		if is_in_circle(circ, r):
			continue
		
		# Form a circumcircle and classify it on left or right side
		cross = _cross_product(px, py, qx, qy, r[0], r[1])
		c = make_circumcircle(p, q, r)
		if c is None:
			continue
		elif cross > 0.0 and (left is None or _cross_product(px, py, qx, qy, c[0], c[1]) > _cross_product(px, py, qx, qy, left[0], left[1])):
			left = c
		elif cross < 0.0 and (right is None or _cross_product(px, py, qx, qy, c[0], c[1]) < _cross_product(px, py, qx, qy, right[0], right[1])):
			right = c
	
	# Select which circle to return
	if left is None and right is None:
		return circ
	elif left is None:
		return right
	elif right is None:
		return left
	else:
		return left if (left[2] <= right[2]) else right


def make_diameter(a, b):
	cx = (a[0] + b[0]) / 2
	cy = (a[1] + b[1]) / 2
	r0 = math.hypot(cx - a[0], cy - a[1])
	r1 = math.hypot(cx - b[0], cy - b[1])
	return (cx, cy, max(r0, r1))


def make_circumcircle(a, b, c):
	# Mathematical algorithm from Wikipedia: Circumscribed circle
	ox = (min(a[0], b[0], c[0]) + max(a[0], b[0], c[0])) / 2
	oy = (min(a[1], b[1], c[1]) + max(a[1], b[1], c[1])) / 2
	ax = a[0] - ox;  ay = a[1] - oy
	bx = b[0] - ox;  by = b[1] - oy
	cx = c[0] - ox;  cy = c[1] - oy
	d = (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by)) * 2.0
	if d == 0.0:
		return None
	x = ox + ((ax*ax + ay*ay) * (by - cy) + (bx*bx + by*by) * (cy - ay) + (cx*cx + cy*cy) * (ay - by)) / d
	y = oy + ((ax*ax + ay*ay) * (cx - bx) + (bx*bx + by*by) * (ax - cx) + (cx*cx + cy*cy) * (bx - ax)) / d
	ra = math.hypot(x - a[0], y - a[1])
	rb = math.hypot(x - b[0], y - b[1])
	rc = math.hypot(x - c[0], y - c[1])
	return (x, y, max(ra, rb, rc))


_MULTIPLICATIVE_EPSILON = 1 + 1e-14

def is_in_circle(c, p):
	return c is not None and math.hypot(p[0] - c[0], p[1] - c[1]) <= c[2] * _MULTIPLICATIVE_EPSILON


# Returns twice the signed area of the triangle defined by (x0, y0), (x1, y1), (x2, y2).
def _cross_product(x0, y0, x1, y1, x2, y2):
	return (x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0)

#----- Other Internal functions -----

def get_circle(geometry):
    all_points = []
    if geometry.geom_type == "MultiPolygon":
        for polygon in geometry:
            x, y = polygon.exterior.coords.xy
            for i in range(len(x)):
                all_points.append((x[i], y[i]))
    else:
        x, y = geometry.exterior.coords.xy
        for i in range(len(x)):
            all_points.append((x[i], y[i]))
    circum_circle = make_circle(all_points)
    radius = circum_circle[2]
    return (math.pi * radius ** 2)


def subset_list(percentage_list, lower_threshold, upper_threshold):
    percentage_intermediate = percentage_list[percentage_list >= lower_threshold]
    return len(percentage_intermediate[percentage_intermediate < upper_threshold])

def votes_to_seats(voteshare, st_dev):
    normalized_voteshare = (voteshare - 0.5) / st_dev
    return norm.cdf(normalized_voteshare)

def cube_law(voteshare):
    cube = (voteshare / (1 - voteshare))**3
    return cube / (1 + cube)



# ---- Load in graph ------
data = gpd.read_file(path_to_shapefile)

if cong:
    chamber = "congressional"
elif sen:
    chamber = "state-senate"
elif house:
    chamber = "state-house"



for district_assignment_col in district_column_names:

    agg_data = data.dissolve(by=district_assignment_col, aggfunc='sum')
    
    # ---- Geographic Metrics ------
    
    # County Splits
    if cong:
        threshold_pop = cong_factor * cong_seat_pop
    else:
        total_pop = sum(data[pop_col])
        threshold_pop = leg_factor * total_pop / num_dist
    
    num_splits = 0
    lower_county_splits = 0
    counties = list(filter(None, pd.unique(data[county_assignment_col])))
    
    for county in counties:
        county_subset = data[data[county_assignment_col] == county]
        split_times = len(pd.unique(county_subset[district_assignment_col]))
        if split_times > 1:
            num_splits = num_splits + 1
        county_pop = sum(county_subset[pop_col])
        if county_pop > threshold_pop:
            lower_county_splits = lower_county_splits + 1
        
    upper_county_splits = len(counties)
    
    # Split pairs
    splits_measures = metrics.calculate_all_metrics(data, plan_col = district_assignment_col, lclty_col=county_assignment_col, pop_col=pop_col)
    split_pairs_measure = splits_measures["split_pairs"]
    
    
        
    # --- Compactness -----
    
    # Polsby-Popper & Reock
    district_labels = agg_data.index
    for district in district_labels:
        area = agg_data.at[district, "geometry"].area
        per = agg_data.at[district, "geometry"].length
        pp = 4 * math.pi * area / per ** 2
        agg_data.at[district, "Polsby-Popper"] = pp
        area_of_bounding_circle = get_circle(agg_data.at[district, "geometry"])
        reock = area / area_of_bounding_circle
        agg_data.at[district, "Reock"] = reock
        
    min_pp = min(agg_data["Polsby-Popper"])
    avg_pp = np.average(agg_data["Polsby-Popper"])
    
    min_reock = min(agg_data["Reock"])
    avg_reock = np.average(agg_data["Reock"])
    
    
    
    # ---- Racial Fairness Metrics (geopd) -----    
    bvap_estimates = agg_data[BVAP_col] / agg_data[VAP_col]
    bvap_estimates_dict = bvap_estimates.to_dict()
    
    hvap_estimates = agg_data[HVAP_col] / agg_data[VAP_col]
    hvap_estimates_dict = hvap_estimates.to_dict()
    
    avap_estimates = agg_data[AVAP_col] / agg_data[VAP_col]
    avap_estimates_dict = avap_estimates.to_dict()
    
    nvap_estimates = agg_data[Native_col] / agg_data[VAP_col]
    nvap_estimates_dict = nvap_estimates.to_dict()
    
    pvap_estimates = agg_data[Pacific_col] / agg_data[VAP_col]
    pvap_estimates_dict = pvap_estimates.to_dict()
    
    agg_data["MVAP"] = agg_data[VAP_col] - agg_data[WVAP_col]
    mvap_estimates = agg_data["MVAP"] / agg_data[VAP_col]
    mvap_estimates_dict = mvap_estimates.to_dict()
    
    vaps_list = [mvap_estimates, bvap_estimates, hvap_estimates, avap_estimates, nvap_estimates, pvap_estimates]
    
    racial_vaps = np.zeros((7,6))
    
    for i in range(6):
        racial_vaps[0, i] = subset_list(vaps_list[i], 0.30, 0.35)
        racial_vaps[1, i] = subset_list(vaps_list[i], 0.35, 0.40)
        racial_vaps[2, i] = subset_list(vaps_list[i], 0.40, 0.45)
        racial_vaps[3, i] = subset_list(vaps_list[i], 0.45, 0.50)
        racial_vaps[4, i] = subset_list(vaps_list[i], 0.50, 0.55)
        racial_vaps[5, i] = subset_list(vaps_list[i], 0.55, 0.60)
        racial_vaps[6, i] = subset_list(vaps_list[i], 0.60, 1.00)
    
    # --- Partisanship -----
    total_dem_votes = []
    total_rep_votes = []
    for i in range(len(election_names)):
        election = election_names[i]
        election_cols = election_columns[i]
        agg_data[election] = agg_data[election_cols[0]] / (agg_data[election_cols[0]] + agg_data[election_cols[1]])
        total_dem_votes.append(sum(agg_data[election_cols[0]]))
        total_rep_votes.append(sum(agg_data[election_cols[1]]))
    
    avg_election_voteshares = pd.DataFrame.mean(agg_data[election_names], axis = 1)
    
    avg_election_dict = avg_election_voteshares.to_dict()
    
    avg_election_votes = np.array(avg_election_voteshares)
    
    total_dem_votes = np.array(total_dem_votes)
    total_rep_votes = np.array(total_rep_votes)
    dem_voteshares_elections = total_dem_votes / (total_dem_votes + total_rep_votes)
    dem_voteshare = np.average(dem_voteshares_elections, axis = 0)
    
    
    # ---- Partisan Fairness Metrics ----
    dem_voteshares = avg_election_votes[avg_election_votes > 0.5]
    rep_voteshares = avg_election_votes[avg_election_votes < 0.5]
    if len(dem_voteshares) == 0:
        avg_dem_win = 0
    else:
        avg_dem_win = np.average(dem_voteshares)
    if len(rep_voteshares) == 0:
        avg_rep_win = 0
    else:
        avg_rep_win = 1 - np.average(rep_voteshares)
    
    # EG assumes equal turnout -- gets a different result than gerrychain
    EG = gm.EG(avg_election_votes)
    
    partisan_bias = gm.partisan_bias(avg_election_votes)
    
    mean_median = gm.mean_median(avg_election_votes)
    
    avgactWins = subset_list(avg_election_votes, 0.5, 1)
    # ---- Competitiveness -----
    competitive_vote_shares = avg_election_dict
    avgactcomp = subset_list(avg_election_votes, comp_lower, comp_higher)
    
    
    # ---- Export to JSON ------
    
        
    plan = {
          "state": state,
          "planName": district_assignment_col,
          "numDists": num_dist,
          "countySplits": num_splits,
          "lowerBoundCountySplits": lower_county_splits,
          "upperBoundCountySplits": upper_county_splits,
          "splitPairs": split_pairs_measure,
          "minPolsby-Popper": min_pp,
          "avgPolsby-Popper": avg_pp,
          "minReock": min_reock,
          "avgReock": avg_reock,
          "ElectionResults": competitive_vote_shares,
          "ActualCompetitive": avgactcomp,
          "BVAPDistricts": bvap_estimates_dict,
          "HVAPDistricts": hvap_estimates_dict,
          "AVAPDistricts": avap_estimates_dict,
          "NVAPDistricts": nvap_estimates_dict,
          "PVAPDistricts": pvap_estimates_dict,
          "MVAPDistricts": mvap_estimates_dict,
          "avgDemWin": avg_dem_win,
          "avgRepWin": avg_rep_win,
          "EG": EG,
          "partisanBias": partisan_bias,
          "meanMedian": mean_median,
          "DemWins": avgactWins
            }
    
    plan_json = {
        "plan": plan
        }
    
    
    with open("{}-{}.json".format(state, district_assignment_col), "w") as data_file:
       json.dump(plan_json, data_file, indent=4)
   

        




