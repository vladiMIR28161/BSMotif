#!/usr/bin/env python3

import re

def extract_min_branch_code(s):
    # Finds all numeric codes in curly brackets and returns the minimum (by hierarchy)
    matches = re.findall(r"\{(.*?)\}", str(s))
    tuples = [tuple((0, int(p)) if p.isdigit() else (1, p) for p in re.findall(r'\d+|[a-z]+', m)) for m in matches]
    return min(tuples) 
    
def extract_branch_code(s):
    # Extracts the numeric code from the Branch column
    m = re.search(r"\{(.*?)\}", str(s))
    return tuple((0, int(p)) if p.isdigit() else (1, p) for p in re.findall(r'\d+|[a-z]+', m.group(1)))
    
def renumber_branch(branch_series):
    # reindexing of ordinal numbers within each branch
    m = re.search(r"(\{.*?\})", branch_series)
    return m.group(1)
