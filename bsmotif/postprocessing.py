#!/usr/bin/env python3

import re

def extract_min_branch_code(s):
    # Finds all numeric codes in curly brackets and returns the minimum (by hierarchy)
    matches = re.findall(r"\{([\d\.]+)\}", str(s))
    tuples = [tuple(map(int, m.split('.'))) for m in matches]
    return min(tuples) 
    
def extract_branch_code(s):
    # Extracts the numeric code from the Branch column
    m = re.search(r"\{([\d\.]+)\}", str(s))
    return tuple(map(int, m.group(1).split('.')))
    
def renumber_branch(branch_series):
    # reindexing of ordinal numbers within each branch
    m = re.search(r"(\{[\d\.]+\})", branch_series)
    branch_prefix = m.group(1)
    return branch_prefix 