#!/usr/bin/env python3

import re

def extract_min_branch_code(s):
    matches = re.findall(r"\{(.*?)\}", str(s))
    tuples = [tuple((0, int(p)) if p.isdigit() else (1, p) for p in re.findall(r'\d+|[a-z]+', m)) for m in matches]
    return min(tuples)

def extract_branch_code(s):
    m = re.search(r"\{(.*?)\}", str(s))
    return tuple((0, int(p)) if p.isdigit() else (1, p) for p in re.findall(r'\d+|[a-z]+', m.group(1)))

def renumber_branch(branch_series):
    m = re.search(r"(\{.*?\})", branch_series)
    return m.group(1)