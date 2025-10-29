#!/usr/bin/env python3

import pandas as pd
import numpy as np

def read_classification(path):
    return pd.read_csv(path, sep='\t', index_col=0)

def read_tomtom(path):
    df = pd.read_csv(path, sep='\t')
    df[['Query_superclass', 'Query_class', 'Query_family', 'Query_subfamily', 'Query_gene',
        'Target_superclass', 'Target_class', 'Target_family', 'Target_subfamily', 'Target_gene',
        'Log10pvalue', 'Score_TF']] = None
    df['pvalue'] = df['pvalue'].astype(float)
    df['Log10pvalue'] = -1 * np.log10(df['pvalue'])
    return df

def save_results(df, path):
    df.to_csv(path, sep='\t')
    return df