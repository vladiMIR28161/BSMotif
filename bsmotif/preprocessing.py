#!/usr/bin/env python3

import pandas as pd

def fill_classification(tomtom, classification):
    for i in pd.concat([tomtom['Query_ID'], tomtom['Target_ID']]).unique():
        superclass = classification.loc[classification['ID'] == i].Superclass.unique()[0]
        usually_class = classification.loc[classification['ID'] == i].Class.unique()[0]
        family = classification.loc[classification['ID'] == i].Family.unique()[0]
        subfamily = classification.loc[classification['ID'] == i].Subfamily.unique()[0]
        gene = classification.loc[classification['ID'] == i].Gene.unique()[0]
        tomtom.loc[tomtom['Query_ID'] == i, ['Query_superclass', 'Query_class', 'Query_family', 'Query_subfamily', 'Query_gene']] = [superclass, usually_class, family, subfamily, gene]
        tomtom.loc[tomtom['Target_ID'] == i, ['Target_superclass', 'Target_class', 'Target_family', 'Target_subfamily', 'Target_gene']] = [superclass, usually_class, family, subfamily, gene]
    return tomtom

def calculate_score_tf(tomtom):
    for i in range(len(tomtom.Query_gene)):
        if pd.isna(tomtom.Score_TF[i]):
            tf_1 = tomtom.Query_gene[i]
            tf_2 = tomtom.Target_gene[i]
            print (f"Calculate score TF {tf_1} {tf_2}")
            local_tomtom = tomtom.loc[(((tomtom.Query_gene == tf_1) & (tomtom.Target_gene == tf_2)) |
                                       ((tomtom.Query_gene == tf_2) & (tomtom.Target_gene == tf_1)))]
            tomtom.loc[(((tomtom.Query_gene == tf_1) & (tomtom.Target_gene == tf_2)) | 
                        ((tomtom.Query_gene == tf_2) & (tomtom.Target_gene == tf_1))), 'Score_TF'] = max(local_tomtom.Log10pvalue)
        else:
            continue
    return tomtom