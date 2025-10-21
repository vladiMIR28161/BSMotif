#!/usr/bin/env python3

from bsmotifs.io import read_classification, read_tomtom, save_results
from bsmotifs.preprocessing import fill_classification, calculate_score_tf
import bsmotifs.hierarchical_classification as hierarchical_classification 
import pandas as pd
import sys

def main():
    if len(sys.argv) != 5:
        print("Usage: bsmotifs <input_classification.tsv> <input_tomtom.tsv> <score.tsv> <classification.xlsx>")
        sys.exit(1)
        
    
    input_classification = sys.argv[1]
    input_tomtom = sys.argv[2]
    result_score = sys.argv[3]
    result_classification = sys.argv[4]

    classification = read_classification(input_classification)
    tomtom = read_tomtom(input_tomtom)

    tomtom = fill_classification(tomtom, classification)
    tomtom = calculate_score_tf(tomtom)
    tomtom = save_results(tomtom, result_score)

    for superclass in tomtom.Query_superclass.unique():
        res_tomtom_filtr = tomtom.loc[tomtom.Query_superclass == superclass].reset_index(drop=True)
        class_ = hierarchical_classification_tf (res_tomtom_filtr, 'Query_superclass', 'Target_superclass', 'Query_class', 'Target_class')
        family = hierarchical_classification_tf (class_, 'Query_class', 'Target_class', 'Query_family', 'Target_family')
        subfamily = hierarchical_classification_tf (family, 'Query_family', 'Target_family', 'Query_subfamily', 'Target_subfamily')
        gene = hierarchical_classification_tf (subfamily, 'Query_subfamily', 'Target_subfamily', 'Query_gene', 'Target_gene')

    hierarchical_classification.results_df.to_excel(result_classification)
    
if __name__ == "__main__":
    main()