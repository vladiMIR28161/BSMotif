#!/usr/bin/env python3

from bsmotif.io import read_classification, read_tomtom, save_results
from bsmotif.preprocessing import fill_classification, calculate_score_tf
from bsmotif.postprocessing import extract_min_branch_code, extract_branch_code, renumber_branch
from pathlib import Path
import bsmotif.hierarchical_classification as hierarchical_classification 
import pandas as pd
import sys
import subprocess
import numpy as np

def main():
    if len(sys.argv) != 4:
        print("Usage: bsmotif <classification.tsv> <pfm_dir> <output_dir>")
        sys.exit(1)
        
    
    input_classification = sys.argv[1]
    pfm_dir = sys.argv[2]
    output_tomtom_dir = sys.argv[3]
    result_score = output_tomtom_dir + "/Score.tsv"
    result_classification = output_tomtom_dir + "/Branches.tsv"
    
    # Define the path to the bash script (in the same folder as main.py)
    script_dir = Path(__file__).resolve().parent
    bash_script = script_dir / "run_tomtom.sh"

    bash_script.chmod(0o755)

    print("\n________Launching Tomtom via run_tomtom.sh________")
    result = subprocess.run(
        [str(bash_script), str(input_classification), str(pfm_dir), str(output_tomtom_dir)],
        cwd=script_dir,
        capture_output=True,
        text=True
    )

    # Checking for errors
    if result.returncode != 0:
        print("Error running run_tomtom.sh:")
        print(result.stderr)
        sys.exit(1)

    print("Tomtom completed successfully!")
    print(result.stdout)
    
    tomtom_result = output_tomtom_dir + "/Tomtom_results.tsv"

    classification = read_classification(input_classification)
    tomtom = read_tomtom(str(tomtom_result))
    print(tomtom)

    tomtom = fill_classification(tomtom, classification)
    tomtom = calculate_score_tf(tomtom)
    tomtom = save_results(tomtom, result_score)

    for superclass in tomtom.Query_superclass.unique():
        res_tomtom_filtr = tomtom.loc[tomtom.Query_superclass == superclass].reset_index(drop=True)
        class_ = hierarchical_classification.hierarchical_classification_tf (res_tomtom_filtr, 'Query_superclass', 'Target_superclass', 'Query_class', 'Target_class')
        family = hierarchical_classification.hierarchical_classification_tf (class_, 'Query_class', 'Target_class', 'Query_family', 'Target_family')
        subfamily = hierarchical_classification.hierarchical_classification_tf (family, 'Query_family', 'Target_family', 'Query_subfamily', 'Target_subfamily')
        gene = hierarchical_classification.hierarchical_classification_tf (subfamily, 'Query_subfamily', 'Target_subfamily', 'Query_gene', 'Target_gene')

    # Sorting results_df
    hierarchical_classification.results_df['Branch_key'] = hierarchical_classification.results_df['Branch'].apply(extract_min_branch_code)
    hierarchical_classification.results_df['List_key'] = hierarchical_classification.results_df['List'].apply(extract_branch_code)

    hierarchical_classification.results_df = hierarchical_classification.results_df.sort_values(
        by=['Branch_key', 'List_key'],
        key=lambda col: col
    ).reset_index(drop=True)

    # Removing old numbers
    hierarchical_classification.results_df['Branch'] = hierarchical_classification.results_df['Branch'].apply(renumber_branch)

    # Add new numbers in the order they appear within each branch.
    hierarchical_classification.results_df['Branch'] = (
        hierarchical_classification.results_df.groupby('Branch')['Branch']
          .transform(lambda s: [f"{s.iloc[0]} {i+1}" for i in range(len(s))])
    )

    # Removing auxiliary columns
    hierarchical_classification.results_df = hierarchical_classification.results_df.drop(columns=['Branch_key', 'List_key'])
    
    hierarchical_classification.results_df.to_csv(result_classification,
        sep='\t',
        index=False,
        encoding='utf-8',
        float_format='%.3f', 
        quoting=3)
    
if __name__ == "__main__":
    main()