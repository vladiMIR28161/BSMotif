### BSMotifs

BSMotifs is a Python package and pipeline for analyzing the similarity of transcription factor (TF) binding site (TFBS) motifs. It integrates the results of motif comparison with the TF classification by DNA-binding domain structure and constructs the tree of similar motifs consisting of different hierarchical levels from class to TF gene.

## Overview

The pipeline takes two input files:

1. Tomtom.tsv — results of pairwise TFBS motif comparisons between TFs (e.g., using Tomtom tool ([Gupta et al., 2007](https://doi.org/10.1186/gb-2007-8-2-r24))).

2. Classification.tsv — TF classification by DNA-binding domains (including superclass, class, family, subfamily, gene).

Outputs:

1. Score.tsv — table with motif similarity scores for each TF pair.

2. Classification.xlsx — hierarchical list of branches of similar TFBS motifs within each classification level (class, family, subfamily, TF).

Example:

* TFs from families XBP1-related {1.1.5} and CREB-related {1.1.7} are significantly similar → the pipeline merges them into one branch within the family level.

This allows identifying clusters of similar motifs within structural levels of DNA-binding domains.

## Installation

``` bash
git clone https://github.com/vladiMIR28161/BSMotifs.git
cd bsmotifs
pip install .
```

After installation, the bsmotifs command becomes available in the terminal.

## Usage

``` bash
bsmotifs input_classification.tsv input_tomtom.tsv score.tsv classification.xlsx
```

Where:

* input_classification.tsv — input table with TF classification

* input_tomtom.tsv — table with Tomtom results

* score.tsv — output file with similarity scores

* classification.xlsx — output file with hierarchical branches

## Example code

``` python
from bsmotifs.io import read_classification, read_tomtom, save_results
from bsmotifs.preprocessing import fill_classification, calculate_score_tf
from bsmotifs.hierarchical_classification import hierarchical_classification_tf
import pandas as pd

classification = read_classification("Classification.tsv")
tomtom = read_tomtom("Tomtom.tsv")

tomtom = fill_classification(tomtom, classification)
tomtom = calculate_score_tf(tomtom)
tomtom = save_results(tomtom, "Score.tsv")

# Hierarchical classification by levels
for superclass in tomtom.Query_superclass.unique():
    res_tomtom_filtr = tomtom.loc[tomtom.Query_superclass == superclass].reset_index(drop=True)
    class_ = hierarchical_classification_tf(res_tomtom_filtr, 'Query_superclass', 'Target_superclass', 'Query_class', 'Target_class')
    family = hierarchical_classification_tf(class_, 'Query_class', 'Target_class', 'Query_family', 'Target_family')
    subfamily = hierarchical_classification_tf(family, 'Query_family', 'Target_family', 'Query_subfamily', 'Target_subfamily')
    gene = hierarchical_classification_tf(subfamily, 'Query_subfamily', 'Target_subfamily', 'Query_gene', 'Target_gene')
```

## Requirements

*Python >= 3.8

*pandas


*numpy



