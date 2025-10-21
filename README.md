### BSMotifs

BSMotifs is a Python package and pipeline for analyzing transcription factor (TF) binding site motif similarity. It integrates motif comparison results with TF classification by DNA-binding domain and constructs hierarchical branches of similar motifs within different classification levels.

## Overview

The pipeline takes two input files:

1. Tomtom.tsv — results of pairwise motif comparisons between transcription factors (e.g., using Tomtom).

2. Classification.tsv — transcription factor classification by DNA-binding domain (including superclass, class, family, subfamily, gene).

Outputs:

1. Score.tsv — table with motif similarity scores for each TF pair.

2. Classification.xlsx — hierarchical list of branches of similar motifs within each classification level (class, family, subfamily, TF).

Example:

* TFs from families XBP1-related {1.1.5} and CREB-related {1.1.7} are significantly similar → the pipeline merges them into one branch within the family level.

This allows identifying clusters of similar motifs within structural levels of DNA-binding domains.

## Installation

``` bash
git clone https://github.com/vladiMIR28161/BSMotifs.git
cd bsmotifs
pip install
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

