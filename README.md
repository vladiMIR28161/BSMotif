### BSMotifs

BSMotif is a Python package representing a pipeline for cluster analysis based on
- the similarity of transcription factors (TF) binding site (TFBS) motifs, e.g., using the [Tomtom motif comparison tool](https://meme-suite.org/meme/tools/tomtom) ([Gupta et al., 2007](https://doi.org/10.1186/gb-2007-8-2-r24));
- TF classification by DNA-binding domain structure from the [TFClass database](http://www.edgar-wingender.de/huTF_classification.html) (Wingender et al., [2013](https://doi.org/10.1093/nar/gks1123), [2015](https://doi.org/10.1093/nar/gku1064), [2018](https://doi.org/10.1093/nar/gku1064)). 
BSMotif constructs the tree of similar TFBS motifs consisting of distinct branches. Within each branch, the majority of TF pairs have significantly similar binding site motifs. Each branch can include one or more sister elementary units of the hierarchy and all its/their lower levels: one or more TFs of the same subfamily, or the whole subfamily, one or several subfamilies of the same family, the whole family, etc. up to the whole class.

## Description

De novo motif search is the main approach to determine the nucleotide specificity of binding of the key regulators of gene transcription, TFs, based on data from massive genome-wide sequencing of their binding site regions in vivo, such as ChIP-seq. The number of motifs of known TFBSs has increased several times in recent years. Due to the similarity in the structure of the DNA-binding domains of TFs, many structurally cognate TFs have similar and sometimes almost indistinguishable binding site motifs. The classification of TFs by the structure of the DNA-binding domains from the TFClass database defines the top levels of the hierarchy (superclasses and classes of TFs) by the structure of these domains, and the next levels (families and subfamilies of TFs) by the alignments of amino acid sequences of domains. However, this classification does not take into account the similarity of TFBS motifs, whereas to identify valid TFs from massive sequencing data of TFBSs, such as ChIP-seq, one has to deal with TFBS motifs rather than TFs themselves.
The BSMotif determines the levels of the hierarchy (classes, families, subfamilies, or TFs), starting from which and lower in the TFClass hierarchy the known TFs become the significantly similar in terms of their TFBS motifs.
The application of BSMotif to the TFBS motif collections from the [Hocomoco](https://hocomoco.autosome.org/) and [Jaspar](https://jaspar.elixir.no/) databases was recently described in Levitsky et al. [2025](https://doi.org/10.AAA)

## Overview

Input data:

1. file Classification.tsv — TF classification by DNA-binding domains (including superclass, class, family, subfamily, gene).

2. path to the directory where the PWM files are stored in .meme format;

3. path to the directory where the output files will be located.

Output data:

1. file Tomtom_results.tsv — the result of executing the Tomtom.

2. file Score.tsv — table with motif similarity scores for each TF pair.

2. file Branches.tsv — list of TF branches of similar TFBS motifs.

Example:

* TFs from families XBP1-related {1.1.5} and CREB-related {1.1.7} are significantly similar → the pipeline merges them into one branch within the family level.

This allows identifying clusters of similar motifs within structural levels of DNA-binding domains.

## Installation

``` bash
git clone https://github.com/vladiMIR28161/BSMotif.git
cd bsmotif
pip install .
```

After installation, the bsmotif command becomes available in the terminal.

## Usage

``` bash
bsmotif <input_classification.tsv> <pfm_dir> <output_dir>
```

## Example code

``` python
from bsmotif.io import read_classification, read_tomtom, save_results
from bsmotif.preprocessing import fill_classification, calculate_score_tf
from bsmotif.hierarchical_classification import hierarchical_classification_tf
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

* Tomtom

* Python >= 3.8

* pandas

* numpy

## References

* Edgar Wingender, Torsten Schoeps, Martin Haubrock, Mathias Krull, Jürgen Dönitz, TFClass: expanding the classification of human transcription factors to their mammalian orthologs, Nucleic Acids Research, Volume 46, Issue D1, 4 January 2018, Pages D343–D347, https://doi.org/10.1093/nar/gkx987

* Edgar Wingender, Torsten Schoeps, Martin Haubrock, Jürgen Dönitz, TFClass: a classification of human transcription factors and their rodent orthologs, Nucleic Acids Research, Volume 43, Issue D1, 28 January 2015, Pages D97–D102, https://doi.org/10.1093/nar/gku1064

* Edgar Wingender, Torsten Schoeps, Jürgen Dönitz, TFClass: an expandable hierarchical classification of human transcription factors, Nucleic Acids Research, Volume 41, Issue D1, 1 January 2013, Pages D165–D170, https://doi.org/10.1093/nar/gks1123

* Gupta, S., Stamatoyannopoulos, J.A., Bailey, T.L. et al. Quantifying similarity between motifs. Genome Biol 8, R24 (2007). https://doi.org/10.1186/gb-2007-8-2-r24

* Ieva Rauluseviciute, Rafael Riudavets-Puig, Romain Blanc-Mathieu, Jaime A Castro-Mondragon, Katalin Ferenc, Vipin Kumar, Roza Berhanu Lemma, Jérémy Lucas, Jeanne Chèneby, Damir Baranasic, Aziz Khan, Oriol Fornes, Sveinung Gundersen, Morten Johansen, Eivind Hovig, Boris Lenhard, Albin Sandelin, Wyeth W Wasserman, François Parcy, Anthony Mathelier, JASPAR 2024: 20th anniversary of the open-access database of transcription factor binding profiles, Nucleic Acids Research, Volume 52, Issue D1, 5 January 2024, Pages D174–D182, https://doi.org/10.1093/nar/gkad1059

* Ilya E Vorontsov, Irina A Eliseeva, Arsenii Zinkevich, Mikhail Nikonov, Sergey Abramov, Alexandr Boytsov, Vasily Kamenets, Alexandra Kasianova, Semyon Kolmykov, Ivan S Yevshin, Alexander Favorov, Yulia A Medvedeva, Arttu Jolma, Fedor Kolpakov, Vsevolod J Makeev, Ivan V Kulakovskiy, HOCOMOCO in 2024: a rebuild of the curated collection of binding models for human and mouse transcription factors, Nucleic Acids Research, Volume 52, Issue D1, 5 January 2024, Pages D154–D163, https://doi.org/10.1093/nar/gkad1077

* Levitsky V. G., Vatolina T. Yu., Raditsa V. V. Relation between the hierarchical classification of transcription factors by the structure of their DNA-binding domains and the variability in the binding site motifs of these factors. Journal of Genetics and Breeding.


