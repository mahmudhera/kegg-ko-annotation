# kegg-ko-annotation

We will annotate KO's of the KEGG database using sourmash and Diamond, and then
compare the results -- see how it goes.

Step 1: restrict the whole functional annotation pipeline within KEGG organisms.
Step 2: run the existing pipeline at https://github.com/KoslickiLab/KEGG_sketching_annotation.

Our main focus, for now, is to implement Step 1.

## Installation

```
conda create -n kegg_ko
conda install -y -n kegg_ko -c bioconda -c conda-forge -c conda --file requirements.txt
conda activate kegg_ko
```

## Workflow

Still not written, work in progress
