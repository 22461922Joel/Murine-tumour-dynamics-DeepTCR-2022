# Murine-tumour-dynamics-DeepTCR-2022
DeepTCR scripts for classification of TCRs

## Citation
Kidman J et al., 2024, Oncoimmunology, Manuscript under review.
Sequencing data available GSE222575. Data will become publicly available once manuscript is published. Links will be updated following publication.

## Introduction
These scripts were adapted from [DeepTCR](sidhomj.github.io/deeptcr/) and are generated to suit bulk TCRB data from murine cancer models with response to immune checkpoint blockade. 
[TCR data](https://github.com/22461922Joel/Murine-tumour-dynamics-DeepTCR-2022/tree/72a0ed060befd0d69cacf64df0a0d2345e268873/data) was generated from BALB/c mice bearing mesothelioma (AB1) or renal cell carcinoma (RENCA). Using a bilaterial tumour model, mice treated with immune checkpoint therapy were identified as responders (RS) or non-rsponders (NR). Complete data information on murine model is described  [Zemek et al., Nat Comms 2020](https://www.nature.com/articles/s41467-022-32567-8)
Outputs generated were used in manuscript (Kidman J, 2024, Oncoimmunology, Manuscript under review) for publication. Refer to publication for full description of methods and results.

## Installation
Refer to [DeepTCR] (sidhomj.github.io/deeptcr/) for complete install information.
See requirements.txt to see package versions used to analyse this dataset.

## Pipeline
model_classifer.py & response_classifer.py must be run first before other scripts
Outputs from model_classifer.py and respones_classifer.py are required for downstream analyses.

## Figure breakdown
The below specifies which script was responsible for generating each figure in the manuscript (Kidman J et al., 2024 Oncoimmunology, under review).
Figure 4 - dynamic_signatures.py
Figure S5 - inf_response.py, inf_model.py, vis_corr_models.py


