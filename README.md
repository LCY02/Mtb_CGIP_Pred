# Deep learning - driven prediction of drug mechanism of action from large-scale chemical-genetic interaction profiles

![](https://github.com/LCY02/Mtb_CGIP_Pred/blob/main/Image/framework.png?raw=true)
> Framework

## File structure
``` bash
├── Codes/
│   ├── Gene_clustering/
│           ├── Hierarchical_clustering/
│                   └── hierarchical_clustering.R
│           └── Homology_search/
│                   ├── Ecoli_homologs_Rproj_020921.Rproj
│                   └── Mtub_Ecoli_homologs_gene.Rmd
│   └── Mtb_pred_ML/
│           ├── main.ipynb
│           └── data_processing.ipynb
├── Data/
│   ├── BLAST/
│           ├── Ecoli_K12_MG1655_GCF_000005845_2_ASM584v2_genomic.gff
│           ├── Ecoli_K12_MG1655_GCF_000005845_2_ASM584v2_protein.faa
│           ├── Ecoli_K12_MG1655_GCF_000005845_2_ASM584v2_protein.faa.phr
│           ├── Ecoli_K12_MG1655_GCF_000005845_2_ASM584v2_protein.faa.pin
│           ├── Ecoli_K12_MG1655_GCF_000005845_2_ASM584v2_protein.faa.psq
│           ├── Mtb_genes_of interest_names.txt
│           ├── Mtuberculosis_H37RvGCF_000195955_2_ASM19595v2_genomic.gff
│           └── Mtuberculosis_H37RvGCF_000195955_2_ASM19595v2_protein.faa
│   ├── CGIP/
│           └── 47K_SMILES_Zscore.csv
│   ├── Johnson_data/
│           └── README.md
│   └── ML_data/
│           ├── Data_split/
│                   ├── Features/
│                           ├── MB_test.npz
│                           ├── MB_train.npz
│                           ├── MB_val.npz
│                           ├── RN_test.npz
│                           ├── RN_train.npz
│                           └── RN_val.npz
│                   ├── Opt_hyperpars/
│                           ├── DMPNN_RN.json
│                           ├── FFN_MB.json
│                           ├── FFN_RN.json
│                           └── MPNN_RN.json
│                   ├── test.csv
│                   ├── train.csv
│                   └── val.csv
│           └── binary_13_clusters.csv
├── Image/
│   └── framework.png
├── Results/
│   ├── Homology_results/
│           ├── Ecoli_vs_Mtuberculosis_outfmt7.tsv
│           ├── Mtb_to_Ecoli_154_Ecoli_gff.txt
│           └── Mtub_to_Ecoli_154_Ecoli_gff.txt
│   ├── Mtb_inhibitors_pred/
│           ├── Mtb_inhibitors_DMPNN_preds.csv
│           ├── Mtb_inhibitors_FFN_MB_preds.csv
│           ├── Mtb_inhibitors_FFN_RN_preds.csv
│           └── Mtb_inhibitors_MPNN_preds.csv
│   └── Trained_model
│           ├── DMPNN_RN_Ensemble_5/
│                   └── ......
│           ├── FFN_MB_Ensemble_5/
│                   └── ......
│           ├── FFN_MB_Ensemble_5/
│                   └── ......
│           └── MPNN_RN_Ensemble_5/
│                   └── ......
├── README.md
```

## Summary
This repository contains scripts, datasets and partial results that support the findings of our study. The Johnson et al. data investigated in this study is publicly accessible on the [web site](https://www.chemicalgenomicsoftb.com/). The code to train the Message Passing Neural Networks is available on [Chemprop](https://github.com/chemprop).
