# TCGA
# BRCA
# case/control
# LightGBM
# linear combination of existing prs
# common samples: 843
# normal samples: 75
# patients: 768
# R2: 0.9208219649580095

"""
This file is dedicate to improving the predictive accuracy of PRS model based on LightGBM  
by the combination of the existing results in each molecular dataset.
There are two kinds of methods to improve the accuracy.
One is to mix the results of each molecuar dataset.
Another is construct the new model which the new feature is the existing predictive PRS of each molecuar dataset.
"""


# -*- coding:utf-8 -*- 
import pandas as pd

# input the each PRS csv file
disease = "BRCA"
type = "case"
model = "LightGBM"
root = "/Users/panjianqiao/Desktop/PRS_data"
methylation_prs = pd.read_csv(root + "/" + disease +"_data" + "/TCGA_" + disease +"_" + type + "_control_methylation_lgb_prs.csv", index_col=0)
mirna_prs = pd.read_csv(root + "/" + disease +"_data" +"/TCGA_" + disease + "_" + type + "_control_miRNA_lgb_prs.csv", index_col=0)
mrna_prs = pd.read_csv(root + "/" + disease +"_data" +"/TCGA_" + disease + "_" + type + "_control_mRNA_lgb_prs.csv", index_col=0)
lncrna_prs = pd.read_csv(root + "/" + disease +"_data" +"/TCGA_" + disease + "_" + type + "_control_lncRNA_lgb_prs.csv", index_col=0)

# check the index of each csv file and select the common samples
methylation_prs_index = list(methylation_prs.index)
mirna_prs_index = list(mirna_prs.index)
mrna_prs_index = list(mrna_prs.index)
lncrna_prs_index = list(lncrna_prs.index)
common_index = set(methylation_prs_index) & set(mirna_prs_index) & set(mrna_prs_index) & set(lncrna_prs_index)
common_index = list(common_index)
print("The number of common samples: ", len(common_index))
normal_samples = 0
patients = 0
for i in range(len(common_index)):
    if common_index[i].split("-")[-1] == "01":  # patient is labelled 01
        patients += 1
    else:  # normal sample is labelled 11
        normal_samples += 1
print("The number of normal samples: ", normal_samples)  # 24
print("The number of patients: ", patients)  # 317

# linear combination of existing prs
prs = pd.concat([methylation_prs, mirna_prs, mrna_prs, lncrna_prs], axis=1, join="inner")
prs = prs.loc[:, ["methy_prs", "mirna_prs", "mrna_prs", "lncrna_prs", "methy_label"]]
la_prs = (prs.loc[:, "methy_prs"] + prs.loc[:, "mirna_prs"] + prs.loc[:, "mrna_prs"] + prs.loc[:, "lncrna_prs"])/4
prs["la_prs"] = la_prs
print(prs)

# R2
correlation = prs.methy_label.corr(prs.la_prs)
print(correlation ** 2)

prs_linear = prs.loc[:, ["methy_label", "la_prs"]]
prs_linear.rename(columns={"methy_label": "label", "la_prs": "clp"}, inplace=True)  # combination_linear_prs -> clp
order = ["clp", "label"]
prs_linear = prs_linear[order]
prs_linear.to_csv(root + "/" + disease + "_data/combination/TCGA_" + disease + "_" + type + "_control_linear_combination.csv")


