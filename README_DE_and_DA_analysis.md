# DE analysis

For Differential Gene Expression (DE) analysis see the jupyter notebook: Downstream_DE_analysis_without_time_as_covariate_v2.ipynb. The conda environment is multiome_NEW_downstream_analysis (see Dependencies/multiome_NEW_downstream_analysis.yml).

# DA analysis

For the Differential Accesibility analysis see:

## 1. Link Peaks to DE genes

$ bash ~/Scripts/Wraper_scripts/133_Link_peaks.sh /group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/ Downstream_analysis

## 2. Analyze DA in those peaks

Jupyter notebook: Downstream_DA_analysis_without_time_as_covariate.ipynb

## 3. For global analysis of DA parallelized by cell type see:

$ bash ~/Scripts/Wraper_scripts/132_PB_DA_CHEK2.sh /group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/ Downstream_analysis