# enrichment_playground
Project for learning basic enrichment analysis

## GSEA implementation
As I wanted to run the analysis in Polyaxon (so that it would be faster), I needed to prepare the data for older R version - 3.6. This means I could not use tmod package. Therefore, the data for the GSEA algorithm is in the `polyaxon.RData` - tmod object is already converted into the data frame. Instead of the `source/gsea.R`, the `source/gsea_polyaxon.R` should be used. The code run for obtaining the results is available in the `job.R`.
