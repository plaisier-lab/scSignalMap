---
layout: default
title: Accessing Outputs
nav_order: 5
---
# Accessing Outputs
 Upon running `run_full_scSignalMap_pipeline`, returned objects should be as expected:

```r
names(results)
```

```r
[1] "LR_interactions"
[2] "de_cond_celltype"
[3] "de_receptors"
[4] "interactions_filtered"
[5] "de_receptors_filtered_and_compared"
[6] "enrichr_results"
```

Users can choose to compile results into a master interaction table using the function `master_interaction_list` below:
```r
master_interaction_list = create_master_interaction_list(
  enrichr_results = results$enrichr_results,
  de_receptors = results$upreg_receptors_filtered_and_compared,
  scSignalMap_data_filtered = results$interactions_filtered)
```

## Accessing Outputs
Accessing results from `run_full_scSignalMap_pipeline` can be accomplished using the following code:
```r
write.csv(results$<output>, "<file_name.csv>")
```
Accessing results from `create_master_interaction_list` can be accomplished using the following code:
```r
write.csv(master_interaction_list, "<file_name.csv>")
```
