# Lineages are named by the detected label column

A run's lineages are always named in the label column detected from the dataset, also when its
terminals are named in another column. The API already did this. The CLI named lineages in the
terminals' column, and because the label column is not part of a run's identity, the same request
could carry different lineage names under one run id. We kept the API's rule because `launch`
shows it, every run it stored already records it, and it gives readable names (a cell type rather
than a cluster number). Naming in the terminals' column, which would echo what was asked, was
rejected.
