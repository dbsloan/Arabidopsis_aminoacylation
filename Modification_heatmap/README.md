# Modification_heatmap

Generation of a heatmap showing base incorporations and deletion as an indication of post-transcriptional base modifications.

Input files with sequence variant data were generated based on workflow from [Edera and Sanchez-Puerta 2021](https://link.springer.com/protocol/10.1007/978-1-0716-0787-9_2). 

Two R scripts run in series to generate plot.

base-modification-profile-processing.revised.R
base-modification-profile-heatmap.revised.R

These are modified versions that collapse individual genes into isodecoder families. The original version (separate by gene) and required files are available in the `original_version_in_bioRxiv_preprint` directory.
