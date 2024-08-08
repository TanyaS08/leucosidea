# Leucosidea data

Repository to look at the effects of *Leucosidea* on the surrounding grassland community. 

## Repo organisation

This repository is contained in an Rproject and should be used as such (this is done by opening the .Rproj folder in RStudio). The scripts are contained in `scripts/` and are numbered sequentially based in the order within the methods in the manuscript. The `scripts/_internals.R` script is called in within scripts and is used to allow for on the fly data manipulation but still keeping the working scripts clear and concise.

Any figures that are generated will be written to `figures/` and data tables to `outputs/`, with the necessary data contained in `data/`.