# TAD Map
## Overview

With the advent of large-scale single-cell transcriptomic and cell type-specific chromatin conformation (HI-C) data, researchers can assess how 3D structure impacts transcription systematically across several cell types. Topologically associating domains (TADs) are contiguous segments of the genome where the genomic elements are in frequent contact with each other. The TAD Map provides a consensus estimate of this layout in human and mouse, aggregated from multiple experimental datasets. We also leverage recently developed single-cell foundation models that provide broad gene-gene coexpression estimates across a variety of cell types and cross-reference them against gene groupings designated by TAD Map. 

## Resources

The TAD scaffold and gene groups are available for download from our [website](https://singhlab.net/tadmap/). We provide several useful APIs to:
- Compute TAD signatures
- Directly retrieve the TAD Map
- Parse results from TAD Map files

For complete documentation, please visit our [documentation page](https://tadmap.readthedocs.io/en/latest/overview.html).

## Figures

In the figures folder, we have many scripts have been used to analyze TAD Map data and to create figures. The following are included:

## Contribute

We encourage you to report issues at our [Github page](https://github.com/rs239/tadmap); you can also create pull reports there to contribute your enhancements.

## Citation

If the TAD Map is useful in your research, please consider citing our preprint [bioRxiv (2021)](https://doi.org/10.1101/2021.10.28.466333).