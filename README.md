
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7034283.svg)](https://zenodo.org/badge/DOI/10.5281/zenodo.7034283.svg)

# islandbirds

Code and data for the paper:

Matthews et al (2022) Threatened and extinct island endemic birds: distribution, threats and functional diversity. Journal of Biogeography, in press

Code:
Two files - main script and source code. The source code needs to be loaded in first (top line of the main script)

Authors: Tom Matthews and Joe Wayman (Figure 7 code)


Data:

All data are saved in a single rds file: IUCN_Data.rds, which is part of this repository. This is a list with 7 elements (more info is provided in the scripts):

1) The extinct species data - species names, IUCN status, traits etc (data frame)

2) extant species data - species names, IUCN status, traits etc (data frame)

3) species IUCN threat data (list)

4) threatened and extinct species geographic data (data frame)

5) null model values for the beak morphospace analysis (using the all species pool) (list)

6) null model values for the beak morphospace analysis (using the island species pool) (list)

7) Hypervolume beta-div and null model results (list)


