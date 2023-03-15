# Range and extinction dynamics of the steppe bison in Siberia: a pattern-oriented modeling approach

Authors: Julia A Pilowsky, Sean Haythorne, Stuart C Brown, Damien A Fordham.

This respository contains the `R` code and some example data (for 7 unique niche samples) to run the steppe bison model `paleopop` simulations. Please note that **not all the data are here** due to Github repository size limits. The data for mean and SD of human abundance, as well as the niche cuts, are available for download [on Zenodo](https://zenodo.org/record/7098687). The human data belong in the `Data` folder while the niche cuts should go in the `k_cuts` folder.

In order to run the example simulations, `paleopop v 2.0.0` or later must be installed. It is available on CRAN.

The code has a number of package dependencies listed below:

- "poems"
- "paleopop"
- "raster"
- "purrr"
- "stringr"
- "dplyr"
- "furrr"
- "data.table"
- "fs"
- "readxl"
- "rgdal"

The scripts are designed to run in parallel across *n* sessions - please set *n* accordingly.

There are two scripts contained in the repository. They need to be run in the following order:

- `bison_simulations.R`
- `bison_summary_metrics.R`

The first script build, runs, and saves the output from a multi-simulation manager (see `paleopop` documentation for details). The second script calculates summary metrics for the simulation outputs. These summary metrics can be compared against observed patterns for model selection by Approximate Bayesian Computation, as described in the main manuscript.

All outputs are saved in the `results/` directory.

Additional niche samples are available on reasonable request.

Directory structure:

```
Mammoth/
  - Data/ # data necessary to build the models
  - k_cuts/ # niche samples required for the models
  - results/ # output directory for simulation results
```

This code is released under the following licence:

GNU GENERAL PUBLIC LICENSE
   Version 3, 29 June 2007

Copyright (C) 2007 Free Software Foundation, Inc. <https://fsf.org/>
Everyone is permitted to copy and distribute verbatim copies
of this license document, but changing it is not allowed.
