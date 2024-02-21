# LostOptionsCommitment

This repository contains the code and data necessary to reproduce the results in the research article [Lost options commitment: how short-term policies affect long-term scope of action](https://doi.org/10.1093/oxfclm/kgae004).

The computations in the mentioned article rely on the [SURFER model](https://gmd.copernicus.org/articles/15/8059/2022/) and make use of world total CO2 emissions for the SSP scenarios used for CMIP6 emission driven experiments. Emission data was originally retrieved from the [SSP](https://tntcat.iiasa.ac.at/SspDb/dsd?Action=htmlpage&page=10) and [RCP](https://tntcat.iiasa.ac.at/RcpDb/dsd?Action=htmlpage&page=welcome) Databases.

The code is written in [The Julia Programming Language](https://julialang.org/)

Several files are included in the repository:
- `README.md`: this file
- `LICENSE`: the license file
- `scenarios/`: a folder containing the SSP emission data required to run the provided code
- `load_ssp_scenarios.jl`: file that generates interpolating functions for the historic and SSP emission scenarios from the data available in the `emissions` folder.
- `SURFER-commit.jl`: the SURFER code with the addition of an extra function that allows to force the system with more sophisticated CO2 emissions and SO2 injections.
- `Lost options commitment.ipynb`: A Jupyter notebook that defines the necessary functions and shows how to reproduce the main results in the [article](https://doi.org/10.1093/oxfclm/kgae004).
