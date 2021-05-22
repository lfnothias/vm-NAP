# Virtual Metabolism Network Annotation Propagation (vm-NAP)

This notebook downloads results of spectral annotations from [classical molecular networking](https://ccms-ucsd.github.io/GNPSDocumentation/networking/) or [feature-based molecular networking ](https://ccms-ucsd.github.io/GNPSDocumentation/featurebasedmolecularnetworking/)job from GNPS [[http://gnps.ucsd.edu](http://gnps.ucsd.edu)] and generate virtual metabolites either with SyGMa or BioTransformer. The resulting candidates can be used for [Network Annotation Propagation](https://ccms-ucsd.github.io/GNPSDocumentation/nap/) on GNPS or with [SIRIUS](https://boecker-lab.github.io/docs.sirius.github.io/install/).

## Running the notebook

View the notebook in [non-interactive view.](https://nbviewer.jupyter.org/github/lfnothias/vm-NAP/blob/main/2105_vm-NAP-GNPS.ipynb)

### One-click cloud usage with Binder [Recommended]

The interative notebook can be accessed via this badge -> [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/lfnothias/vm-NAP/main?urlpath=lab/tree/2105_vm-NAP-GNPS.ipynb). 

Note that for large computation (many compounds, many metabolization steps), this usage might be limited (2GB of max RAM).

### Local installation

Install locally in conda with:

>Download the present repository.

>In the terminal, navigate to the repository folder.

> Install the environment with:
`conda env vm-NAP --file environment.yml`


It might requires other depencies related ipython and jupter notebooks.

### Using vm-NAP with Network Annotation Propagation

See the documentation for custom database in [NAP](https://ccms-ucsd.github.io/GNPSDocumentation/nap/#structure-database) and how to run Network Annotation Propagation (NAP) on GNPS [https://ccms-ucsd.github.io/GNPSDocumentation/nap/#structure-database](https://ccms-ucsd.github.io/GNPSDocumentation/nap/#structure-database).

### Using vm-NAP With SIRIUS

See the documentation to generate the SIRIUS [custom database here](https://boecker-lab.github.io/docs.sirius.github.io/cli-standalone/#custom-database-tool).