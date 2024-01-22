# Virtual Metabolism Network Annotation Propagation (vm-NAP)

This notebook downloads results of spectral annotations from [classical molecular networking](https://ccms-ucsd.github.io/GNPSDocumentation/networking/) or [feature-based molecular networking ](https://ccms-ucsd.github.io/GNPSDocumentation/featurebasedmolecularnetworking/)job from GNPS [[http://gnps.ucsd.edu](http://gnps.ucsd.edu)] and generate virtual metabolites either with SyGMa or BioTransformer. The resulting candidates can be used for [Network Annotation Propagation](https://ccms-ucsd.github.io/GNPSDocumentation/nap/) on GNPS or with [SIRIUS](https://boecker-lab.github.io/docs.sirius.github.io/install/).


## Running vm-NAP

### vm-NAP web-app (cloud-based)

Click on the following link to launch the vm-NAP web-app. 
Note that this is a streamlit temporary instance with limited ressources.


### Local installation

Install locally in conda with:

>Download the present repository.

>In the terminal, navigate to the repository folder.

> Install the environment with:
`conda env create --file environment.yml`

> Initiate the environment:
`conda activate vm-NAP`

#### vm-NAP web-app locally

> Start the streamlit app with:

```
streamlit run vm_NAP_streamlit.py --server.port 8501 --server.address 0.0.0.0
```
> or

```
streamlit run vm_NAP_streamlit.py --server.port 8501 --server.address localhost
```

#### vm-NAP commandline:

> Representative command for the python script:

```
python src/vm_NAP_processing.py --job_id='bbee697a63b1400ea585410fafc95723' --run_sygma --run_biotransformer --sirius_input_file 'input/compound_identifications.tsv' --debug --max_compounds_debug=3
```

> Running this for help:

```
python src/vm_NAP_processing.py --help
```


### Jupyter notebook on Binder
The interactive notebook can be accessed via this badge -> [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/lfnothias/vm-NAP/main?urlpath=lab/tree/2401_vm-NAP-demo-notebook.ipynb)

Alternative - The interactive notebook can be accessed via this badge an gesis server-> [![Binder](https://mybinder.org/badge_logo.svg)](https://notebooks.gesis.org/binder/v2/gh/lfnothias/vm-NAP/main?urlpath=lab/tree/home/jovyan/2401_vm-NAP-demo-notebook.ipynb)

Note that this is also a temporary instance with limited ressources.
### Using vm-NAP with Network Annotation Propagation

See the documentation for custom database in [NAP](https://ccms-ucsd.github.io/GNPSDocumentation/nap/#structure-database) and how to run Network Annotation Propagation (NAP) on GNPS [https://ccms-ucsd.github.io/GNPSDocumentation/nap/#structure-database](https://ccms-ucsd.github.io/GNPSDocumentation/nap/#structure-database).

### Using vm-NAP With SIRIUS

See the documentation to generate the SIRIUS [custom database here](https://boecker-lab.github.io/docs.sirius.github.io/cli-standalone/#custom-database-tool).

## Spectral library requirement

*IMPORTANT: Note that only spectral annotations that have a valid InChI or SMILES identifier will be considered downstream. If the annotations you are interested in don't have an identifier in the library, go back to the GNPS library entry, update the entry by adding an identifier, and rerun your GNPS job*
