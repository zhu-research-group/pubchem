# PubChem

# Tutorials

### AID Class
```python

# import AID class
from aids import AID

# instantiate class with AID number of interest
aid = AID(1)

# get interested property from XML
# examples, 'Name' or 'Description'
aid.get_property('Name')
```
['NCI human tumor cell line growth inhibition assay. Data for the NCI-H23 Non-Small Cell Lung cell line']

### Command Line Script to Generate and Plot Bioprofile

This script will pull bioassay information from PubChem for a specified list of compounds, make a bioprofile, and plot a clustered heatmap.

Running this code requires an Anaconda distribution for Python. To create the necessary Python environment, run the following command in a Windows command window or Bash shell:

```conda env create -f environment.yml```

If that causes an error, you can instead run the following command:

```conda create -n pubchem pandas requests seaborn```

After the necessary dependencies are installed, you must activate this environment before running these scripts. To do that, run the following command.

On Windows:
```activate pubchem```

On Linux/OS:
```source activate pubchem```

This script requires a .csv file containing a type of identifier for each compound of interest (PubChem CID, InChIKey, Canonical SMILES, CASRN).

This script accepts five parameters:
1. `-c`: Indicates whether or not you are running on the cluster (kestrel or Amarel); Type in 'y' if you are using the cluster and 'n' otherwise.
2. `-df`: Name of .csv file containing identifiers
3. `-ev`: Name of environment variable containing files. If this isn't specified, the script will assume the files are in the same directory as the scripts.
4. `-i`: Type of identifier used
5. `-ma`: Minimum number of active responses across compounds to keep an assay in the bioprofile

An example code to run this script on a local PC on a file `my_cids.csv` containing PubChem CIDs for my compounds and keeping assays with at least 10 active responses is:

```python bioprofile.py -c n -df my_cids -i CID -ma 10```

This script creates three files in the same directory as the input file:
1. `(your datafile name)_bioprofile.csv`: A matrix of 1, 0, and -1 corresponding to active, inconclusive/missing, and inactive responses, respectively
2. `(your datafile name)_bioprofile_heatmap.png`: An image file of a heatmap based on the bioprofile that is clustered by compounds and assays
3. `(your datafile name)_skipped_cpds.csv`: A list of compounds that were skipped because of a web error or inability to find them on PubChem with the listed identifier (as a .csv file; if no compounds were skipped, this file will be absent.)
