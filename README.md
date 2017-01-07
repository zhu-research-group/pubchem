# pubchem

### Tutorial 
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
