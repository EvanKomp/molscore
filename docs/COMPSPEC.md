# Component specification

`molscore.data.DataHandler`
Handles registering and retrieving of data and managing metadata for filesystem database.
*Inputs*: Directory to operate in, whether or not to restart
__Methods__
- `register_data`: Commits a dataset to the filesystem database
- - *Inputs*: A dataset
- - *Outputs*: An identifier for the dataset in the filesystem
- `get_data`: loads a dataset from the database
- - *Inputs*: Name or id of dataset
- - *Outputs*: Dataset
__Attributes__
- `root`: str, the location of the filesystem database
- `metadata`: dict, data associated with the dataset, such as names, id numbers, and location of datasets
- `initialized`: bool, whether the dataset has been initialized already
- `dataset_names`: list, names of registered datasets

`molscore.data.Dataset`
One set of molecules.
*Inputs*: An iterable of molecules to consider, the data handler to be used, otherwise the default will be used, a name for the dataset
__Methods__
- `register`: Have the data handler commit the data to the filesystem database.
- `get_default_handler`: Return the data handler to be used if none are specified.
__Attributes__
- `name`: str, name of this dataset
- `handler`: DataHandler, The data handler this dataset is assigned to
- `data`: Iterable, The iterable of molecules making up this dataset
- `raise_smiles_errors`: bool, whether to raise errors when converting inputs or to remove offending data

`molscore.map.Mapper`
Transforms one or more datasets into a dimensionality that allows them to be compared
*Inputs*: Datasets to use, mapping parameters
__Methods__
- `map`: Map the data
- - *Inputs*: optional reference map to use instead of fitting to data
- - *Outputs*: list, mapped representation of the datasets
- `register`: Save the mapped representation of the data
- `get_default_handler`: Return the data handler to be used if none are specified.
__Attributes__
- `datasets`: list of datasets, the datasets in the map
- `reference_map`: Object, The fit reference mapper used to transform the data if not fit directly
- `mapped_datasets`: list, mapped representations of data

`molscore.analyze.Scorer`
Compare mapped datasets
*Inputs*: maps to use as reference
__Methods__
- `score_breadth`: returns the breadth score for a 
- `score_depth`: returns the depth score for a dataset
