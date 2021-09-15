# `molscore` use cases

1. *Create a global molecule space.*
- __Components__: `DataHandler`
- User initializes a database of SMILES strings that act as a "global" smiles space

2. *Create a set of data (SMILES) cannonicalized and formated.*
- __Components__: `Dataset`
- User prodived an iterable of molecule strings, and they are formated into the packages's api

3. *Add data to the global molecule space*
- __Components__: `Dataset.register`
- Formatted data is registerd to the database to be part of the global molecule space.

4. *Map a set of one or more datasets to low dimensional molecule space*
- __Components__: `Mapper`
- A dataset is transformed to low dimensional molecule space

5. *Compare a set of data against another set of data*
- __Components__: `Scorer`
- Metrics between two sets of data are compared