from __future__ import annotations

import json
import os
import shutil
import logging

from typing import Union, Iterable, Type

import numpy
import rdkit.Chem # type: ignore

import molscore

logger = logging.getLogger(__name__) 

class DataHandler:
    """Organizes interacts with many datasets of reactions on file.
    
    This class initializes a filesystem database at a specified directory.
    It then oversees adding datasets of molecules to the database, removing
    them, and updating them.
    
    Parameters
    ----------
    root : str
        The filepath to initialize the directory.
    restart : bool, default False
        Completely restart the state of the database - this will delete
        everything currently stored at `root`
        
    Attributes
    ----------
    initialized : bool
        Whether there is an initialized database at `self.root`.
        This will always return true after construction.
    metadata : dict
        Metadata for the database. Contents:
            'name': dict of {dataset_id, dataset_name}
    dataset_ids : list of int
        All dataset ids currently stored
    dataset_names : list of str
        All dataset names currently stored
    """
    def __init__(self, root: str = './molscore_data', restart: bool = False):
        # ensure no slash to save hassle later
        raise BaseException('ran handler')
        if root.endswith('/'):
            root = root[:-1]
        self.root = root
        
        if restart:
            shutil.rmtree(root, ignore_errors=True)
            logging.info(f'Existing database at `{root}` removed,')
            self._initialize()
        elif not self.initialized:
            self._initialize()
        else:
            file = open(self.root+'/metadata.json', 'r')
            self.metadata = json.load(file)
            file.close()
            logging.info(f'Existing database loaded at `{root}`.')
        return
    
    @property
    def initialized(self):
        """bool : the directory at `self.root` is initialized"""
        if os.path.exists(self.root+'/metadata.json'):
            metadata = True
        else:
            metadata = False
        return metadata
    
    def _initialize(self):
        """Performs conception of new database at `self.root`."""
        # create the folder
        os.mkdir(self.root)
        os.mkdir(self.root+'/data')
        # and the metadata
        self.metadata = {
            'names':{}
        }
        self._save_metadata()
        logging.info(f'Created database at `{self.root}`.')
        return
    
    def _save_metadata(self):
        """Update the file system with the metadata in this class."""
        file = open(self.root+'/metadata.json', 'w')
        json.dump(self.metadata, file)
        file.close()
        logging.debug(f'Metadata updated. Current state: \n{self.metadata}')
        return
    
    @property
    def dataset_ids(self):
        """list of int: ids of all datasets on file"""
        return list(self.metadata['names'].keys())
    
    @property
    def dataset_names(self):
        """list of str: names of all datasets on file"""
        return list(self.metadata['names'].values())
    
    def register_dataset(self, dataset: Dataset):
        """Register a dataset to this handler.
        
        Assigns an id and saves the dataset information to the filesystem.
        
        Parameters
        ----------
        dataset: Dataset
            The dataset to register.
            
        Returns
        -------
        int: the dataset id
        """
        # get newest id to assign to data
        if len(self.dataset_ids) < 1:
            dataset_id = 0
        else:
            largest_id = max(self.dataset_ids)
            # get smallest id not in use
            dataset_id = min(
                set(
                    range(largest_id+2)
                ) - set(
                    self.dataset_ids
                )
            )
        
        # if the dataset has a name, keep track of it
        if dataset.name == None:
            name = dataset_id
        else:
            if dataset.name in self.dataset_names:
                raise ValueError(
                    f'A dataset with name `{dataset.name}` already exists. Unregister that\
 dataset to continue with current registration. If you wish to update a\
 registered dataset, use `update_dataset` instead. '
                )
            name = dataset.name
        
        # update metadata
        dataset._save_data_to_file(self._get_path_to_dataset(dataset_id))
        self.metadata['names'][dataset_id] = name
        logging.info(f'Dataset registered with id and name: `{dataset_id}`, `{name}`.')
        logging.debug(f'Data in dataset `{dataset_id}`: \n {dataset.data}')
        self._save_metadata()
        return dataset_id
    
    def update_dataset(self, dataset: Dataset):
        """Update a specified dataset's info in the filesystem.
        
        Parameters
        ----------
        dataset: Dataset
            dataset to update, must be in the database already
        """
        if dataset._dataset_id is None or dataset._dataset_id not in self.dataset_ids:
            raise ValueError(f'Cannot update dataset, it is not registered.')
        else:
            dataset._save_data_to_file(self._get_path_to_dataset(dataset._dataset_id))
            if dataset.name is not None and dataset.name != self.metadata['names'][dataset._dataset_id]:
                self.metadata['names'][dataset._dataset_id] = dataset.name
                self._save_metadata()
            logging.info(f'Dataset `{dataset._dataset_id}` updated.')
            logging.debug(f'Data in dataset `{dataset._dataset_id}`: \n {dataset.data}')
        return
    
    def unregister_dataset(self, identifier: Union[int, str]):
        """Remove a dataset by name or id from the filesystem.
        
        This will look for ids before looking for names. If a dataset name that
        is an integer is passed, and it is the same as a dataset id, the dataset
        whose id matches `identifier` will be selected.
        
        Parameters
        ----------
        identifier: int or str
            The dataset id or dataset name to remove.
        """
        dataset_id = self._get_dataset_id_from_id_or_name(identifier)
        
        # remove from registry and the file system
        os.remove(self.root+f'/data/{dataset_id}.npy')
        logging.info(
            f'Dataset with id and name: `{dataset_id}`, `{self.metadata["names"][dataset_id]}` has been unregistered.')
        del self.metadata['names'][dataset_id]
        self._save_metadata()
        return
    
    def load_dataset(self, identifier: Union[str, int]):
        """Initialize a dataset that is saved to the file system.
        
        This will look for ids before looking for names. If a dataset name that
        is an integer is passed, and it is the same as a dataset id, the dataset
        whose id matches `identifier` will be selected.
        
        Parameters
        ----------
        identifier: str or int
            The dataset id or dataset name to load
        """
        dataset = Dataset.load(self, identifier)
        return dataset
    
    def _get_dataset_id_from_id_or_name(self, identifier: Union[str, int]):
        """Searches the database for a dataset id.
        
        Looks first in ids, then in names, and if it finds a dataset,
        returns the id.
        
        Parameters
        ----------
        identifier: str or int
            The dataset id or dataset name to look for
            
        Returns
        -------
        int: the dataset id found
        """
        # first try ids
        try:
            dataset_id = int(identifier)
        except:
            # could not convert input to int, must be the name of the dataset
            if identifier not in self.dataset_names:
                raise ValueError(
                    f'Could not interperet {identifier} as a dataset id or name'
                )
            else:
                dataset_id = self.dataset_ids[self.dataset_names.index(identifier)]
        
        if dataset_id not in self.dataset_ids:
            raise ValueError(
                f'Dataset id {dataset_id} not in the registry of datasets: {self.dataset_ids}'
            )
        return dataset_id
    
    def _get_path_to_dataset(self, dataset_id: int):
        """Retrieves the filepath that contains the molecule data for a dataset
        id."""
        return self.root+f'/data/{dataset_id}.npy'

class Dataset:
    """Represents a set of molecules.
    
    Molecules represented in the dataset will be canonicalized.
    
    Parameters
    ----------
    data: Iterable
        An iterable of strings representing molecules
    name: str, optional
        A name to give to this dataset.
    data_handler: DataHandler, optional
        The handler operating the filesystem to use. The default handler will
        will be used if none are passes.
    raise_smiles_errors: bool, default True
        Whether to raise errors when canonicalization fails. If False, skip
        failures.
        
    Attributes
    ----------
    name: str
        Name of this dataset
    data: Iterable[str]
        The canonicalized molecules in the dataset.
    raise_smiles_errors: bool
        Whether to raise errors when canonicalization fails. If False, skip
        failures.
    """
    def __str__(self):
        rep = "Dataset instance"
        if self.name is not None:
            rep = rep + f': `self.name`'
        rep += '\n\t'
        string_data = '\n\t'.join(list(self.data))
        rep += string_data
        return rep
    
    def __init__(
        self,
        data: Iterable[str],
        name: str = None,
        data_handler: DataHandler = None,
        raise_smiles_errors: bool = True
    ):
        self.name = name
        if data_handler is None:
            data_handler = self.get_default_handler()
        if not isinstance(data_handler, DataHandler):
            raise ValueError(
                f'`data_handler` must be of type DataHandler,\
 not {type(data_handler)}'
            )
        self.handler = data_handler
        
        # set the initial data id
        self._dataset_id = None
        
        # format and set the data
        self.raise_smiles_errors = raise_smiles_errors
        self.data = data
        return
    
    def register(self):
        """Register this dataset to filesystem.
        
        Uses the handler assigned to this instance.
        """
        self._dataset_id = self.handler.register_dataset(self)
        return
    
    def update(self):
        """Update this dataset in the filesystem.
        
        Uses the handler assigned to this instance.
        """
        self.handler.update_dataset(self)
        return
    
    def unregister(self):
        """Remove this dataset from the filesystem.
        
        Uses the handler assigned to this instance.
        """
        self.handler.unregister_dataset(self._dataset_id)
        return
    
    @property
    def name(self):
        """str : the name of this dataset"""
        return self._name
    
    @name.setter
    def name(self, new_name):
        if type(new_name) == str or new_name is None:
            self._name = new_name
        else:
            raise ValueError(
                f'Name must be of type `str` not {type(new_name)}')
        return
    
    @property
    def data(self):
        """Iterable[str]: The canonicalized molecules in the dataset."""
        return self._data
    
    @data.setter
    def data(self, d):
        d = numpy.array(d).reshape(-1)
        try:
            d = d.astype(str)
        except:
            raise ValueError(
                'All data inputs must be able to be converted to string'
            )
        
        cansmiles = []
        for instring in d:
            try:
                mol = rdkit.Chem.MolFromSmiles(instring)
                smiles = rdkit.Chem.MolToSmiles(mol)
                cansmiles.append(smiles)
            except:
                if not self.raise_smiles_errors:
                    logging.info(f'"{instring}" could not be connonicalized, skipping')
                    pass
                else:
                    raise ValueError(f'"{instring}" could not be connonicalized')
        self._data = numpy.array(cansmiles)
        return
    
    @staticmethod
    def get_default_handler():
        """Return the default handler that will be used if none are given."""
        return molscore.DEFAULT_HANDLER
        
    @classmethod
    def load(cls, handler: DataHandler, identifier: Union[str, int]):
        """Construction method to load a dataset from filesystem.
        
        This will look for ids before looking for names. If a dataset name that
        is an integer is passed, and it is the same as a dataset id, the dataset
        whose id matches `identifier` will be selected.
        
        Parameters
        ----------
        handler: DataHandler
            The handler to load the datsaet from.
        identifier: str or int
            The dataset id or dataset name to look for.
            
        Returns
        -------
        Dataset: the loaded dataset
        """
        # get the info on the dataset
        dataset_id = handler._get_dataset_id_from_id_or_name(identifier)
        path_to_load = handler._get_path_to_dataset(identifier)
        # load it
        data = numpy.load(path_to_load)
        name = handler.metadata['names'][dataset_id]
        dataset = cls(data=data, name=name, data_handler=handler)
        return dataset
        
    def _save_data_to_file(self, path: str):
        """Save the molecules in the dataset to a desired filepath as numpy.
        
        Parameters
        ----------
        path: str
            The path to save the data to.
        """
        numpy.save(path, self.data)
        return