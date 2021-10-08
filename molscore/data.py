from __future__ import annotations

import json
import os
import shutil
import logging

from typing import Union, Iterable, Type

import numpy
import rdkit.Chem # type: ignore

logger = logging.getLogger(__name__) 

class DataHandler:
    
    def __init__(self, root: str = './molscore_data', restart: bool = False):
        # ensure no slash to save hassle later
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
        """Whether the directory at this location is initialized already."""
        if os.path.exists(self.root+'/metadata.json'):
            metadata = True
        else:
            metadata = False
        return metadata
    
    def _initialize(self):
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
        file = open(self.root+'/metadata.json', 'w')
        json.dump(self.metadata, file)
        file.close()
        logging.debug(f'Metadata updated. Current state: \n{self.metadata}')
        return
    
    @property
    def dataset_ids(self):
        return list(self.metadata['names'].keys())
    
    @property
    def dataset_names(self):
        return list(self.metadata['names'].values())
    
    def register_dataset(self, dataset: Dataset):
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
        if dataset._dataset_id is None or dataset._dataset_id not in self.dataset_ids:
            raise ValueError(f'Cannot update dataset, it is not registered.')
        else:
            dataset._save_data_to_file(self._get_path_to_dataset(dataset._dataset_id))
            logging.info(f'Dataset `{dataset._dataset_id}` updated.')
            logging.debug(f'Data in dataset `{dataset._dataset_id}`: \n {dataset.data}')
        return
    
    def unregister_dataset(self, identifier: Union[int, str]):
        dataset_id = self._get_dataset_id_from_id_or_name(identifier)
        
        # remove from registry and the file system
        os.remove(self.root+f'/data/{dataset_id}.npy')
        logging.info(
            f'Dataset with id and name: `{dataset_id}`, `{self.metadata["names"][dataset_id]}` has been unregistered.')
        del self.metadata['names'][dataset_id]
        self._save_metadata()
        return
    
    def load_dataset(self, identifier: Union[str, int]):
        dataset = Dataset.load(self, identifier)
        return dataset
    
    def _get_dataset_id_from_id_or_name(self, identifier):
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
        return self.root+f'/data/{dataset_id}.npy'

class Dataset:
    
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
        data: Iterable,
        name: str = None,
        data_handler: DataHandler = None,
        raise_smiles_errors: bool = True
    ):
        self.name = name
        if data_handler is None:
            data_handler = self.get_default_handler()
        self.handler = data_handler
        
        # set the initial data id
        self._dataset_id = None
        
        # format and set the data
        self.raise_smiles_errors = raise_smiles_errors
        self.data = data
        return
    
    def register(self):
        self._dataset_id = self.handler.register_dataset(self)
        return
    
    def update(self):
        self.handler.update_dataset(self)
        return
    
    def unregister(self):
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
        return DataHandler()
        
    @classmethod
    def load(cls, handler: DataHandler, identifier: Union[str, int]):
        # get the info on the dataset
        dataset_id = handler._get_dataset_id_from_id_or_name(identifier)
        path_to_load = handler._get_path_to_dataset(identifier)
        # load it
        data = numpy.load(path_to_load)
        name = handler.metadata['names'][dataset_id]
        dataset = cls(data=data, name=name, data_handler=handler)
        return dataset
        
    def _save_data_to_file(self, path: str):
        numpy.save(path, self.data)
        return