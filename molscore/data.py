"""Handling and structure of data.

"""
import json
import os
import shutil
import warnings

import numpy
import rdkit.Chem # type: ignore

from typing import Union, Iterable, Type

class DataHandler:
    
    def __init__(self, root: str = './molscore_data', restart: bool = False):
        # ensure no slash to save hassle later
        if root.endswith('/'):
            root = root[:-1]
        self.root = root
        
        if restart:
            shutil.rmtree(root, ignore_errors=True)
            self._initialize()
        elif not self.initialized:
            self._initialize()
        else:
            file = open(self.root+'/metadata.json', 'r')
            self.metadata = json.load(file)
            file.close()
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
        os.makedir(self.root)
        os.makedir(self.root+'/data')
        # and the metadata
        self.metadata = {
            'names':{}
        }
        self._save_metadata()
        return
    
    def _save_metadata(self):
        file = open(self.root+'/metadata.json', 'w')
        json.dump(self.metadata, file)
        file.close
        return
    
    @property
    def dataset_ids(self):
        return list(self.metadata['names'].values())
    
    def resgister_data(self, dataset: 'Dataset'):
        # get newest id to assign to data
        if len(self.dataset_ids) < 1:
            dataset_id = 0
        else:
            dataset_id = int(max(self.dataset_ids)+1)
        
        # if the dataset has a name, keep track of it
        if dataset.name == None:
            name = dataset_id
        else:
            name = dataset.name
        numpy.save(self.root+f'/data/{dataset_id}.npy', dataset.data)
        
        # update metadata
        self.metadata['names'][name] = dataset_id
        self._save_metadata()
        return dataset_id
    
    
class Dataset:
    
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
        self._data_id = None
        
        # format and set the data
        self.raise_smiles_errors = raise_smiles_errors
        self.data = data
        return
    
    def register(self):
        self._data_id = self.handler.resgister_dataset(self)
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
                    warnings.warn(f'"{instring}" could not be connonicalized, skipping')
                    pass
                else:
                    raise ValueError(f'"{instring}" could not be connonicalized')
        self._data = numpy.array(cansmiles)
        return
    
    @staticmethod
    def get_default_handler():
        return DataHandler()
        
    