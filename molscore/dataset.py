import logging

import numpy
import rdkit.Chem # type: ignore

from typing import Union, Iterable, Type

logger = logging.getLogger(__name__) 

from molscore.handler import DataHandler

class Dataset:
    
    def __str__(self):
        rep = "Dataset instance"
        if self.name is not None:
            rep = rep + f': `self.name`'
        rep += '\n'
        string_data = '\n'.join(list(self.data))
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
                    logger.info(f'"{instring}" could not be connonicalized, skipping')
                    pass
                else:
                    raise ValueError(f'"{instring}" could not be connonicalized')
        self._data = numpy.array(cansmiles).reshape(-1,1)
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
        dataset = cls(data=data, name=name)
        return dataset
        
    def _save_data_to_file(self, path: str):
        numpy.save(path, self.data)
        return