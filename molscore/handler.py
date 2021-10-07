import json
import os
import shutil
import logging

from typing import Union, Iterable, Type

from molscore.dataset import Dataset

logger = logging.getLogger(__name__) 


class DataHandler:
    
    def __init__(self, root: str = './molscore_data', restart: bool = False):
        # ensure no slash to save hassle later
        if root.endswith('/'):
            root = root[:-1]
        self.root = root
        
        if restart:
            shutil.rmtree(root, ignore_errors=True)
            logger.info(f'Existing database at `{root}` removed,')
            self._initialize()
        elif not self.initialized:
            self._initialize()
        else:
            file = open(self.root+'/metadata.json', 'r')
            self.metadata = json.load(file)
            file.close()
            logger.info(f'Existing database loaded at `{root}`.')
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
        logger.info(f'Created database at `{self.root}`.')
        return
    
    def _save_metadata(self):
        file = open(self.root+'/metadata.json', 'w')
        json.dump(self.metadata, file)
        file.close()
        logger.debug(f'Metadata updated. Current state: \n{self.metadata}')
        return
    
    @property
    def dataset_ids(self):
        return list(self.metadata['names'].keys())
    
    @property
    def dataset_names(self):
        return list(self.metadata['names'].values())
    
    def register_data(self, dataset: Dataset):
        # get newest id to assign to data
        if len(self.dataset_ids) < 1:
            dataset_id = 0
        else:
            largest_id = max(self.dataset_ids)
            # get smallest id not in use
            dataset_id = min(
                set(
                    range(largest_id+1)
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
        logger.info(f'Dataset registered with id and name: `{dataset_id}`, `{name}`.')
        logger.debug(f'Data in dataset `{dataset._dataset_id}`: \n {dataset.data}')
        self._save_metadata()
        return dataset_id
    
    def update_dataset(self, dataset: Dataset):
        if dataset._dataset_id is None or dataset._dataset_id not in self.dataset_ids:
            raise ValueError(f'Cannot update dataset, it is not registered.')
        else:
            dataset._save_data_to_file(self._get_path_to_dataset(dataset._dataset_id))
            logger.info(f'Dataset `{dataset._dataset_id}` updated.')
            logger.debug(f'Data in dataset `{dataset._dataset_id}`: \n {dataset.data}')
        return
    
    def unregister_dataset(self, identifier: Union[int, str]):
        dataset_id = _get_dataset_id_from_id_or_name(identifier)
        
        # remove from registry and the file system
        del self.metadata['names'][dataset_id]
        os.remove(self.root+f'/data/{dataset_id}.npy')
        logger.info(f'Dataset with id and name: `{dataset_id}`, `{name}` has been unregistered.')
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