"""This module handles saving, retrieving, and organizing data."""
import numpy as np
import os
import pytest
import shutil

import unittest.mock as mock

@pytest.fixture
def DATASET():
    import molscore.data
    dataset = molscore.data.Dataset([
        'C',
        'O=C=O'
    ])
    return dataset

class TestDataHandler:
    
    def test___init__(self, working_test_dir):
        """this is the saving of data to a specified location."""
        import molscore.data
        
        
        handler = molscore.data.DataHandler(root=working_test_dir+'/root')
        assert os.path.exists(working_test_dir+'/root/metadata.json'),\
            "Did not create metadata."
        assert os.path.exists(working_test_dir+'/root/data'),\
            "Did not create data dir."
        assert handler.root == working_test_dir+'/root',\
            "Not keeping track of root dir"
        #test the restart keyword
        file = open(working_test_dir+'/root/data/tmp.txt', 'w')
        file.close()
        handler = molscore.data.DataHandler(root=working_test_dir+'/root', restart=True)
        assert not os.path.exists(working_test_dir+'/root/data/tmp.txt'),\
            "Restart did not delete old data."
        
        # see if it will load back up
        handler.metadata['thing'] = 1
        handler._save_metadata()
        handler = molscore.data.DataHandler(root=working_test_dir+'/root')
        assert 'thing' in handler.metadata,\
            "Did not load metadata from previously"
        return
    
    def test__initialize(self, working_test_dir):
        import molscore.data
        handler = molscore.data.DataHandler(root=working_test_dir+'/root')
        assert handler.initialized, "init not recognized"
        shutil.rmtree(handler.root)
        assert not handler.initialized,\
            "handler did not recognize that it was not initialized"
        handler._initialize()
        assert os.path.exists(working_test_dir+'/root/metadata.json'),\
            "_initialize did not create metadata"
        return
    
    @mock.patch('molscore.data.DataHandler._save_metadata')
    def test_register_dataset(self, mocked_metadata_saver, working_test_dir, DATASET):
        import molscore.data
        handler = molscore.data.DataHandler(root=working_test_dir+'/root', restart=True)
        DATASET.name = 'a set'
        handler.register_dataset(DATASET)
        assert  DATASET.dataset_id == 0,\
            "Incorrect first dataset id returned"
        DATASET.name = None
        mocked_metadata_saver.reset_mock()
        handler.register_dataset(DATASET)
        assert DATASET.dataset_id == 1,\
            "Incorrect second dataset id returned"
        assert handler.metadata['names'] == {0: 'a set', 1: 1},\
            "register not recorded in metadata"
        assert mocked_metadata_saver.called,\
            "Save metadata not called"
        assert os.path.exists(working_test_dir+'/root/data/0.npy')
        return
    
    @mock.patch('molscore.data.DataHandler._save_metadata')
    def test_update_dataset(self, mocked_metadata_saver, working_test_dir, DATASET):
        import molscore.data
        handler = molscore.data.DataHandler(root=working_test_dir+'/root', restart=True)
        DATASET.name = 'a set'
        handler.register_dataset(DATASET)
        mocked_metadata_saver.reset_mock()
        
        # rename and change data to resave
        DATASET.name = 'new name'
        DATASET.data = np.append(DATASET.data, ['CCC'])
        handler.update_dataset(DATASET)
        assert mocked_metadata_saver.called,\
            "Save metadata not called"
        assert handler.metadata['names'][0] == 'new name',\
            "New name not saved"
        assert len(np.load(working_test_dir+'/root/data/0.npy')) == 3,\
            "New data not saved"
        
        # see if it can catch not a registered dataset
        DATASET._dataset_id = None
        with pytest.raises(ValueError):
            handler.update_dataset(DATASET)
        return
    
    @mock.patch('molscore.data.DataHandler._save_metadata')
    def test_unregister_dataset(self, mocked_metadata_saver, working_test_dir, DATASET):
        import molscore.data
        handler = molscore.data.DataHandler(root=working_test_dir+'/root', restart=True)
        DATASET.name = 'a set'
        handler.register_dataset(DATASET)
        mocked_metadata_saver.reset_mock()
        
        # just unregister it and check
        handler.unregister_dataset(0)
        assert mocked_metadata_saver.called,\
            "Save metadata not called"
        assert len(handler.metadata['names']) == 0,\
            "dataset was not removed from metadata"
        assert not os.path.exists(working_test_dir+'/root/data/0.npy'),\
            "dataset was not removed from disk"
        return
    
    @mock.patch('molscore.data.Dataset.load')
    def test_load_dataset(self, mocked_loader, working_test_dir, DATASET):
        import molscore.data
        handler = molscore.data.DataHandler(root=working_test_dir+'/root', restart=True)
        DATASET.name = 'a set'
        handler.register_dataset(DATASET)

        dataset = handler.load_dataset('a set')
        assert mocked_loader.calledwith('a set'),\
            "Did not call loader"
        return
    
    def test_hidden_methods(self, working_test_dir, DATASET):
        import molscore.data
        handler = molscore.data.DataHandler(root=working_test_dir+'/root', restart=True)
        print(handler.metadata)
        DATASET.name = 'a name'
        handler.register_dataset(DATASET)
        
        # getting id from name or id
        assert handler._getdataset_id_from_id_or_name(0) == 0,\
            "cannot get id from id"
        assert handler._getdataset_id_from_id_or_name('a name') == 0,\
            "cannot get id from name"
        
        # get the path from id
        assert handler._get_path_to_dataset(0) == working_test_dir+'/root/data/0.npy',\
            "wrong path returned"
        
        #saving metadata
        os.remove(working_test_dir+'/root/metadata.json')
        handler._save_metadata()
        assert os.path.exists(working_test_dir+'/root/metadata.json'),\
            "Did not save metadata"
        return
        
        
class TestDataset:
    
    def test___init__(self, working_test_dir):
        import molscore.data
        handler = molscore.data.DataHandler(root=working_test_dir+'/root', restart=True)
        
        dataset = molscore.data.Dataset(
            ['C', 'CC'],
            data_handler=handler,
            name='name')
        assert len(dataset.data) == 2, "did not store data"
        assert dataset.handler is handler, "not using correct handler"
        assert dataset.name == 'name'
        
        # now check default handler
        dataset = molscore.data.Dataset(
            ['C', 'CC'])
        assert dataset.handler is dataset.get_default_handler() is molscore.DEFAULT_HANDLER,\
            "Did not use default handler"
        return
    
    @mock.patch('molscore.data.DataHandler.register_dataset')
    def test_register(self, mocked_register, working_test_dir, DATASET):
        import molscore.data
        DATASET.register()
        mocked_register.assert_called_with(DATASET)
        return
    
    @mock.patch('molscore.data.DataHandler.update_dataset')
    def test_update(self, mocked_update, working_test_dir, DATASET):
        import molscore.data
        DATASET.register()
        DATASET.update()
        mocked_update.assert_called_with(DATASET)
        return
    
    @mock.patch('molscore.data.DataHandler.unregister_dataset')
    def test_unregister(self, mocked_unregister, working_test_dir, DATASET):
        import molscore.data
        DATASET.register()
        DATASET.unregister()
        mocked_unregister.assert_called_with(DATASET.dataset_id)
        return
    
    def test_load(self, working_test_dir, DATASET):
        import molscore.data
        handler = molscore.data.DataHandler(root=working_test_dir+'/root', restart=True)
        DATASET.handler = handler
        DATASET.register()
        new_dataset = molscore.data.Dataset.load(handler, 0)  
        assert new_dataset == DATASET,\
            "loaded dataset is not the same"
        return
    
    def test__save_data_to_file(self, working_test_dir, DATASET):
        DATASET._save_data_to_file(working_test_dir+'/test.npy')
        assert os.path.exists(working_test_dir+'/test.npy'),\
            "Did not save data array"
        return