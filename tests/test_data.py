"""This module handles saving, retrieving, and organizing data."""
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
        handler = molscore.data.DataHandler(root=working_test_dir+'/root')
        DATASET.name = 'a set'
        assert handler.register_dataset(DATASET) == 0,\
            "Incorrect first dataset id returned"
        DATASET.name = None
        mocked_metadata_saver.reset_mock()
        assert handler.register_dataset(DATASET) == 1,\
            "Incorrect second dataset id returned"
        assert handler.metadata['names'] == {0: 'a set', 1: 1},\
            "register not recorded in metadata"
        assert mocked_metadata_saver.called,\
            "Save metadata not called"
        return