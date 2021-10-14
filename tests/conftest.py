import json
import os
import pytest
import tempfile

@pytest.fixture(scope='session', autouse=True)
def working_test_dir():
    # create a produce a temporary directory to use for everything
    tmp_working_dir = tempfile.TemporaryDirectory()
    yield tmp_working_dir.name
    # delete it at the end of the session
    tmp_working_dir.cleanup()
    return

@pytest.fixture(scope='session', autouse=True)
def patched_config_file(working_test_dir):
    
    # now we have to manually modify the config file and replace it later
    this_dir, this_filename = os.path.split(__file__)
    config_path = os.path.join(this_dir, "../molscore/config.json")
    file = open(config_path, 'r')
    config_save = json.load(file)
    file.close()
    
    file = open(config_path, 'w')
    file.write('{"DEFAULT_DATABASE_ROOT": "'+str(working_test_dir)+'/data"}')
    file.close()
    yield None
    # now we have to save the old one back
    file = open(config_path, 'w')
    json.dump(config_save, file)
    file.close()
    return 
    