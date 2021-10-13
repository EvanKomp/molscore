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
def patched_config_file(working_test_dir, session_mocker):
    print(working_test_dir)
    file = open(working_test_dir+'/config.json', 'w')
    file.write('\{"DEFAULT_DATABASE_ROOT": "'+str(working_test_dir)+'/data"\}')
    file.close()
    
    # patch this path to the package path
    session_mocker.patch('molscore.config.CONFIG_PATH', working_test_dir+'/config.json')
    session_mocker.patch('molscore._initial_config', {"DEFAULT_DATABASE_ROOT": f"{working_test_dir}/data"})
    return 
    