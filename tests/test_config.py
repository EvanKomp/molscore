"""Ensures that the config for the package is handled correctly."""
import io

import pytest
import unittest.mock as mock


fake_io_file1 = io.StringIO('{"DEFAULT_DATABASE_ROOT": "./molscore_data"}')
fake_io_file2 = io.StringIO('')
fake_io_file2.close = lambda: None

def test_local_variables():
    import molscore.config
    """ensure the config variables are initialized"""
    assert type(molscore.config._initial_config) == dict,\
        "Config should be loaded as dict"
    assert type(molscore.config.VALID_CONFIG_PARAMETERS) == list,\
        "Config param options should be list of str."
    return


def test__check_config():
    import molscore.config
    """Another layer of protection for ensuring configs are handled."""
    # good dict, valid config params
    # nothing should happen, eg nothing raised
    good = {'DEFAULT_DATABASE_ROOT': './'}
    molscore.config._check_config(good)
    
    # bad dict, incorrect value
    bad = {'DEFAULT_DATABASE_ROOT': 5}
    with pytest.raises(TypeError):
        molscore.config._check_config(bad)
    # bad dict, folder does not exist
    bad = {'DEFAULT_DATABASE_ROOT': './not_a_real_folder_hopefully/my_data'}
    with pytest.raises(ValueError):
        molscore.config._check_config(bad)
    return

@mock.patch('molscore._set_globals')
@mock.patch('molscore.config.open', side_effect=[fake_io_file1, fake_io_file2])
def test_update(mocked_open, mocked_global_setter):
    import molscore.config
    """Update config without actually doing so."""
    molscore.config.update('DEFAULT_DATABASE_ROOT', './new_location')
    
    assert mocked_global_setter.called,\
        """Did not update global variables"""
    assert fake_io_file2.getvalue() == '{"DEFAULT_DATABASE_ROOT": "./new_location"}',\
        "Did not save to file the new variable"
    return
    