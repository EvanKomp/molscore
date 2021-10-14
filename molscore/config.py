"""Handles configuration arguments by the user."""
import json
import os

import molscore

this_dir, this_filename = os.path.split(__file__)
CONFIG_PATH = os.path.join(this_dir, "config.json")

# load the current config
file = open(CONFIG_PATH, 'r')
_initial_config = json.load(file)
file.close()
VALID_CONFIG_PARAMETERS = list(_initial_config.keys())


def _check_config(config: dict):
    """A check to see that incoming parameters are valid"""
    if type(config['DEFAULT_DATABASE_ROOT']) != str:
        raise TypeError(
            f'Cannot set `DEFAULT_DATABASE_ROOT`, value must be string, not\
 {type(config["DEFAULT_DATABASE_ROOT"])}')
    if not os.path.exists(
        '/'.join(
            config['DEFAULT_DATABASE_ROOT'].split('/')[:-1]
        )
    ):
        raise ValueError(
            f'Cannot set database root to {config["DEFAULT_DATABASE_ROOT"]},\
 parent directory does not exist.'
        )
    return

def help():
    print(f"""Call update(parameter, new_value) to update config parameters.
    
    
    Valid parameters:
    {VALID_CONFIG_PARAMETERS}
    """)

def update(key: str, value):
    """Update a configuration parameter.
    
    Parameters
    ----------
    key : str
        The config parameter from to change
    value
        The value to change the parameter to.
    """
    # load the current config
    file = open(CONFIG_PATH, 'r')
    config = json.load(file)
    file.close()
    
    if key not in list(config.keys()):
        raise ValueError(
            f'`{key}` not a valid config parameter. Parameters available:\n{list(config.keys())}'
        )
    config[key] = value
    _check_config(config)
    
    # now save the new config to file
    file = open(CONFIG_PATH, 'w')
    json.dump(config, file)
    file.close()
    
    # now update the globals used by the package
    molscore._set_globals(config)
    return