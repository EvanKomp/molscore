__version__='0.0.1'
from .config import _initial_config

from .data import Dataset, DataHandler
print(_initial_config)

# Load the user configurations
def _set_globals(config: dict):
    """Set the global config parameters."""
    for key, value in config.items():
        globals()[key] = value
        
    # Initialize the default handler
    default_handler = DataHandler(
        root = DEFAULT_DATABASE_ROOT,
        restart=False
    )
    globals()['DEFAULT_HANDLER'] = default_handler
    return
_set_globals(_initial_config)