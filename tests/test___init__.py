"""Checks global variables initialization."""
import pytest


def test_globals(working_test_dir):
    import molscore
    molscore._set_globals(
        {'DEFAULT_DATABASE_ROOT': f'{working_test_dir}/another_data'})
    assert molscore.DEFAULT_HANDLER.root == f'{working_test_dir}/another_data',\
        "Did not change default handler"
    return
