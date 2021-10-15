import setuptools
from os import path
import molscore

here = path.abspath(path.dirname(__file__))
AUTHORS = """
Evan Komp
"""

LICENSE = open('LICENSE').read()

# Get the long description from the README file
with open(path.join(here, 'README.md')) as f:
    long_description = f.read()

if __name__ == "__main__":
    setuptools.setup(
        name='molscore',
        version=molscore.__version__,
        author=AUTHORS,
        project_urls={
            'Source': 'https://github.com/EvanKomp/molscore/',
        },
        description='Package to score molecular datasets.',
        long_description=long_description,
        include_package_data=False,
        keywords=[
            'Machine Learning',
            'Chemical Engineering','Chemistry', 
        ],
        license=LICENSE,
        packages=setuptools.find_packages(exclude="tests"),
        scripts = [], #if we want to include shell scripts we make in the install
        install_requires=[
            'numpy', 
        ],
        extras_require={
            'tests': [
                'pytest',
                'pytest-mock',
                'coverage',
                'flake8',
                'flake8-docstrings',
                'yapf',
                'mypy',
                'typed_ast'
            ],
        },
        classifiers=[
            'Development Status :: 1 - Planning',
            'Environment :: Console',
            'Operating System :: OS Independant',
            'Programming Language :: Python',
            'Topic :: Scientific/Engineering',
        ],
        zip_safe=False,
    )
