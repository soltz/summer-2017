from setuptools import setup

setup(
    name='HEPDataManager',
    version='0.1',
    description='Manager for local storage and loading of HEP data',
    author='Matthias Heinz',
    author_email='heinz.30@osu.edu',
    install_requires=[
        'future',
        'PyYAML'
    ],
    packages=[
        'hepdata_manager',
        'test'
    ]
)