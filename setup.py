from setuptools import setup

setup(name='pysimm',
    version='0.1',
    description='Python Simulation Interface for Molecular Modeling',
    url='http://github.com/polysimtools/pysimm',
    author='Michael E. Fortunato, Coray M. Colina',
    author_email='mefortunato@pysimm.org, colina@chem.ufl.edu',
    license='MIT',
    packages=['pysimm', 'pysimm.apps', 'pysimm.forcefield'],
    include_package_data=True,
    scripts=['bin/pysimm'],
    zip_safe=False)
