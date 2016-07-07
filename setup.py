from setuptools import setup

setup(name='pysimm',
    version='0.1',
    description='Python Simulation Interface for Molecular Modeling',
    url='http://github.com/polysimtools/pysimm',
    author='Michael Fortunato',
    author_email='mefortunato@pysimm.org',
    license='MIT',
    packages=['pysimm', 'pysimm.apps'],
    package_data={'pysimm': ['forcefields/*.xml']},
    scripts=['bin/pysimm'],
    zip_safe=False)
