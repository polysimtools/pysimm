import pysimm
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pysimm",
    version=pysimm.__version__,
    author="Mike Fortunato",
    author_email="mef231@gmail.com",
    description="python simulation interface for molecular modeling",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/polysimtools/pysimm/",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    include_package_data=True,
    package_data={
        'pysimm': [
            'data/*json',
            'data/*/*json',
            'data/*/*.xml',
            'data/*/*.inp'
        ]
    }
)
