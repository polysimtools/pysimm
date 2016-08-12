#!/bin/bash

sphinx-apidoc -fe -o . ../../pysimm/
make html
