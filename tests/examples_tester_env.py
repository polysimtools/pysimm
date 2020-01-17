from os import path as osp
from os import chdir
from glob import glob
from imp import load_source

import unittest
import shutil
import re


class AbstractExamplesTestCase(unittest.TestCase):
    test_folder = 'Examples_to_test'

    @classmethod
    def setUpClass(cls):
        loc_path = osp.dirname(osp.relpath(__file__))
        tmp_data_dir = osp.join(loc_path, 'Examples_to_test')
        shutil.copytree(osp.join(loc_path, osp.pardir, 'Examples'), tmp_data_dir)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(osp.join(osp.dirname(osp.relpath(__file__)), 'Examples_to_test'))

    def setUp(self):
        self.loc_path = osp.dirname(osp.abspath(__file__))
        self.tmp_data_dir = osp.join(self.loc_path, 'Examples_to_test')

    def path_generator(self, path_mask, **kwargs):
        return glob(osp.join(self.tmp_data_dir, path_mask, kwargs.get('script_mask', '*.py')))

    def run_example(self, script_path):
        is_passed = True
        chdir(osp.dirname(script_path))
        src = load_source('create', osp.basename(script_path))
        try:
            src.run(test=True)
        except:
            is_passed = False
        chdir(self.loc_path)
        return is_passed


def example_tests_sort(x, y):
    tmp = re.search('\d+', x)
    a = 0
    if tmp:
        a = int(tmp.group(0))

    tmp = re.search('\d+', y)
    b = 0
    if tmp:
        b = int(tmp.group(0))
    return a - b
