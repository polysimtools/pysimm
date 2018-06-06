import unittest
import shutil
import glob
import imp
import os
import re
from os import path as osp


class ExamplesTestCase(unittest.TestCase):

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
        return glob.glob(osp.join(self.tmp_data_dir, path_mask, kwargs.get('script_mask', '*.py')))

    def run_example(self, script_path):
        is_passed = True
        os.chdir(osp.dirname(script_path))
        src = imp.load_source('create', osp.basename(script_path))
        try:
            src.run(test=True)
        except:
            is_passed = False
        os.chdir(self.loc_path)
        return is_passed

    def test_example1(self):
        for p in self.path_generator(osp.join('01_methane', '*')):
            self.assertEqual(self.run_example(p), True)

    def test_example2(self):
        for p in self.path_generator(osp.join('02_methanol', '*')):
            self.assertEqual(self.run_example(p), True)

    def test_example3(self):
        for p in self.path_generator(osp.join('03_benzene', '*')):
            self.assertEqual(self.run_example(p), True)

    def test_example4(self):
        for p in self.path_generator(osp.join('04_polyethylene', '*')):
            self.assertEqual(self.run_example(p), True)

    def test_example5(self):
        for p in self.path_generator(osp.join('05_polyethylene-co-styrene', '*')):
            self.assertEqual(self.run_example(p), True)

    def test_example6(self):
        for p in self.path_generator(osp.join('06_polymethyl_methacryalte_multi', '*')):
            self.assertEqual(self.run_example(p), True)

    def test_example7(self):
        for p in self.path_generator(osp.join('07_lammps_simulation', '')):
            self.assertEqual(self.run_example(p), True)

    def test_example8(self):
        for p in self.path_generator(osp.join('08_ethanol_acetone_mixture', '')):
            self.assertEqual(self.run_example(p), True)

    def test_example9(self):
        for p in self.path_generator(osp.join('09_cassandra_simulations', '*')):
            self.assertEqual(self.run_example(p), True)

    def test_example10(self):
        for p in self.path_generator(osp.join('10_mof_swelling'), script_mask='run.py'):
            self.assertEqual(self.run_example(p), True)

if __name__ == '__main__':
    my_tl = unittest.TestLoader()
    my_tl.sortTestMethodsUsing = lambda x, y: int(re.search('\d+', x).group(0)) - int(re.search('\d+', y).group(0))
    unittest.main(testLoader=my_tl)
