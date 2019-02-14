import unittest
import os
from os import path as osp
from testing import AbstractExamplesTestCase
from testing import example_tests_sort


class LammpsExamplesTestCase(AbstractExamplesTestCase):

    def test_lammps_binary(self):
        binary = os.environ.get('LAMMPS_EXEC')
        self.assertIsNotNone(binary)

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


if __name__ == '__main__':
    my_tl = unittest.TestLoader()
    my_tl.sortTestMethodsUsing = example_tests_sort
    unittest.TextTestRunner(buffer=True, verbosity=2).run(my_tl.loadTestsFromTestCase(LammpsExamplesTestCase))

