import unittest
import os
from os import path as osp
from testing import AbstractExamplesTestCase
from testing import example_tests_sort


class CassExamplesTestCase(AbstractExamplesTestCase):

    def test_cass_binary(self):
        binary = os.environ.get('CASSANDRA_EXEC')
        self.assertIsNotNone(binary)

    def test_example9(self):
        for p in self.path_generator(osp.join('09_cassandra_simulations', '*')):
            self.assertEqual(self.run_example(p), True)

    def test_example10(self):
        for p in self.path_generator(osp.join('10_mof_swelling'), script_mask='run.py'):
            self.assertEqual(self.run_example(p), True)


if __name__ == '__main__':
    my_tl = unittest.TestLoader()
    my_tl.sortTestMethodsUsing = example_tests_sort
    unittest.TextTestRunner(buffer=True, verbosity=2).run(my_tl.loadTestsFromTestCase(CassExamplesTestCase))

