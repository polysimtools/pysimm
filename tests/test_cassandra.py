import unittest
from os import path as osp
from pysimm import system, cassandra


class CassandraTestCase(unittest.TestCase):

    def test_unwrap(self):

        # Read moleculear system from the file and decorate the frame  with fixed property
        test_sst = system.read_lammps(osp.join('test_data', 'wrapped_co2.lmps'))
        css = cassandra.Cassandra(test_sst)

        for p in css.system.particles:
            if p.tag > 3500:
                p.is_fixed = False

        css.unwrap_gas()
        self.assertIsNotNone(css)


if __name__ == '__main__':
    unittest.main()
