import unittest
from os import path as osp
from pysimm import system, cassandra


class CassandraTestCase(unittest.TestCase):

    def test_gas_unwrap(self):

        # Read molecular system from the file and decorate the frame  with fixed property
        test_sst = system.read_lammps(osp.join('test_data', 'wrapped_co2.lmps'))
        css = cassandra.Cassandra(test_sst)
        for p in css.system.particles:
            if p.tag > 3500:
                p.is_fixed = False

        # Assert that gas molecules are wrapped
        bnd_lng = []
        for ind in [3501, 3504, 3507]:
            for i in [1, 2]:
                bnd_lng.append(
                    (css.system.particles[ind].x - css.system.particles[ind + i].x) ** 2 +
                    (css.system.particles[ind].y - css.system.particles[ind + i].y) ** 2 +
                    (css.system.particles[ind].z - css.system.particles[ind + i].z) ** 2)
        self.assertFalse(all([(1.344 < el) and (el < 1.346) for el in bnd_lng]))

        css.unwrap_gas()

        # Assert that gas molecules are unwrapped
        bnd_lng = []
        for ind in [3501, 3504, 3507]:
            for i in [1, 2]:
                bnd_lng.append(
                    (css.system.particles[ind].x - css.system.particles[ind + i].x) ** 2 +
                    (css.system.particles[ind].y - css.system.particles[ind + i].y) ** 2 +
                    (css.system.particles[ind].z - css.system.particles[ind + i].z) ** 2)
        self.assertTrue(all([(1.344 < el) and (el < 1.346) for el in bnd_lng]))


if __name__ == '__main__':
    unittest.main()
