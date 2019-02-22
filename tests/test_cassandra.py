import unittest
from os import path as osp
from pysimm import system, cassandra
import pytest


class CassandraTestCase(unittest.TestCase):

    data_path = 'test_data'

    def test_gas_unwrap(self):

        # Read molecular system from the file and decorate the frame  with fixed property
        test_sst = system.read_lammps(osp.join(self.data_path, 'wrapped_co2.lmps'))
        css = cassandra.Cassandra(test_sst)
        for p in css.system.particles:
            p.is_fixed = False

        # Assert that gas molecules are wrapped
        bnd_lng = []
        for p in test_sst.particles:
            if p.type.name == 'C':
                for i in [1, 2]:
                    bnd_lng.append(
                        (css.system.particles[p.tag].x - css.system.particles[p.tag + i].x) ** 2 +
                        (css.system.particles[p.tag].y - css.system.particles[p.tag + i].y) ** 2 +
                        (css.system.particles[p.tag].z - css.system.particles[p.tag + i].z) ** 2)
        self.assertFalse(all([(1.344 < el) and (el < 1.346) for el in bnd_lng]))

        css.unwrap_gas()

        # Assert that gas molecules are unwrapped
        bnd_lng = []
        for p in test_sst.particles:
            if p.type.name == 'C':
                for i in [1, 2]:
                    bnd_lng.append(
                        (css.system.particles[p.tag].x - css.system.particles[p.tag + i].x) ** 2 +
                        (css.system.particles[p.tag].y - css.system.particles[p.tag + i].y) ** 2 +
                        (css.system.particles[p.tag].z - css.system.particles[p.tag + i].z) ** 2)
        self.assertTrue(all([(1.344 < el) and (el < 1.346) for el in bnd_lng]))

    def test_nonames_mc(self):
        sst = system.System()
        sst.dim = system.Dimension(dx=40, dy=40, dz=40, center=[0, 0, 0])
        sst.forcefield = 'trappe/amber'
        css = cassandra.Cassandra(sst)
        my_gcmc_props = css.read_input(osp.join(self.data_path, 'props.inp'))

        specie = system.read_lammps(osp.join(self.data_path, 'toluene_nonames.lmps'))
        specie.forcefield = 'trappe/amber'

        with pytest.raises(SystemExit) as err_back:
            css.add_gcmc(species=specie, is_rigid=True, max_ins=200,
                         chem_pot=-30.34, out_folder=self.data_path, **my_gcmc_props)
            assert err_back.type == SystemExit


if __name__ == '__main__':
    unittest.main()
