import unittest
import os
from pysimm import system
from pysimm import forcefield


class SystemReadTestCase(unittest.TestCase):
    TEST_DATA_PATH = 'test_data'

    def test_read_xyz(self):
        file_path = os.path.join(self.TEST_DATA_PATH, 'testfile.xyz_input.xyz')

        with open(file_path, 'r') as pntr:
            stream = pntr.read()
            # Asserting that at least something was read
            self.assertIsNotNone(stream)
            part_len = int(stream.split('\n')[0])
            test_sst = system.read_xyz(file_path)

            # Asserting that number of particles read coinside with XYZs control number
            self.assertEqual(len(test_sst.particles), part_len)

    def test_read_cml(self):
        file_path = os.path.join(self.TEST_DATA_PATH, 'testfile.cml_input.cml')
        test_sst = system.read_cml(file_path)

        # Asserting that number of particles and bonds read from .cml test file is greater than 0
        self.assertGreater(len(test_sst.particles), 0)
        self.assertGreater(len(test_sst.bonds), 0)

    def test_read_pdb(self):
        file_path = os.path.join(self.TEST_DATA_PATH, 'testfile.pdb_input.pdb')
        test_sst = system.read_pdb(file_path)

        # Asserting that number of particles read from .pdb test-file is greater than 0
        self.assertGreater(len(test_sst.particles), 0)

    def test_read_v2000mol(self):
        # check the v2000 format mol file
        file_path = os.path.join(self.TEST_DATA_PATH, 'testfile.molv2000_input.mol')
        test_sst = system.read_mol(file_path, version='V2000')
        self.assertGreater(len(test_sst.particles), 0)
        self.assertGreater(len(test_sst.bonds), 0)

    def test_read_v3000mol(self):
        # check the v2000 format mol file
        file_path = os.path.join(self.TEST_DATA_PATH, 'testfile.molv3000_input.mol')
        test_sst = system.read_mol(file_path, version='V3000')
        self.assertGreater(len(test_sst.particles), 0)
        self.assertGreater(len(test_sst.bonds), 0)

    def test_read_cdjson(self):
        # check the reading of the json file
        file_path = os.path.join(self.TEST_DATA_PATH, 'testfile.cd_json_input.json')
        test_sst = system.read_chemdoodle_json(file_path)
        self.assertGreater(len(test_sst.particles), 0)
        self.assertGreater(len(test_sst.bonds), 0)

    def test_read_pubchem_smiles(self):
        test_sst = system.read_pubchem_smiles('CC1=CC=C(C=C1)C')
        self.assertGreater(len(test_sst.particles), 0)
        self.assertGreater(len(test_sst.bonds), 0)

    def test_read_pubchem_cid(self):
        test_sst = system.read_pubchem_cid('7809')
        self.assertGreater(len(test_sst.particles), 0)
        self.assertGreater(len(test_sst.bonds), 0)

    def test_read_ac(self):
        file_path = os.path.join(self.TEST_DATA_PATH, 'testfile.ac_input.ac')
        test_sst = system.read_ac(file_path)
        self.assertGreater(len(test_sst.particles), 0)
        self.assertGreater(len(test_sst.bonds), 0)

    def test_read_prepc(self):
        file_path = os.path.join(self.TEST_DATA_PATH, 'testfile.prepc_input.prepc')
        test_sst = system.read_prepc(file_path)
        self.assertGreater(len(test_sst.particles), 0)

    # Test for reading lammps file of class-II FF
    def test_read_lammps(self):
        file_path = os.path.join(self.TEST_DATA_PATH, 'testfile_class2FF.lmps')
        test_sst = system.read_lammps(file_path)
        self.assertGreater(len(test_sst.particles), 0)
        self.assertGreater(len(test_sst.bonds), 0)

    # Test for both read_lammps and connected read_lammps_trj
    def test_read_lammpstrj(self):
        file_path = os.path.join(self.TEST_DATA_PATH, 'testfile.lammpstrj.$$')
        test_sst = system.read_lammps(file_path.replace('$$', 'lmps'))

        # Verifying that the system was read from the lammps file
        self.assertGreater(len(test_sst.particles), 0)
        self.assertGreater(len(test_sst.bonds), 0)

        # Checking that read_lammpstrj actually updates the coordinates of the atoms in the system
        test_sst.read_lammpstrj(file_path.replace('$$', 'dump'), frame=1)
        tmp_coords_before = [test_sst.particles[1].x, test_sst.particles[1].y, test_sst.particles[1].z]
        test_sst.read_lammpstrj(file_path.replace('$$', 'dump'), frame=3)
        tmp_coords_after = [test_sst.particles[1].x, test_sst.particles[1].y, test_sst.particles[1].z]

        self.assertNotEqual(tmp_coords_before[0], tmp_coords_after[0]) and \
        self.assertNotEqual(tmp_coords_before[1], tmp_coords_after[1]) and \
        self.assertNotEqual(tmp_coords_before[2], tmp_coords_after[2])


class SystemWriteTestCase(unittest.TestCase):

    def setUp(self):
        self.sst = system.read_pubchem_cid(6360)
        self.sst.dim = system.Dimension(xlo=-10, ylo=-10, zlo=-10, xhi=10, yhi=10, zhi=10)
        self.sst.apply_forcefield(forcefield.Pcff(), charges='gasteiger')
        self.fname_temp = 'testfile.output.{}'
        self.fname = ''

    def test_write_lammps(self):
        self.fname = self.fname_temp.format('lmps')
        self.sst.write_lammps(self.fname)
        self.assertTrue(os.path.isfile(self.fname))
        self.assertGreater(os.path.getsize(self.fname), 0)

    def test_write_lammps_empty(self):
        for p in self.sst.particles:
            self.sst.particles.remove(p.tag, update=False)
        self.sst.remove_spare_bonding()
        self.fname = self.fname_temp.format('empty.lmps')
        self.sst.write_lammps(self.fname)
        self.assertTrue(os.path.isfile(self.fname))
        self.assertGreater(os.path.getsize(self.fname), 0)

    def test_write_lammps_mol(self):
        self.fname = self.fname_temp.format('lmps')
        self.sst.write_lammps_mol(self.fname)
        self.assertTrue(os.path.isfile(self.fname))
        self.assertGreater(os.path.getsize(self.fname), 0)

    def test_write_mol(self):
        self.fname = self.fname_temp.format('mol')
        self.sst.write_mol(self.fname)
        self.assertTrue(os.path.isfile(self.fname))
        self.assertGreater(os.path.getsize(self.fname), 0)

    def test_write_pdb(self):
        self.fname = self.fname_temp.format('pdb')
        self.sst.write_pdb(self.fname)
        self.assertTrue(os.path.isfile(self.fname))
        self.assertGreater(os.path.getsize(self.fname), 0)

    def test_write_cssr(self):
        self.fname = self.fname_temp.format('cssr')
        self.sst.write_cssr(self.fname)
        self.assertTrue(os.path.isfile(self.fname))
        self.assertGreater(os.path.getsize(self.fname), 0)

    def tearDown(self):
        if len(self.fname) > 0:
            os.remove(self.fname)


if __name__ == '__main__':
    unittest.main()
