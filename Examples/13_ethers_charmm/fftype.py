from pysimm import lmps
from pysimm import system
from pysimm import forcefield

sst = system.read_mol('ethylpropylether.mol')

print(len(sst.particles))

for b in sst.bonds:
    b.order = 1

sst.apply_forcefield(forcefield.Charmm(), charges='gasteiger')

sst.write_lammps('test_out.lmps')
