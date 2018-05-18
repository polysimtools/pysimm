from pysimm import system, lmps, forcefield

def run(test=False):
  # use a smiles string to query the pubchem search database and read the mol file returned from the http request
  # if cannot get to internet, read local cached response from pubchem
  try:
    s = system.read_pubchem_smiles('c1=cc=cc=c1')
  except:
      import os
      s = system.read_mol(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir, 'c1=cc=cc=c1.mol'))
  
  # the resulting system (benzene) has alternating double bonds
  # we want pysimm to recognize the ring as aromatic, so we define each bond in the ring to be bond order 'A'
  for b in s.bonds:
    if b.a.elem=='C' and b.b.elem=='C':
      b.order='A'
  
  # the resulting system has sufficient information to type with a forcefield, here we will use the Dreiding force field
  # we will also determine partial charges using the gasteiger algorithm
  s.apply_forcefield(forcefield.Dreiding(), charges='gasteiger')
  
  # we'll perform a 2 step energy minimization using the steepest decent and conjugate gradient algorithms in LAMMPS
  lmps.quick_min(s, min_style='sd', name='min_sd')
  lmps.quick_min(s, min_style='cg', name='min_cg')
  
  # write a few different file formats
  s.write_xyz('benzene.xyz')
  s.write_yaml('benzene.yaml')
  s.write_lammps('benzene.lmps')
  s.write_chemdoodle_json('benzene.json')

if __name__ == '__main__':
    run()