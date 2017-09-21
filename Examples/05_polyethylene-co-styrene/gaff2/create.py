from pysimm import system, lmps, forcefield
from pysimm.apps.random_walk import copolymer
from pysimm.models.monomers.gaff2.pe import monomer as pe_monomer
from pysimm.models.monomers.gaff2.ps import monomer as ps_monomer

def run(test=False):
    # we'll make a polyethylene monomer and a polystyrene monomer from the pysimm models database
    pe = pe_monomer()
    ps = ps_monomer()
    
    # we'll instantiate a GAFF2 forcefield object for use later
    f = forcefield.Gaff2()
    
    # the monomers do not have any charges, so we will derive partial charges using the gasteiger algorithm
    pe.apply_charges(f, charges='gasteiger')
    ps.apply_charges(f, charges='gasteiger')
    
    # run the copolymer random walk method with 10 total repeat units, using an alternating pattern
    polymer = copolymer([pe, ps], 10, pattern=[1, 1], forcefield=f)
    
    # write a few different file formats
    polymer.write_xyz('polymer.xyz')
    polymer.write_yaml('polymer.yaml')
    polymer.write_lammps('polymer.lmps')
    polymer.write_chemdoodle_json('polymer.json')
    
if __name__ == '__main__':
    run()