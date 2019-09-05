from pysimm import cassandra
from pysimm import system
from os import path as osp
import numpy
import pyiast
import pandas as pd

# Gas names as they will be referred through simulations
gas_names = ['ch4', 'co2']

# Corresponding mole fractions of gases that will allow us to calculate their partial pressures through Dalton's law
mol_frac = [0.5, 0.5]

# Calibrated previously functional forms of chemical potentials of gases for GCMC simulations as a functions of pressure
chem_pots = [lambda x: 2.4153 * numpy.log(x) - 36.722,
             lambda x: 2.40 * numpy.log(x) - 40.701]

# Root directory for some data (For PySIMM examples it is )
data_dir = osp.join('..', '09_cassandra_simulations', 'gcmc')

# Setup of adsorbate model
gases = []
for gn in gas_names:
    gases.append(system.read_lammps(osp.join(data_dir, gn + '.lmps')))
    gases[-1].forcefield = 'trappe/amber'

# Setup of adsorbent model
frame = system.read_lammps('pim.lmps')
frame.forcefield = 'trappe/amber'

# Setup of the GCMC simulations
css = cassandra.Cassandra(frame)
sim_settings = css.read_input('run_props.inp')

# This function in given context will
def calculate_isotherm_point(gas_name, press):
    run_fldr = osp.join(gas_name, str(press))
    idx = gas_names.index(gas_name)
    # sim_settings.update({'Run_Name':  'gcmc'})
    css.add_gcmc(species=gases[idx], is_new=True, chem_pot=chem_pots[idx](press),
                 out_folder=run_fldr, props_file='gcmc.inp', **sim_settings)
    css.run()
    full_prp = css.run_queue[0].get_prp()
    return numpy.average(full_prp[3][int(len(2 * full_prp[3]) / 3):])

# Calculation of adsorption isotherms for pure CH4 and CO2 gases for further usage in IAST simulations.
# This is the **MOST TIME CONSUMING STEP** in this example, if you want to skip it switch the key is_simulated to False
# The IAST will be done using PyIAST package, thus isotherms are wrapped into the corresponding object
gas_press = [0.1, 1, 5, 10, 25, 50]
lk = 'Loading(mmol/g)'
pk = 'Pressure(bar)'
isotherms = []
for gn in gas_names:
    loadings = []
    for p in gas_press:
        data = calculate_isotherm_point(gn, p)
        loadings.append(data)
    isotherms.append(pyiast.ModelIsotherm(pd.DataFrame(zip(gas_press, loadings),
                                          columns=[pk, lk]),  loading_key=lk, pressure_key=pk,
                                          model='BET', optimization_method='L-BFGS-B'))

# The PyIAST run for calculating of
# Initial guesses of adsorbed mole fractions do span broad range of values, because PyIAST might not find
#  solution at certain values of mole fractions and through an exception
guesses = [[a, 1 - a] for a in numpy.linspace(0.01, 0.99, 50)]
for in_g in guesses:
    mix_isotherm = []
    try:
        for p in gas_press:
            mix_isotherm.append(list(pyiast.iast(numpy.array([p] * 2), isotherms,
                                                 verboseflag=False, adsorbed_mole_fraction_guess=in_g)))
        mix_isotherm = numpy.array(mix_isotherm)
        break
    except:
        print('Initial guess {:} had failed to converge'.format(in_g))
        continue

print(mix_isotherm)


'''
prefix = gas_names[0] + '.prod'
prod_sim_settings = css.read_input('prod_props.inp')
prod_sim_settings.update({'Run_Name': prefix + '.gcmc',
                          'Start_Type': {'start_type': 'checkpoint', 'file_name': equil_prefix + '.gcmc.chk'}})
css.add_gcmc(species=gases[0], is_new=True, chem_pot=chem_pots[0],
             props_file=prefix + '.gcmc.inp', **prod_sim_settings)
css.run()
'''

