from pysimm import cassandra
from pysimm import system
from os import path as osp
import numpy
import re
from matplotlib import pyplot as mplp

try:
    import pyiast
    import pandas
except ImportError:
    print('Either PyIAST or Pandas (that is PyIAST dependence) packages are not installed or cannot be found by this '
          'Python setup. \nThis example will not work properly. \nExiting...')
    exit(1)


def run(test=False):
    # Set to False if you **do not** want to recalculate pure gas adsorption isotherms
    is_simulate_loadings = False
    loadings_file = 'loadings.dat'

    # Option to draw the isotherms: it is either **'ToFile'** or **'ToScreen'** (case insensitive).
    # Any other value will be interpreted as no graphics
    graphing = 'none'

    # Gas names as they will be referred through simulations
    gas_names = ['ch4', 'co2']

    # Corresponding mole fractions of gases that will allow us to calculate their partial pressures through the Dalton's law
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

    # Constant for loadings calculations
    molec2mmols_g = 1e+3 / frame.mass

    # Setup of the GCMC simulations
    css = cassandra.Cassandra(frame)
    sim_settings = css.read_input('run_props.inp')

    # This function in given context will calculate the loading from short GCMC simulations
    def calculate_isotherm_point(gas_name, press):
        run_fldr = osp.join(gas_name, str(press))
        idx = gas_names.index(gas_name)
        # sim_settings.update({'Run_Name':  'gcmc'})
        css.add_gcmc(species=gases[idx], is_new=True, chem_pot=chem_pots[idx](press),
                     out_folder=run_fldr, props_file='gcmc.inp', **sim_settings)
        css.run()
        full_prp = css.run_queue[0].get_prp()
        return molec2mmols_g * numpy.average(full_prp[3][int(len(2 * full_prp[3]) / 3):])

    # This function in given context will load the pre-calculated loading value from previously done GCMC simulations
    def load_isotherm_point(gas_name, press):
        with open(loadings_file, 'r') as pntr:
            stream = pntr.read()
            tmp = stream.split('\n' + gas_name)[1]
            idx = re.search('[a-zA-Z]|\Z', tmp)
            value = re.findall('\n{:}\s+\d+\.\d+'.format(press), tmp[:idx.start()])[0]
            return float(re.split('\s+', value)[-1])

    # Calculation of adsorption isotherms for pure CH4 and CO2 gases for further usage in IAST simulations.
    # This is the **MOST TIME CONSUMING STEP** in this example, if you want to skip it switch the key is_simulated to False
    # The IAST will be done using PyIAST package, thus isotherms are wrapped into the corresponding object
    gas_press = [0.1, 1, 5, 10, 25, 50]
    lk = 'Loading(mmol/g)'
    pk = 'Pressure(bar)'
    isotherms = []
    loadings = dict.fromkeys(gas_names)
    for gn in gas_names:
        loadings[gn] = []
        for p in gas_press:
            if is_simulate_loadings:
                data = calculate_isotherm_point(gn, p)
            else:
                data = load_isotherm_point(gn, p)
            loadings[gn].append(data)
        isotherms.append(pyiast.ModelIsotherm(pandas.DataFrame(zip(gas_press, loadings[gn]), columns=[pk, lk]),
                                              loading_key=lk, pressure_key=pk, model='BET', optimization_method='Powell'))

    # The PyIAST run for calculating of mixed adsorption isotherm
    # Initial guesses of adsorbed mole fractions do span broad range of values, because PyIAST might not find
    #  solution at certain values of mole fractions and through an exception
    guesses = [[a, 1 - a] for a in numpy.linspace(0.01, 0.99, 50)]
    for in_g in guesses:
        mix_loadings = []
        try:
            for p in gas_press:
                mix_loadings.append(list(pyiast.iast(p * numpy.array(mol_frac), isotherms,
                                                     verboseflag=False, adsorbed_mole_fraction_guess=in_g)))
            mix_loadings = numpy.array(mix_loadings)
            break
        except:
            print('Initial guess {:} had failed to converge'.format(in_g))
            continue

    mix_loadings = numpy.sum(mix_loadings, axis=1)
    mix_isotherm = pyiast.ModelIsotherm(pandas.DataFrame(zip(gas_press, mix_loadings), columns=[pk, lk]),
                                        loading_key=lk, pressure_key=pk, model='BET', optimization_method='L-BFGS-B')


    # Output: Graphing of constructed isotherms
    def _plot_isotherms(ax, loc_gp, loc_isoth, loc_mix_load, loc_mix_isoth):
        rng = numpy.linspace(min(loc_gp), max(loc_gp), 100)
        ax.plot(loc_gp, loadings[gas_names[0]], 'og', lw=2.5, label='{:} loadings'.format(gas_names[0].upper()))
        ax.plot(rng, [loc_isoth[0].loading(t) for t in rng], '--g', lw=2, label='BET fit of {:} loadings'.format(gas_names[0].upper()))
        ax.plot(loc_gp, loadings[gas_names[1]], 'or', lw=2.5, label='{:} loadings'.format(gas_names[1].upper()))
        ax.plot(rng,  [loc_isoth[1].loading(t) for t in rng], '--r', lw=2, label='BET fit of {:} loadings'.format(gas_names[1].upper()))
        ax.plot(loc_gp, loc_mix_load, 'ob', lw=2.5, label='1-to-1 mixture loadings')
        ax.plot(rng, [loc_mix_isoth.loading(t) for t in rng], '--b', lw=2,  label='BET fit of 1-to-1 mixture loadings')
        ax.set_xlabel('Gas pressure [bar]', fontsize=20)
        ax.set_ylabel('Loading [mmol / g]', fontsize=20)
        ax.tick_params(axis='both', labelsize=16)
        ax.grid(True)
        ax.legend(fontsize=16)
        mplp.tight_layout()

    if graphing.lower() == 'tofile':
        fig, axs = mplp.subplots(1, 1, figsize=(10, 5))
        _plot_isotherms(axs, gas_press, isotherms, mix_loadings, mix_isotherm)
        mplp.savefig('pim1_mix_adsorption.png', dpi=192)
    elif graphing.lower() == 'toscreen':
        mplp.figure()
        axs = mplp.gca()
        _plot_isotherms(axs, gas_press, isotherms, mix_loadings, mix_isotherm)
        mplp.show()

    with open('iast_loadings.dat', 'w') as pntr:
        pntr.write('{:}\t\t{:}\n'.format(pk, lk))
        pntr.write('{:}-{:} 1-to-1\n'.format(gas_names[0], gas_names[1]))
        for p, ml in zip(gas_press, mix_loadings):
            pntr.write('{:}\t\t{:}\n'.format(p, ml))

if __name__ == '__main__':
    run(False)

