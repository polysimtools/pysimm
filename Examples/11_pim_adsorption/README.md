Example 11: Simulation of binary mixture adsorption isotherm of CO<sub>2</sub>-CH<sub>4</sub> ideal gas mixture with CASSNADRA and PyIAST
=========================================================================================================================================  
by Alexander Demidov
    
### Requirements
 * [PyIAST](https://github.com/CorySimon/pyIAST) -- Python package for simulation of Ideal Adsorption Solution Theory. <br>
     Can be easily installed with pip. For other installation options please refer to [PyIAST documentation](https://pyiast.readthedocs.io/en/latest/).  
 * [Pandas](https://pandas.pydata.org/) -- requirement of PyIAST used for advanced work with tabulated data (the package also can be easily installed with pip).
 
 
 
### Summary

Example shows the automated setup of multiple GCMC simulations using CASSANDRA followed up by processing of the simulations 
results and gathering them into adsorption isotherms. The isotherms are build using abilities of PyIAST package, furthermore, 
PyIAST is used to calculate the adsorption isotherms of the 1-to-1 mixture of the gases.

### Initialization

Top of the file there are few variables that can change the behavior of the script: ```is_simulate_loadings``` will 
enable/disable the Monte Carlo simulations that is by far the most time-consuming part of this example. If MC simulations are 
switched off the code will use pre-calculated loadings from the ```loading_file```.
The ```graphing``` option will switch between script saving graphics to file or to the screen or skipping the graphics step.
  
```python
is_simulate_loadings = True
loadings_file = 'loadings.dat'

graphing = 'ToFile'
```

Next section is the description of adsorbate system (names of gases, their mole fractions, and values of chemical potential for each gas)

```python
gas_names = ['ch4', 'co2']
mol_frac = [0.5, 0.5]
chem_pots = [lambda x: 2.4153 * numpy.log(x) - 36.722,
             lambda x: 2.40 * numpy.log(x) - 40.701]
```

Finally, similar to previous PySIMM:CASSANDRA examples we set the simulation system: create PySIMM system for each 
gas adding them to ```gases``` list, and a system for the framework which is wrapped into CASSANDRA system ```css```. 

```python
gases = []
for gn in gas_names:
    gases.append(system.read_lammps(osp.join(data_dir, gn + '.lmps')))
    gases[-1].forcefield = 'trappe/amber'

frame = system.read_lammps('pim.lmps')
frame.forcefield = 'trappe/amber'

css = cassandra.Cassandra(frame)
sim_settings = css.read_input('run_props.inp')
```


### Loadings simulation
The adsorption of 2 different gases will be calculated at 6 different pressures (0.1, 1, 5, 10, 25, and 50 bars) 
and will be added to the ```loadings``` dictionary with the key corresponding to the gas name. 

```python
gas_press = [0.1, 1, 5, 10, 25, 50]
loadings = dict.fromkeys(gas_names)
``` 
 
The loadings will be either simulated with CASSANDRA that are fully automated with PySIMM, or loaded from a file 
that contains pre-simulated values. Either way obtained loadings will be wrapped into the PyIAST ```ModelIsotherm``` class 
that describes simulated loadings using Bernauer-Emmet-Taller (BET) adsorption model.

```python
isotherms.append(pyiast.ModelIsotherm(pd.DataFrame(zip(gas_press, loadings[gn]), columns=[pk, lk]), 
                                      loading_key=lk, pressure_key=pk, model='BET', optimization_method='L-BFGS-B'))
```
 
### Binary mixture adsorption isotherm calculation
 
The PyIAST further in this example can be used for calculation of the binary mixture adsorption isotherm. 
For this type of calculations, PyIAST needs adsorption isotherms of pure gases and series of partial gas pressures.

```python
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
```

Please note, that for stable convergence of nonlinear solver in PyIAST core it is better to set good initial guesses 
of adsorbed mole fractions (e.g. if gas A is known to be much less preferable for adsorption than gas B it is logic 
to set initial guess of mole fractions to [0.05, 0.95] for A and B correspondingly). 
If the initial guess is and nonlinear solver did not converge, the PyIAST will throw an error. 
To avoid those possible crushes the try-catch block is set to iterate over the range of initial guesses, and use the ones that had converged.

### Loadings simulation
The example will print out the figure with simulated loadings for both pure gases, and for the 1-to-1 gas mixture 
as well as the approximation of all adsorption by the BET model.
The figure will be printed to the file if ```graphing``` option at the beginning of the script is set to `'ToFile'`, or will be shown on the screen
 if ```graphing``` set to `ToScreen`, otherwise the graphing part will be omitted.

 
 
 
 
 
 
 
 