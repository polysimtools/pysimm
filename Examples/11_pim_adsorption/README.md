Example 11: Binary mixture adsorption isotherm of CO<sub>2</sub>-CH<sub>4</sub> gas mixture with CASSANDRA and PyIAST
=========================================================================================================================================  
by Alexander Demidov
    
### Requirements
 * Note: in case you work with multiple Python versions, all additional Python packages (particularly PyIAST, and Pandas) should be installed for Python 2.7+.   
 * [PyIAST](https://github.com/CorySimon/pyIAST) -- Python package for simulation of Ideal Adsorption Solution Theory. <br>
     Can be easily installed with pip. For other installation options please refer to [PyIAST documentation](https://pyiast.readthedocs.io/en/latest/).  
 * [Pandas](https://pandas.pydata.org/) -- requirement of PyIAST used for advanced work with tabulated data (the package also can be easily installed with pip).
 
 
 
### Summary

This example shows the automated setup of multiple GCMC simulations using CASSANDRA followed up by gathering simulated loadings 
into adsorption isotherms. Calculated loadings are fitted with Bernauer-Emmet-Taller (BET) adsorption model using PyIAST package. 
Furthermore, PyIAST is used to calculate the adsorption isotherm of the 1-to-1 mixture of the gases for comparison with the simulation results.

### Initialization

At the top of the file there are a few variables that can be selected to change the execution of the script: ```is_simulate_loadings``` will 
enable/disable the Monte Carlo simulations that is by far the most time-consuming part of this example. If MC simulations are 
switched off the code will use pre-calculated gas loadings from the ```loadings_file```.
The ```graphing``` option will switch between script saving graphics to file or to the screen or skipping the graphics step.
  
```python
is_simulate_loadings = True
loadings_file = 'loadings.dat'

graphing = 'ToFile'
```

Next section is the description of the adsorbate system (names of gases, their mole fractions, and values of chemical potential for each gas)

```python
gas_names = ['ch4', 'co2']
mol_frac = [0.5, 0.5]
chem_pots = [lambda x: 2.4153 * numpy.log(x) - 36.722,
             lambda x: 2.40 * numpy.log(x) - 40.701]
```

Finally, similar to previous PySIMM:CASSANDRA examples we set the simulation system: create a PySIMM system for each 
gas adding them to the ```gases``` list, and a PySIMM system for the framework which is sent to constructor of CASSANDRA system ```css```. 

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


### Gas loadings simulation
The adsorption of 2 different gases will be calculated at 6 different pressures (0.1, 1, 5, 10, 25, and 50 bars) 
and will be added to the ```loadings``` dictionary with the key corresponding to the gas name. 

```python
gas_press = [0.1, 1, 5, 10, 25, 50]
loadings = dict.fromkeys(gas_names)
``` 
 
The loadings will be either simulated with CASSANDRA that are fully automated with PySIMM, or loaded from a file 
that contains pre-simulated values. Either way obtained loadings will be wrapped into the PyIAST ```ModelIsotherm``` class 
that describes fitted loadings using Bernauer-Emmet-Taller (BET) adsorption model.

```python
isotherms.append(pyiast.ModelIsotherm(pd.DataFrame(zip(gas_press, loadings[gn]), columns=[pk, lk]), 
                                      loading_key=lk, pressure_key=pk, model='BET', optimization_method='Powell'))
```
 
### Binary mixture adsorption IAST calculation
 
PyIAST can be used to calculate the "IAST" of the binary mixture adsorption isotherm. 
For this type of calculations, PyIAST requires adsorption isotherms of pure gases and series of partial gas pressures.

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
If the initial guess is inapropriate and nonlinear solver did not converge, the PyIAST will raise an error. 
The try-catch block is set to iterate over the range of initial guesses, and use the ones that had converged.

### Example output
The example will print out the figure with simulated gas loadings for both pure gases, and for the 1-to-1 gas mixture 
as well as the approximation of gas adsorption by the BET model.
The figure will be printed to the file if the ```graphing``` option at the beginning of the script is set 
to `'ToFile'`, or will be shown on the screen if ```graphing``` set to `ToScreen`, otherwise the graphing part will be omitted.

 
 
 
 
 
 
 
 