Example 12: Tacticity control (polystyrene)
=========================================================================================================================================  
by Sibo Lin
    
### Requirements.   
 * [NumPy](https://numpy.org/) -- Can be easily installed with pip. 
 * [SciPy](https://www.scipy.org/scipylib/index.html) -- Can be easily installed with pip. 
 
 
### Summary

This example shows the ability of random_walk_tacticity function to build syndiotactic, isotactic, and mixed tacticity polystyrene chains.

### Differences from standard random_walk function

Compared to the standard random_walk function, the random_walk_tacticity function requires a few more arguments. In this specific example:

```python
polymer = random_walk_tacticity(A, 10, forcefield=f, capped=True, tacticity='isotactic', sim='no')
```

The *tacticity* of the polymer is a float, between 0 and 1. The special keywords "isotatic" and "syndiotactic" are 
mapped to values of 1 and 0, respectively. A value of 0.4 would give 40% isotactic insertions (and 60% syndiotactic insertions).

For small monomers such as styrene, isotatic polymer chains can be built up by simply (a) copying the last monomer, (b) 
translating this copy, and (d) removing the capping atoms and defining a new bond between the copy and the previous 
monomer. However, with very large monomers (potentially sterically bulky styrenes), such an algorithm would result in 
hardcore overlaps. So a step must be inserted between (b) and (d) to avoid hardcore overlaps: (c) a *rotation* of 180 
degrees is applied to the new monomer along its bond to the previous monomer. For small monomers or purely syndiotactic 
chains, rotation does not need to be specified (and will default to 0).

Syndiotactic insertions are performed by (a) copying the last monomer, (b) translating this copy, (c) reflecting this 
copy, and (d) removing the capping atoms and defining a new bond between the copy and the previous monomer. 
Not that to define the mirror plane for step (c), one of the atoms of the monomer must be defined as the "mirrorAtom".

