Example 12: Tacticity control (syndiotactic polystyrene)
=========================================================================================================================================  
by Sibo Lin
    
### Requirements
 * Note: in case you work with multiple Python versions, all additional Python packages should be installed for Python 3+.   
 * [NumPy](https://numpy.org/) -- Can be easily installed with pip. 
 * [SciPy](https://www.scipy.org/scipylib/index.html) -- Can be easily installed with pip. 
 
 
 
### Summary

This example shows the ability of random_walk_tacticity function to build syndiotactic vinyl-type polynorbornene chains.

### Differences from standard random_walk function

Compared to the standard random_walk function, the random_walk_tacticity function requires a few more arguments. In this specific example:

`polymer = random_walk_tacticity(A, 10, forcefield=f,capped=True,tacticity='isotactic',sim='no')`

The monomer (A) must be *capped* (meaning that the head and tail atoms of the monomer are bonded to sacrificial atoms that will be removed when adding the monomer to a growing oligomer).

The *tacticity* of the polymer is a float, between 0 and 1. The special keywords "isotatic" and "syndiotactic" are mapped to values of 1 and 0, respectively. A value of 0.4 would give 40% syndiotactic insertions and 60% isotactic insertions. 

For small monomers such as styrene, isotatic polymer chains can be built up by simply (a) copying the last monomer, (b) translating this copy, and (d) removing the capping atoms and defining a new bond between the copy and the previous monomer. However, with very large monomers (potentially sterically bulky styrenes), such an algorithm would result in hardcore overlaps. So a step must be inserted between (b) and (d) to avoid hardcore overlaps: (c) a *rotation* of 180 degrees is applied to the new monomer along its bond to the previous monomer. For small monomers or purely syndiotactic chains, rotation does not need to be specified (and will default to 0).

Syndiotactic insertions are performed by (a) copying the last monomer, (b) translating this copy, (c) reflecting this copy, and (d) removing the capping atoms and defining a new bond between the copy and the previous monomer. Not that to define the mirror plane for step (c), one of the atoms of the monomer must be defined as the "mirrorAtom".
