Example 12: Tacticity control (vinyl-type polynorbornene)
========================================================================================================================  
by Sibo Lin
    
### Requirements.   
 * [NumPy](https://numpy.org/) -- Can be easily installed with pip. 
 * [SciPy](https://www.scipy.org/scipylib/index.html) -- Can be easily installed with pip. 
 
 
### Summary

This example shows the ability of random_walk_tacticity function to build syndiotactic and isotactic vinyl-type polynorbornene chains.

### Differences from standard random_walk function

Compared to the standard random_walk function, the random_walk_tacticity function requires a few more arguments. In this specific example:

```python
polymer = random_walk_tacticity(A, 10, forcefield=f, capped=True, tacticity='isotactic', sim='no')
```

The monomer (A) must be *capped* (meaning that the head and tail atoms of the monomer are bonded to sacrificial atoms 
that will be removed when adding the monomer to a growing oligomer).

The *tacticity* of the polymer is a float, between 0 and 1. The special keywords "isotatic" and "syndiotactic" are 
mapped to values of 1 and 0, respectively. A value of 0.4 would give 40% isotactic insertions (and 60% syndiotactic insertions).

For small monomers such as styrene, isotatic polymer chains can be built up by simply (a) copying the last monomer, (b) 
translating this copy, and (d) removing the capping atoms and defining a new bond between the copy and the previous monomer. 
However, with very large monomers (potentially sterically bulky styrenes), such an algorithm would result in hardcore overlaps. 
So a step must be inserted between (b) and (d) to avoid hardcore overlaps: (c) a *rotation* of 180 degrees is applied to the 
new monomer along its bond to the previous monomer. For small monomers or purely syndiotactic chains, rotation does 
not need to be specified (and will default to 0).

Syndiotactic insertions are performed by (a) copying the last monomer, (b) translating this copy, (c) reflecting this 
copy, and (d) removing the capping atoms and defining a new bond between the copy and the previous monomer. Not that to 
define the mirror plane for step (c), one of the atoms of the monomer must be defined as the "mirrorAtom".

If the "sim" argument is set to a string such as "no", then the chain is simply built up without any simulations. 
Short chains build up this way later be equilibrated with the equil module. Alternatively, if the "sim" argument is 
not specified in the call to random_walk_tacticity, then a simulation occurs as the chain grows. In this case, the 
optional "md_spacing" argument specifies how many monomers to insert between molecular dynamics relaxtion steps. 
For monomers that result in torturous paths (such as vinyl-type polynorbornenes), inserting new monomers may result in 
hardcore overlaps. For instance, the monomer could be generated threaded through in the middle of a pre-existing aromatic 
ring, in such a way that could not be later fixed by molecular dynamics. The "error_check" parameter, if set to True, 
will trigger a procedure to check for hardcore overlaps. If an overlap is found, 
the last monomer is artificially shrunken, and gradually re-expanded to normal scale while allowing the nearby atoms 
to get pushed out of the way, preventing the new monomer from being threaded through pre-existing aromatic rings. 
For this example monomer, growing a 50-mer with relaxation between every monomer insertion usually results in an 
invalid oligomer with hardcore overlaps, but the error_check procedure does allow such chains to be generated in most cases:

```python 
polymer = random_walk_tacticity(A, 50, forcefield=f, capped=True, tacticity='syndiotactic', error_check=True, rotation=180, md_spacing=1)
```


