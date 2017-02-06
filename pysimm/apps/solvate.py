import os
import json
import math
from pysimm import calc
from pysimm import system, lmps, forcefield


def solvate(s, **kwargs):
    s.wrap()
    wat = water_box()
    water_mol = water()
    buffer_dist = kwargs.get('buffer_dist') if kwargs.get('buffer_dist') is not None else 12.0
    closeness = kwargs.get('closeness') if kwargs.get('closeness') is not None else 1.0
    s.set_box(buffer_dist)
    
    print('adding waters to box with dimension:')
    print(s.dim.xlo, s.dim.xhi)
    print(s.dim.ylo, s.dim.yhi)
    print(s.dim.zlo, s.dim.zhi)
    
    empty_wat = system.System()
    empty_wat.dim = s.dim.copy()
    
    n_wat_x = int(math.ceil(s.dim.xhi/wat.dim.dx - 0.5))
    n_wat_y = int(math.ceil(s.dim.yhi/wat.dim.dy - 0.5))
    n_wat_z = int(math.ceil(s.dim.zhi/wat.dim.dz - 0.5))
    
    n_water_added = 0
    for nx in range(-n_wat_x, n_wat_x+1):
        for ny in range(-n_wat_y, n_wat_y+1):
            for nz in range(-n_wat_z, n_wat_z+1):
                w = wat.copy()
                for p in w.particles:
                    p.x += nx*w.dim.dx
                    p.y += ny*w.dim.dy
                    p.z += nz*w.dim.dz
                for m in w.molecules:
                    too_close = False
                    outside = False
                    for sp in s.particles:
                        for wp in m.particles:
                            if (wp.x > empty_wat.dim.xhi - 0.2 or wp.x < empty_wat.dim.xlo + 0.2 or 
                                wp.y > empty_wat.dim.yhi - 0.2 or wp.y < empty_wat.dim.ylo + 0.2 or
                                wp.z > empty_wat.dim.zhi - 0.2 or wp.z < empty_wat.dim.zlo + 0.2):
                                outside = True
                                break
                            if calc.distance(sp, wp) <= closeness:
                                too_close = True
                                break
                        if too_close or outside:
                            break
                    if not too_close and not outside:
                        new_w = water_mol.copy()
                        for p, np in zip(m.particles, new_w.particles):
                            np.set(x=p.x, y=p.y, z=p.z)
                        empty_wat.add(new_w)
                        n_water_added += 1
            
            print ' waters added so far ......', n_water_added
                    
    s.add(empty_wat, change_dim=False)
    if buffer_dist:
        s.set_box(0.5)

def water_box():
    wat = water()
    dim = 18.774349
    box = system.System()
    box.dim = system.Dimension(center=True, dx=dim, dy=dim, dz=dim)
    cfile = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                         os.pardir, os.pardir, 'dat', 'water', 'water_coords.json')
    with file(cfile) as f:
        water_coords = json.loads(f.read())

    c = water_coords['216box']        
    for m in range(216):
        w = wat.copy()
        for p in w.particles:
            x, y, z = c.pop(0)
            p.set(x=x, y=y, z=z)
        box.add(w, change_dim=False)
    
    return box
            

def water():
    s = system.read_pubchem_smiles('O')
    s.apply_forcefield(forcefield.Tip3p())
    
    return s


