import requests
import re
from StringIO import StringIO
from pysimm import system, lmps, forcefield

try:
    import pandas as pd
except ImportError:
    pd = None

# Check whether the pandas installed or not
if not pd:
    print('The script requires pandas to be installed. Exiting...')
    exit(1)

# Requesting the XYZ of the unit MOF cell from the web-resource
resp = requests.get('https://raw.githubusercontent.com/WMD-group/BTW-FF/master/structures/IRMOF-14.xyz')
xyz = StringIO(resp.text)

# Parsing the text stream to form the xyz-like pandas table
df = pd.read_table(xyz, sep='\s+', names=['tag', 'type', 'x', 'y', 'z'], usecols=[0, 1, 2, 3, 4], skiprows=1)
# Retyping the atom names
df['type'] = df['type'].map(lambda vr: vr[0] if vr[0] != 'Z' else 'Zn')

# Writing XYZ
with file('irmof-14_clean.xyz', 'w') as f:
    f.write(str(len(df)) + '\nThis is the place for the header of your XYZ file\n')
    df[['type', 'x', 'y', 'z']].to_csv(f, sep='\t', header=False, index=False)

# Initial setup of the pysimm System with MOF
s = system.System()
tmp = resp.text.encode('ascii', 'ignore').split('\n')
for line in tmp[1:-1]:
    data = line.split()
    tag, ptype, x, y, z, restof = data[:6]
    elem = re.sub('\d+', '', ptype)
    bonds = map(int, data[6:])
    p = system.Particle(tag=int(tag), elem=elem, type_name=ptype, x=float(x), y=float(y), z=float(z), bonds=bonds)
    s.particles.add(p)

for p in s.particles:
    for pb in p.bonds:
        if p.tag < pb:
            s.bonds.add(system.Bond(a=p, b=s.particles[pb]))
s.add_particle_bonding()

# Assign Dreiding forcefield parameters to the atoms of the structure
f = forcefield.Dreiding()
o_3 = s.particle_types.add(f.particle_types.get('O_3')[0].copy())
o_r = s.particle_types.add(f.particle_types.get('O_R')[0].copy())
c_r = s.particle_types.add(f.particle_types.get('C_R')[0].copy())
zn = s.particle_types.add(f.particle_types.get('Zn')[0].copy())
h_ = s.particle_types.add(f.particle_types.get('H_')[0].copy())
for p in s.particles:
    if p.elem == 'O':
        if p.bonds.count == 4:
            p.type = o_3
        else:
            p.type = o_r
    if p.elem == 'C':
        p.type = c_r
    if p.elem == 'Zn':
        p.type = zn
    if p.elem == 'H':
        p.type = h_
f.assign_btypes(s)
f.assign_atypes(s)
f.assign_dtypes(s)
f.assign_itypes(s)

# Assign the calculation box size assuming it is cubic
cc_bnd_lngth = 1.363
dim = cc_bnd_lngth / 2 + max(df['x'].values) - min(df['x'].values)
s.dim = system.Dimension(dx=dim, dy=dim, dz=dim, center=[dim/2, dim/2, dim/2])
s.forcefield = 'dreiding-lj'
s.pair_style = 'lj'
s.bond_style = 'harmonic'
s.angle_style = 'harmonic'
s.dihedral_style = 'harmonic'
s.improper_style = 'harmonic'

s.write_lammps('irmof-14.lmps')

