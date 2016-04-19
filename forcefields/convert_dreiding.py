#!/usr/bin/env python

##############################
#
#
#  forcefield.pcff Module
#
#
##############################

import sys
from xml.dom import minidom
from xml.etree import ElementTree as ET
from polysimtools.system import *
from polysimtools.forcefield import dreiding


def convert(file,out):
  ff=dreiding()
  f=open(file)
  line=f.readline()
  while line:
    line=line.split()
    if len(line)>0 and line[0]=='ATOMTYPES':
      line=f.readline()
      while line.split()[0]!='END':
        desc=line.split('!')[-1].strip()
        line=line.split()
        if line[1]!='X':
          try:
            ff.ptypes[line[0]].mass=float(line[2])
            ff.ptypes[line[0]].elem=line[1]
	    ff.ptypes[line[0]].desc=desc
          except:
            ff.ptypes[line[0]]=ParticleType(mass=float(line[2]),elem=line[1],name=line[0],desc=desc)
        line=f.readline()
    elif len(line)==1 and line[0]=='DIAGONAL_VDW':
      line=f.readline()
      while line.split()[0]!='END':
        line=line.split()
        try:
          ff.ptypes[line[0]].sigma=float(line[2])
          ff.ptypes[line[0]].epsilon=float(line[3])
        except:
          ff.ptypes[line[0]]=ParticleType(sigma=float(line[2]),epsilon=float(line[3]),name=line[0])
        line=f.readline()
    elif len(line)==1 and line[0]=='BOND_STRETCH':
      line=f.readline()
      while line.split()[0]!='END':
        line=line.split()
        if len(line)>3:
          name=line[0]+','+line[1]
          k=float(line[3])/2
          r0=float(line[4])
          try:
            ff.btypes[name].k=k
            ff.btypes[name].r0=r0
          except:
            ff.btypes[name]=BondType(k=k,r0=r0,name=name)
        line=f.readline()
    elif len(line)==1 and line[0]=='ANGLE_BEND':
      line=f.readline()
      while line.split()[0]!='END':
        line=line.split()
        if len(line)>4:
          name=line[0]+','+line[1]+','+line[2]
          k=float(line[4])/2
          theta0=float(line[5])
          try:
            ff.atypes[name].k=k
            ff.atypes[name].theta0=theta0
          except:
            ff.atypes[name]=AngleType(k=k,theta0=theta0,name=name)
        line=f.readline()
    elif len(line)==1 and line[0]=='TORSIONS':
      line=f.readline()
      while line.split()[0]!='END':
        line=line.split()
        if len(line)>3:
          name=line[0]+','+line[1]+','+line[2]+','+line[3]
          hyb_m1=2 if line[1][2]=='R' else int(line[1][2])
          hyb_m2=2 if line[2][2]=='R' else int(line[2][2])
          k=float(line[5])/2/hyb_m1/hyb_m2
          n=int(float(line[6]))
          d=-1*int(float(line[7]))
          try:
            ff.dtypes[name].k=k
	    ff.dtypes[name].n=n
	    ff.dtypes[name].d=d
          except:
            ff.dtypes[name]=DihedralType(k=k,n=n,d=d,name=name)
          if name=='X,C_2,C_3,X':
            print line
            print k,n,d
        line=f.readline()
    elif len(line)==1 and line[0]=='INVERSIONS':
      line=f.readline()
      while line.split()[0]!='END':
        if line[0]=='#': line=f.readline(); continue
        line=line.split()
        if len(line)>5:
          name=line[0]+','+line[1]+','+line[2]+','+line[3]
          k=float(line[5])/2
          x0=float(line[6])
          try:
            ff.itypes[name].k=k
            ff.itypes[name].x0=x0
          except:
            ff.itypes[name]=ImproperType(k=k,x0=x0,name=name)
        line=f.readline()
    else: line=f.readline()
  f.close()

  ff.write(out)

def assign_itypes(l,f=None):
  itype_names=set()
  for p in l.particles.values():
    name=l.particle_types[p.type].name.split('@')[-1]
    try:
      f.itypes[name+',X,X,X']
      itype_names.add(name+',X,X,X')
      bonded_to=[]
      for b in l.bonds.values():
        if b.a==p.tag: bonded_to.append(b.b)
	elif b.b==p.tag: bonded_to.append(b.a)
      l.impropers[l.nimpropers+1]=Improper(tag=l.nimpropers,type=0,a=p.tag,b=bonded_to[0],c=bonded_to[1],d=bonded_to[2])
      l.nimpropers+=1
    except: pass
  itype_names=list(itype_names)
  l.nimproper_types=len(itype_names)
  l.improper_style='harmonic'
  for i in range(l.nimproper_types):
    l.improper_types[i+1]=ImproperType(tag=i+1,name=itype_names[i],k=f.itypes[itype_names[i]].k,x0=f.itypes[itype_names[i]].x0)
  for i in l.impropers.values():
    for it in l.improper_types.values():
      if l.particle_types[l.particles[i.a].type].name==it.name.split(',')[0]:
        i.type=it.tag
  return l


if __name__=='__main__':
  convert(sys.argv[1],sys.argv[2])
