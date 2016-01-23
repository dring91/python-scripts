# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 09:37:28 2015

@author: David
"""
"""
def updateVal(dic,key,item):
  
  # try to pop a key off the dict to modify the value
  # if no error, update the value
  # if error, then add the key and value pair
  try:
    val = dic.pop(key)
    val.update(item)
    dic.update({key:val})
  except KeyError:
    dic.update({key:set(item)})
  
  return dic
  
def makeNestedDict(atoms):
  
  # Nesting scheme
  # atomtype -> molecule_id -> coords -> [x,y,z]
  atomtype = {}
  mol_id = {}
  coords = {}
  for atom in atoms:
    atomtype = updateVal(atomtype,atom[2],atom[1:2])
    mol_id = updateVal(mol_id,atom[1],atom[0:1])
    coords.update({atom[0]:atom[3:]})
    
  return atomtype, mol_id, coords
  """
def updateItem(dict,key,val):

  try:
    item = dict.pop(key)
    item.extend(val)
  except KeyError:
    item = val
    
  return item
  
def updateVal(oDict,oKey,iKey,val):
  
  # try to pop a key off the dict to modify the value
  # if no error, update the value
  # if error, then add the key and value pair
  try:
    iDict = oDict.pop(oKey)
    item = updateItem(iDict,iKey,val)
    iDict.update({iKey:item})
  except KeyError:
    iDict = {iKey:val}
    
  oDict.update({oKey:iDict})
  
  return oDict
  
def makeNestedDict(atoms):
  
  # Nesting scheme
  # atomtype -> molecule_id -> coords -> [x,y,z]
  info = {}
  coords = {}
  for atom in atoms:
    info = updateVal(info,atom[2],atom[1],atom[0:1])
    coords.update({atom[0]:atom[3:]})
    
  return info, coords
 
def printDict(info,coords):
  
  for types,mols in info.iteritems():
    for mol,atoms in mols.iteritems():
      for atom in atoms:
        print types,mol,atom,coords[atom][0],coords[atom][1],coords[atom][2]

def iterNested(info, coords):
  # loop over each atomtype
  for atype,mols in info.iteritems():
    # loop over each molecule in an atomtype
    for mol,atoms in mols.iteritems():
      # for each atom, return the coordinates
      for i, atom in enumerate(atoms):
        yield atype, mol, i+1, coords[atom]
        
def iterMol(info, coords):
  # loop over each atomtype
  for atype,mols in info.iteritems():
    # loop over each molecule in an atomtype
    for mol,atoms in mols.iteritems():
      # for each atom, return the coordinates
      yield atype, mol
