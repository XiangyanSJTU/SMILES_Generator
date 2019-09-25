from create.public import *

def create_benzene_1(n_heavy):
    growing_from_smiles(n_heavy, 'c1ccccc1', category='benzene')

def create_naphthalene_1(n_heavy):
    growing_from_smiles(n_heavy, 'c1ccc2c(c1)cccc2', category='naphthalene')


