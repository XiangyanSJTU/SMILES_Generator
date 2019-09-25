from app.models import *
import json

import argparse
parser = argparse.ArgumentParser(description='This is a code to dump smiles list in a txt file')
parser.add_argument('-c', '--category', type=str, help='The category of molecule to be dump')
parser.add_argument('-n', '--nheavy', type=int, help='The cutoff of heavy atoms', default=0)
opt = parser.parse_args()

mols_all = session.query(Task)
if opt.category =='all':
    mols = mols_all
else:
    mols = mols_all.filter(Task.category == opt.category)

if opt.nheavy!=0:
    mols = mols.filter(Task.n_heavy <= opt.nheavy)

f = open('smiles.txt', 'w')
for mol in mols:
    smiles_list = json.loads(mol.smiles_list)
    for smiles in smiles_list:
        f.write('%s\n' % (smiles))
