#!/usr/bin/env python3
# coding=utf-8

import pybel, json, sys, pandas as pd, os

os.chdir('..')
sys.path.append('.')

from smiles import *
from app.models import *
from rdkit.Chem import AllChem as Chem


import argparse
parser = argparse.ArgumentParser(description='This is a code to dump smiles list in a txt file')
parser.add_argument('-n', '--nheavy', type=int, help='The cutoff of heavy atoms', default=0)
parser.add_argument('-c', '--cutoff', type=float, help='The cutoff of similarity', default=0)
opt = parser.parse_args()


def similarity_comparison(smiles1, smiles2):
    from rdkit.Chem import AllChem as Chem
    from rdkit import DataStructs
    rdk_mol1 = Chem.MolFromSmiles(smiles1)
    fp1 = Chem.GetMorganFingerprintAsBitVect(rdk_mol1, 2)
    rdk_mol2 = Chem.MolFromSmiles(smiles2)
    fp2 = Chem.GetMorganFingerprintAsBitVect(rdk_mol2, 2)
    # print(smiles1, smiles2)
    return DataStructs.DiceSimilarity(fp1, fp2)

# a switch function of similarity_comparison
def similarity_score(smiles1, smiles2):
    a = similarity_comparison(smiles1, smiles2)
    if a < 0.6:
        return 0.
    else:
        return (a - 0.6) / 0.4

def get_similarity_score(smiles, smiles_list):
    score = 0.
    for s in smiles_list:
        score += similarity_score(smiles, s)
    return score

def is_similar(smiles, smiles_list, cutoff):
    score = 0.
    for s in smiles_list:
        score += similarity_score(smiles, s)
        if score > cutoff:
            return True
    return False

cation_smiles_list = []
anion_smiles_list = ['[O-][Cl](=O)(=O)=O', 'F[P-](F)(F)(F)(F)F', '[Cl-]', '[S-]C#N', 'N#C[N-]C#N', 'N#C[C-](C#N)C#N',
                     'N#C[B-](C#N)(C#N)C#N', '[O-][N+](=O)[O-]', 'F[B-](F)(F)F']

cation_category = ['imidazolium-v1', 'piperidinium-v1', 'pyrrolidinium-v1', 'amide-0h', 'amide-1h', 'pyridinium-v1', 'phosphonium']
anion_category = ['FluoroAlkylAcetate', 'FluoroAlkylSulfonate', 'FluoroAlkylSulfonylImide']

tasks_all = session.query(Task)
tasks = tasks_all.filter(Task.category.in_(cation_category)).filter(Task.n_heavy <= 15)
for task in tasks:
    smiles_list = json.loads(task.smiles_list)
    cation_smiles_list += smiles_list

tasks = session.query(Task).filter(Task.category.in_(anion_category)).filter(Task.n_heavy <= opt.nheavy)
for task in tasks:
    smiles_list = json.loads(task.smiles_list)
    anion_smiles_list += smiles_list

os.chdir('IonicLiquid')
f = open('cation_%.2f.txt' % (opt.cutoff), 'w')
smiles_list = []
for i, smiles in enumerate(cation_smiles_list):
    sys.stdout.write('\r%i / %i / %i' % (len(smiles_list), i, len(cation_smiles_list)))
    if not is_similar(smiles, smiles_list, opt.cutoff):
        smiles_list.append(smiles)
        f.write('%s\n' % (smiles))
cation_smiles_list = smiles_list
f.close()

f = open('anion_%.2f.txt' % (opt.cutoff), 'w')
smiles_list = []
for i, smiles in enumerate(anion_smiles_list):
    sys.stdout.write('\r%i / %i / %i' % (len(smiles_list), i, len(anion_smiles_list)))
    if not is_similar(smiles, smiles_list, opt.cutoff):
        smiles_list.append(smiles)
        f.write('%s\n' % (smiles))
anion_smiles_list = smiles_list
f.close()

f = open('il_%.2f.txt' % (opt.cutoff), 'w')
smiles_list = []
for i, a in enumerate(cation_smiles_list):
    sys.stdout.write('\r%i / %i / %i' % (len(smiles_list), i, len(cation_smiles_list)))
    for b in anion_smiles_list:
        smiles = a + '.' + b
        if not is_similar(smiles, smiles_list, opt.cutoff):
            smiles_list.append(smiles)
            f.write('%s\n' % (smiles))
f.close()



