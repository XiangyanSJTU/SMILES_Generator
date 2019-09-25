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




