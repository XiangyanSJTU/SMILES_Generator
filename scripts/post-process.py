from app.models import *
from smiles import *
import sys

import argparse
parser = argparse.ArgumentParser(description='This is a post-process code')
opt = parser.parse_args()


def generate_molecule_table_from_task():
    tasks_all = session.query(Task).filter(Task.charge == 0).filter(Task.id > 971)
    tasks = tasks_all
    for task in tasks:
        task.generate_Molecule_table()
        # task.generate_Molecule_table(smiles_unique_check=True)

def update_similarity010_information():
    mols = session.query(Molecule).filter(Molecule.remark == None).filter(Molecule.similarity_010 != None)
    smiles_list_010 = []
    for mol in mols:
        if mol.similarity_010 == 'training':
            smiles_list_010.append(mol.smiles)
    mols = session.query(Molecule).filter(Molecule.remark == None).filter(Molecule.similarity_010 == None).limit(100000)
    for i, mol in enumerate(mols):
        sys.stdout.write('\r %i / %i' % (i, mols.count()))
        smiles = mol.smiles
        if is_similar(smiles, smiles_list_010, 0.1):
            mol.similarity_010 = 'validation'
        else:
            mol.similarity_010 = 'training'
            smiles_list_010.append(smiles)
        if i % 100 == 0:
            session.commit()
    session.commit()
def update_similarity030_information():
    mols = session.query(Molecule).filter(Molecule.remark == None).filter(Molecule.similarity_030 != None)
    smiles_list_030 = []
    for mol in mols:
        if mol.similarity_030 == 'training':
            smiles_list_030.append(mol.smiles)
    mols = session.query(Molecule).filter(Molecule.remark == None).filter(Molecule.similarity_030 == None).limit(100000)
    for i, mol in enumerate(mols):
        sys.stdout.write('\r %i / %i' % (i, mols.count()))
        smiles = mol.smiles
        if is_similar(smiles, smiles_list_030, 0.3):
            mol.similarity_030 = 'validation'
        else:
            mol.similarity_030 = 'training'
            smiles_list_030.append(smiles)
        if i % 100 == 0:
            session.commit()
    session.commit()
def update_similarity050_information():
    mols = session.query(Molecule).filter(Molecule.remark == None).filter(Molecule.similarity_050 != None)
    smiles_list_050 = []
    for mol in mols:
        if mol.similarity_050 == 'training':
            smiles_list_050.append(mol.smiles)
    mols = session.query(Molecule).filter(Molecule.remark == None).filter(Molecule.similarity_050 == None).limit(100000)
    for i, mol in enumerate(mols):
        sys.stdout.write('\r %i / %i' % (i, mols.count()))
        smiles = mol.smiles
        if is_similar(smiles, smiles_list_050, 0.5):
            mol.similarity_050 = 'validation'
        else:
            mol.similarity_050 = 'training'
            smiles_list_050.append(smiles)
        if i % 100 == 0:
            session.commit()
    session.commit()
def update_similarity080_information():
    mols = session.query(Molecule).filter(Molecule.remark == None).filter(Molecule.similarity_080 != None)
    smiles_list_080 = []
    for mol in mols:
        if mol.similarity_080 == 'training':
            smiles_list_080.append(mol.smiles)
    mols = session.query(Molecule).filter(Molecule.remark == None).filter(Molecule.similarity_080 == None).limit(100000)
    for i, mol in enumerate(mols):
        sys.stdout.write('\r %i / %i' % (i, mols.count()))
        smiles = mol.smiles
        if is_similar(smiles, smiles_list_080, 0.8):
            mol.similarity_080 = 'validation'
        else:
            mol.similarity_080 = 'training'
            smiles_list_080.append(smiles)
        if i % 100 == 0:
            session.commit()
    session.commit()
def update_similarity100_information():
    mols = session.query(Molecule).filter(Molecule.remark == None).filter(Molecule.similarity_100 != None)
    smiles_list_100 = []
    for mol in mols:
        if mol.similarity_100 == 'training':
            smiles_list_100.append(mol.smiles)
    mols = session.query(Molecule).filter(Molecule.remark == None).filter(Molecule.similarity_100 == None).limit(100000)
    for i, mol in enumerate(mols):
        sys.stdout.write('\r %i / %i' % (i, mols.count()))
        smiles = mol.smiles
        if is_similar(smiles, smiles_list_100, 1.0):
            mol.similarity_100 = 'validation'
        else:
            mol.similarity_100 = 'training'
            smiles_list_100.append(smiles)
        if i % 100 == 0:
            session.commit()
    session.commit()
def update_similarity150_information():
    mols = session.query(Molecule).filter(Molecule.remark == None).filter(Molecule.similarity_150 != None)
    smiles_list_150 = []
    for mol in mols:
        if mol.similarity_150 == 'training':
            smiles_list_150.append(mol.smiles)
    mols = session.query(Molecule).filter(Molecule.remark == None).filter(Molecule.similarity_150 == None).limit(100000)
    for i, mol in enumerate(mols):
        sys.stdout.write('\r %i / %i' % (i, mols.count()))
        smiles = mol.smiles
        if is_similar(smiles, smiles_list_150, 1.5):
            mol.similarity_150 = 'validation'
        else:
            mol.similarity_150 = 'training'
            smiles_list_150.append(smiles)
        if i % 100 == 0:
            session.commit()
    session.commit()
def update_similarity200_information():
    mols = session.query(Molecule).filter(Molecule.remark == None).filter(Molecule.similarity_200 != None)
    smiles_list_200 = []
    for mol in mols:
        if mol.similarity_200 == 'training':
            smiles_list_200.append(mol.smiles)
    mols = session.query(Molecule).filter(Molecule.remark == None).filter(Molecule.similarity_200 == None).limit(100000)
    for i, mol in enumerate(mols):
        sys.stdout.write('\r %i / %i' % (i, mols.count()))
        smiles = mol.smiles
        if is_similar(smiles, smiles_list_200, 2.0):
            mol.similarity_200 = 'validation'
        else:
            mol.similarity_200 = 'training'
            smiles_list_200.append(smiles)
        if i % 100 == 0:
            session.commit()
    session.commit()

def get_training_info():
    mols = session.query(Molecule).filter(Molecule.remark == None).filter(Molecule.similarity_010 == 'training')
    print(mols.count())

generate_molecule_table_from_task()
'''

for i in range(60):
    get_training_info()
    update_similarity030_information()
    update_similarity050_information()
    update_similarity010_information()
    update_similarity080_information()
    update_similarity100_information()
    update_similarity150_information()
    update_similarity200_information()
    '''


