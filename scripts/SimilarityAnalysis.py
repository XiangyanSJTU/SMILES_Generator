#!/usr/bin/env python3
# coding=utf-8

import math
import sys

sys.path.append('..')
from app.models import *
from sqlalchemy import func

import argparse

parser = argparse.ArgumentParser(description='This is a code to generate similarity information')
parser.add_argument('--useChirality', default=False, help='use chirality to do the simularity analysis')
parser.add_argument('--clearinfo', default=False, help='set all similarity info False.')
parser.add_argument('--similarity', default=False, help='Generate Similarity information')
opt = parser.parse_args()
if opt.clearinfo:
    stereos = session.query(StereoIsomer)
    n = 10000
    N = math.ceil(stereos.count() / n)
    for j in range(N):
        stereos = session.query(StereoIsomer).filter(StereoIsomer.id >= j * n).filter(StereoIsomer.id < (j + 1) * n)
        for i, stereo in enumerate(stereos):
            if i % 1000 == 0:
                sys.stdout.write('\r %i / %i' % (i + j * n, N * n))
            stereo.training = False
            stereo.similarity_010 = False
            stereo.similarity_020 = False
            stereo.similarity_030 = False
            stereo.similarity_050 = False
            stereo.similarity_080 = False
            stereo.similarity_100 = False
            stereo.similarity_150 = False
            stereo.similarity_200 = False
        session.commit()
elif opt.similarity:
    if opt.useChirality:
        stereos = session.query(StereoIsomer)
        n = 1000
        N = math.ceil(stereos.count() / n)
        # similarity_cutoff_list = [0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0]
        similarity_cutoff_list = [0.5]
        smiles_dict = {}  # {0.1: []}
        for cutoff in similarity_cutoff_list:
            smiles_dict[cutoff] = []
        for j in range(N):
            stereos = session.query(StereoIsomer).filter(StereoIsomer.id >= j * n).filter(StereoIsomer.id < (j + 1) * n)
            for i, stereo in enumerate(stereos):
                sys.stdout.write('\r %i / %i' % (i + j * n, N * n))
                for cutoff in similarity_cutoff_list:
                    if not is_similar(stereo.smiles, smiles_dict[cutoff], cutoff=cutoff):
                        smiles_dict[cutoff].append(stereo.smiles)
                        if cutoff == 0.1:
                            stereo.similarity_010 = True
                            session.commit()
                        elif cutoff == 0.2:
                            stereo.similarity_020 = True
                            session.commit()
                        elif cutoff == 0.3:
                            stereo.similarity_030 = True
                            session.commit()
                        elif cutoff == 0.5:
                            stereo.similarity_050 = True
                            session.commit()
                        elif cutoff == 0.8:
                            stereo.similarity_080 = True
                            session.commit()
                        elif cutoff == 1.0:
                            stereo.similarity_100 = True
                            session.commit()
                        elif cutoff == 1.5:
                            stereo.similarity_150 = True
                            session.commit()
                        elif cutoff == 2.0:
                            stereo.similarity_200 = True
                            session.commit()
                    else:
                        if cutoff == 0.1:
                            stereo.similarity_010 = False
                            session.commit()
                        elif cutoff == 0.2:
                            stereo.similarity_020 = False
                            session.commit()
                        elif cutoff == 0.3:
                            stereo.similarity_030 = False
                            session.commit()
                        elif cutoff == 0.5:
                            stereo.similarity_050 = False
                            session.commit()
                        elif cutoff == 0.8:
                            stereo.similarity_080 = False
                            session.commit()
                        elif cutoff == 1.0:
                            stereo.similarity_100 = False
                            session.commit()
                        elif cutoff == 1.5:
                            stereo.similarity_150 = False
                            session.commit()
                        elif cutoff == 2.0:
                            stereo.similarity_200 = False
                            session.commit()
    else:
        mols = session.query(Molecule)
        n = 1000
        N = math.ceil(mols.count() / n)
        similarity_cutoff_list = [0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0]
        # similarity_cutoff_list = [0.5]
        smiles_dict = {}  # {0.1: []}
        for cutoff in similarity_cutoff_list:
            smiles_dict[cutoff] = []
        for j in range(N):
            mols = session.query(Molecule).filter(Molecule.id >= j * n).filter(Molecule.id < (j + 1) * n)
            for i, mol in enumerate(mols):
                if i % 100 == 0:
                    sys.stdout.write('\r %i / %i' % (i + j * n, N * n))
                for cutoff in similarity_cutoff_list:
                    if not is_similar(mol.smiles, smiles_dict[cutoff], cutoff=cutoff):
                        smiles_dict[cutoff].append(mol.smiles)
                        stereo = mol.stereoisomers.first()
                        if cutoff == 0.1:
                            stereo.similarity_010 = True
                        elif cutoff == 0.2:
                            stereo.similarity_020 = True
                        elif cutoff == 0.3:
                            stereo.similarity_030 = True
                        elif cutoff == 0.5:
                            stereo.similarity_050 = True
                        elif cutoff == 0.8:
                            stereo.similarity_080 = True
                        elif cutoff == 1.0:
                            stereo.similarity_100 = True
                        elif cutoff == 1.5:
                            stereo.similarity_150 = True
                        elif cutoff == 2.0:
                            stereo.similarity_200 = True
                session.commit()

    f = open('similarity_analysis.txt', 'w')
    tasks = session.query(Task).filter(Task.n_heavy > 5).order_by(Task.n_heavy)
    for task in tasks:
        f.write('%i %s\n' % (task.n_heavy, ' '.join(list(map(str, task.get_mol_number_list(similarity=True))))))
else:
    f = open('similarity_analysis.txt', 'w')
    tasks = session.query(Task).filter(Task.n_heavy > 5).order_by(Task.n_heavy)
    for task in tasks:
        sys.stdout.write('\rtask id = %i\n' % task.id)
        f.write('%i %s\n' % (task.n_heavy, ' '.join(list(map(str, task.get_mol_number_list(similarity=False))))))
