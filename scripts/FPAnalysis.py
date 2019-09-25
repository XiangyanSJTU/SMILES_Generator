#!/usr/bin/env python3
# coding=utf-8

import sys

sys.path.append('..')
from app.models import *
import argparse, math

parser = argparse.ArgumentParser(description='This is a code to generate fingerprint information in database')
parser.add_argument('-n', '--nheavy', type=int, default=0, help='Only analyze molecules with given heavy atoms')
opt = parser.parse_args()


def get_n_mols(mols, nheavy):
    count = 0
    for mol in mols:
        if mol.task.n_heavy == nheavy:
            count += 1
    return count


f = open('fp_analysis.txt', 'w')
stereos = session.query(StereoIsomer).limit(100000)
n = 10000
N = math.ceil(stereos.count() / n)
if Config.MORGAN:
    morgan_identifier_list = []
    morgan_fp_list = []
    morgan_identifier_count_dict = {}
    morgan_identifier_count_list = []
    for j in range(N):
        stereos = session.query(StereoIsomer).filter(StereoIsomer.id >= j * n).filter(StereoIsomer.id < (j + 1) * n)
        for i, stereo in enumerate(stereos):
            if opt.nheavy != 0 and stereo.molecule.task.n_heavy != opt.nheavy:
                continue
            if stereo.molecule.task.n_heavy < 6:
                continue
            sys.stdout.write('\rmorgan %i / %i' % (i + j * n, N * n))
            morgan_fp = json.loads(stereo.morgan_fp)
            if morgan_fp not in morgan_fp_list:
                morgan_fp_list.append(morgan_fp)
            for identifier in morgan_fp.keys():
                if identifier not in morgan_identifier_list:
                    morgan_identifier_list.append(identifier)
                if morgan_identifier_count_dict.get(identifier) is None:
                    morgan_identifier_count_dict[identifier] = 1
                else:
                    morgan_identifier_count_dict[identifier] += 1
    for key in morgan_identifier_count_dict.keys():
        morgan_identifier_count_list.append([key, morgan_identifier_count_dict[key]])
    morgan_identifier_count_list.sort(key=lambda x: x[1])
    mols = session.query(Molecule)
    f.write('morgan_fp: %i distinct fp / %i distinct smiles (no chirality)\n' % (len(morgan_fp_list),
                                                                                 get_n_mols(mols, opt.nheavy)))
    f.write('morgan_fp: fp length = %i\n' % (len(morgan_identifier_list)))
    for identifier, count in morgan_identifier_count_list:
        f.write('morgan_fp: identifier = %s, count = %i\n' % (identifier, count))
if Config.TOPOLOGICAL:
    rdk_identifier_list = []
    rdk_fp_list = []
    rdk_identifier_count_dict = {}
    rdk_identifier_count_list = []
    for j in range(N):
        stereos = session.query(StereoIsomer).filter(StereoIsomer.id >= j * n).filter(StereoIsomer.id < (j + 1) * n)
        for i, stereo in enumerate(stereos):
            if opt.nheavy != 0 and stereo.molecule.task.n_heavy != opt.nheavy:
                continue
            if stereo.molecule.task.n_heavy < 6:
                continue
            sys.stdout.write('\rrdk %i / %i' % (i + j * n, N * n))
            rdk_fp = json.loads(stereo.rdk_fp)
            if rdk_fp not in rdk_fp_list:
                rdk_fp_list.append(rdk_fp)
            else:
                print(stereo.smiles)
            for identifier in rdk_fp.keys():
                if identifier not in rdk_identifier_list:
                    rdk_identifier_list.append(identifier)
                if rdk_identifier_count_dict.get(identifier) is None:
                    rdk_identifier_count_dict[identifier] = 1
                else:
                    rdk_identifier_count_dict[identifier] += 1
    for key in rdk_identifier_count_dict.keys():
        rdk_identifier_count_list.append([key, rdk_identifier_count_dict[key]])
    rdk_identifier_count_list.sort(key=lambda x: x[1])
    mols = session.query(Molecule)
    f.write('\nrdk_fp: %i distinct fp / %i distinct smiles (no chirality)\n' % (len(rdk_fp_list),
                                                                                get_n_mols(mols, opt.nheavy)))
    f.write('rdk_fp: fp length = %i\n' % (len(rdk_identifier_list)))
    for identifier, count in rdk_identifier_count_list:
        f.write('rdk_fp: identifier = %s, count = %i\n' % (identifier, count))
if Config.PAIR:
    pair_identifier_list = []
    pair_fp_list = []
    pair_identifier_count_dict = {}
    pair_identifier_count_list = []
    for j in range(N):
        stereos = session.query(StereoIsomer).filter(StereoIsomer.id >= j * n).filter(StereoIsomer.id < (j + 1) * n)
        for i, stereo in enumerate(stereos):
            if opt.nheavy != 0 and stereo.molecule.task.n_heavy != opt.nheavy:
                continue
            if stereo.molecule.task.n_heavy < 6:
                continue
            sys.stdout.write('\rpair %i / %i' % (i + j * n, N * n))
            pair_fp = json.loads(stereo.pair_fp)
            if pair_fp not in pair_fp_list:
                pair_fp_list.append(pair_fp)
            for identifier in pair_fp.keys():
                if identifier not in pair_identifier_list:
                    pair_identifier_list.append(identifier)
                if pair_identifier_count_dict.get(identifier) is None:
                    pair_identifier_count_dict[identifier] = 1
                else:
                    pair_identifier_count_dict[identifier] += 1
    for key in pair_identifier_count_dict.keys():
        pair_identifier_count_list.append([key, pair_identifier_count_dict[key]])
    pair_identifier_count_list.sort(key=lambda x: x[1])
    mols = session.query(Molecule)
    f.write('\npair_fp: %i distinct fp / %i distinct smiles (no chirality)\n' % (len(pair_fp_list),
                                                                                 get_n_mols(mols, opt.nheavy)))
    f.write('pair_fp: fp length = %i\n' % (len(pair_identifier_list)))
    for identifier, count in pair_identifier_count_list:
        f.write('pair_fp: identifier = %s, count = %i\n' % (identifier, count))
if Config.TORSION:
    torsion_identifier_list = []
    torsion_fp_list = []
    torsion_identifier_count_dict = {}
    torsion_identifier_count_list = []
    for j in range(N):
        stereos = session.query(StereoIsomer).filter(StereoIsomer.id >= j * n).filter(StereoIsomer.id < (j + 1) * n)
        for i, stereo in enumerate(stereos):
            if opt.nheavy != 0 and stereo.molecule.task.n_heavy != opt.nheavy:
                continue
            if stereo.molecule.task.n_heavy < 6:
                continue
            sys.stdout.write('\rtorsion %i / %i' % (i + j * n, N * n))
            torsion_fp = json.loads(stereo.torsion_fp)
            if torsion_fp not in torsion_fp_list:
                torsion_fp_list.append(torsion_fp)
            for identifier in torsion_fp.keys():
                if identifier not in torsion_identifier_list:
                    torsion_identifier_list.append(identifier)
                if torsion_identifier_count_dict.get(identifier) is None:
                    torsion_identifier_count_dict[identifier] = 1
                else:
                    torsion_identifier_count_dict[identifier] += 1
    for key in torsion_identifier_count_dict.keys():
        torsion_identifier_count_list.append([key, torsion_identifier_count_dict[key]])
    torsion_identifier_count_list.sort(key=lambda x: x[1])
    mols = session.query(Molecule)
    f.write('\ntorsion_fp: %i distinct fp / %i distinct smiles (no chirality)\n' % (len(torsion_fp_list),
                                                                                    get_n_mols(mols, opt.nheavy)))
    f.write('torsion_fp: fp length = %i\n' % (len(torsion_identifier_list)))
    for identifier, count in torsion_identifier_count_list:
        f.write('torsion_fp: identifier = %s, count = %i\n' % (identifier, count))
