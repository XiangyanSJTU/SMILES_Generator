#!/usr/bin/env python3
# coding=utf-8

import math
import sys

sys.path.append('..')

from app.models import *
from sqlalchemy import func
import argparse

parser = argparse.ArgumentParser(
    description='This is a code to generate smiles list used as training set in machine learning')
parser.add_argument('-t', '--type', type=str, help='The selection rule of training set')
parser.add_argument('--similarity', type=float, help='The cutoff value of FP-similarity rule', default=-1.0)

opt = parser.parse_args()

if opt.type == 'FP-similarity':
    if opt.similarity in [0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0]:
        f1 = open('train_%.2f.txt' % opt.similarity, 'w')
        f2 = open('validate_%.2f.txt' % opt.similarity, 'w')
        stereos = session.query(StereoIsomer)
        if opt.similarity == 0.1:
            sim = StereoIsomer.similarity_010
        elif opt.similarity == 0.2:
            sim = StereoIsomer.similarity_020
        elif opt.similarity == 0.3:
            sim = StereoIsomer.similarity_030
        elif opt.similarity == 0.5:
            sim = StereoIsomer.similarity_050
        elif opt.similarity == 0.8:
            sim = StereoIsomer.similarity_080
        elif opt.similarity == 1.0:
            sim = StereoIsomer.similarity_100
        elif opt.similarity == 1.5:
            sim = StereoIsomer.similarity_150
        elif opt.similarity == 2.0:
            sim = StereoIsomer.similarity_200
        else:
            sys.exit()
        stereos_train = stereos.filter(sim == True)
        for stereo in stereos_train:
            f1.write('%s\n' % stereo.smiles)
        stereos_validate = stereos.filter(sim == False).order_by(func.random()).limit(stereos_train.count())
        for stereo in stereos_validate:
            f2.write('%s\n' % stereo.smiles)
    else:
        print(opt.similarity)
elif opt.type == 'FP-distance':
    stereos = session.query(StereoIsomer)
    stereos_train = stereos.filter(StereoIsomer.training == True)
    file = open('training.txt', 'w')
    for stereo in stereos_train:
        file.write('%s\n' % stereo.smiles)
    file.close()
    file = open('validation.txt', 'w')
    for stereo in stereos.filter(StereoIsomer.training == False).order_by(func.random()).limit(stereos_train.count()):
        file.write('%s\n' % stereo.smiles)
