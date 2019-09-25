#!/usr/bin/env python3
# coding=utf-8

import sys
import math
sys.path.append('..')
from app.models import *
import argparse

parser = argparse.ArgumentParser(description='This is a code to generate final fingerprint')
parser.add_argument('--maxPath', type=int, help='The maxPath of topological fingerprint', default=7)
parser.add_argument('--cutoff2', type=float, help='cutoff^2 prefactor', default=1)

opt = parser.parse_args()


def is_train(fp, fp_list, cutoff2):
    for _fp in fp_list:
        if (fp-_fp).dot(fp-_fp) < cutoff2:
            return False
    return True


def get_rdk_fp_dict():
    hash_list = []
    if os.path.exists('fp.idx'):
        for line in open('fp.idx', 'r').readlines():
            hash_list.append(int(line))
        hash_list.sort()
    else:
        tasks = session.query(Task).filter(Task.n_heavy <= opt.maxPath + 1)
        for task in tasks:
            for mol in task.molecules:
                for stereo in mol.stereoisomers:
                    topol_fp = json.loads(stereo.rdk_fp)
                    for hash in topol_fp.keys():
                        if int(hash) not in hash_list:
                            hash_list.append(int(hash))
        hash_list.sort()
        f = open('fp.idx', 'w')
        for hash in hash_list:
            f.write('%i\n' % hash)

    unfound_hash = [False] * len(hash_list)
    training_fp_list = []
    stereos = session.query(StereoIsomer)
    n = 1000
    N = math.ceil(stereos.count() / n)
    for j in range(N):
        stereos = session.query(StereoIsomer).filter(StereoIsomer.id >= j * n).filter(StereoIsomer.id < (j + 1) * n)
        for i, stereo in enumerate(stereos):
            sys.stdout.write('\r %i / %i' % (i + j * n, N * n))
            topol_fp = json.loads(stereo.rdk_fp)
            fp_list = []
            for i, hash in enumerate(hash_list):
                if topol_fp.get(str(hash)) is None:
                    fp_list.append(0)
                else:
                    fp_list.append(topol_fp.get(str(hash)))
                    if not unfound_hash[i]:
                        unfound_hash[i] = True
                        stereo.training = True

            stereo.final_fp = json.dumps(fp_list)
            fp = np.array(fp_list)
            if stereo.training:
                training_fp_list.append(fp)
            else:
                if is_train(fp, training_fp_list, cutoff2=opt.cutoff2 * 1.15**(stereo.molecule.task.n_heavy*2)):
                    training_fp_list.append(fp)
                    stereo.training = True
                else:
                    stereo.training = False
    session.commit()


get_rdk_fp_dict()
