#!/usr/bin/env python3
# coding=utf-8

import sys

sys.path.append('..')
from config import *

sys.path.append(Config.MS_TOOLS_DIR)
from create.alkane import *
from create.alkene import *
from create.alkyne import *
from create.cations import *
from create.anions import *
from create.hydrocarbon_ring import *
from create.benzene import *
from create.alcohol import *
from create.ketone import *
from create.ether import *
from create.ammounia import *
from app.models import *

import argparse

parser = argparse.ArgumentParser(description='This is a code to generate smiles list and save all infomation in \
                                 Task table in smiles.sqlite')
parser.add_argument('-n', '--nheavy', type=int, help='The cutoff of heavy atoms', default=0)
parser.add_argument('--nheavy_nonCarbon', type=int, help='The cutoff of non-carbon heavy atoms', default=0)
opt = parser.parse_args()


def generate_ion(n_heavy):
    create_imidazolium(n_heavy)
    create_piperidinium(n_heavy)
    create_pyrrolidinium(n_heavy)
    create_amide_0h(n_heavy)
    create_amide_1h(n_heavy)
    create_amide_2h(n_heavy)
    create_pyridinium(n_heavy)
    create_phosphonium(n_heavy)
    create_AlkylAcetate(n_heavy)
    create_FluoroAlkylAcetate(n_heavy)
    create_AlkylSulfonate(n_heavy)
    create_FluoroAlkylSulfonate(n_heavy)
    create_AlkylSulfonylImide(n_heavy)
    create_FluoroAlkylSulfonylImide(n_heavy)


def generate_CH(n_heavy):
    create_alkane(n_heavy)
    create_alkene_1(n_heavy)
    create_alkene_2(n_heavy)
    create_alkene_3(n_heavy)

    create_alkyne_1(n_heavy)
    create_alkyne_2(n_heavy)
    create_alkyne_3(n_heavy)
    create_alkene_1_alkyne_1(n_heavy)
    create_alkene_1_alkyne_2(n_heavy)
    create_alkene_2_alkyne_1(n_heavy)

    create_3ring_1(n_heavy)
    create_4ring_1(n_heavy)
    create_5ring_1(n_heavy)
    create_6ring_1(n_heavy)
    create_7ring_1(n_heavy)
    create_8ring_1(n_heavy)
    create_9ring_1(n_heavy)
    create_10ring_1(n_heavy)
    create_11ring_1(n_heavy)
    create_12ring_1(n_heavy)
    create_13ring_1(n_heavy)
    create_benzene_1(n_heavy)
    create_naphthalene_1(n_heavy)

    create_ring_ring(n_heavy)
    create_ring_ring_1(n_heavy)
    create_ring_ring_2(n_heavy)
    create_ring_alkene_1(n_heavy)
    create_ring_alkene_2(n_heavy)


if __name__ == '__main__':
    if not os.path.exists(os.path.join('..', 'database', 'smiles.sqlite')):
        create_table()
    n_heavy = opt.nheavy
    nheavy_nonCarbon = opt.nheavy_nonCarbon
    if 'alkane' in Config.classes:
        create_alkane(n_heavy)
    if 'alcohol' in Config.classes:
        for i in range(1, nheavy_nonCarbon + 1):
            create_alcohol(n_heavy, i)
    if 'ketone' in Config.classes:
        for i in range(1, nheavy_nonCarbon + 1):
            create_ketone(n_heavy, i)
    if 'ether' in Config.classes:
        for i in range(1, nheavy_nonCarbon + 1):
            create_ether(n_heavy, i)
    if 'ammounia' in Config.classes:
        for i in range(1, nheavy_nonCarbon + 1):
            create_ammounia(n_heavy, i)
    # generate_CH()
    # generate_ion(opt.nheavy)
