#!/usr/bin/env python3
# coding=utf-8

import sys

sys.path.append('..')

from app.models import *

import argparse

parser = argparse.ArgumentParser(description='This is a code to generate smiles list and its stereoisomers in Task and\
                                 Stereoisomer table in smiles.sqlite')
opt = parser.parse_args()


def generate_molecule_table_from_task():
    tasks_all = session.query(Task).filter(Task.charge == 0).filter(Task.n_heavy > 5)
    tasks = tasks_all
    for task in tasks:
        if task.molecules.count() != task.molecule_number:
            task.generate_Molecule_table(smiles_unique_check=True)


generate_molecule_table_from_task()
