import json, sys, pandas as pd
from app.models import *
from smiles import *

import argparse
parser = argparse.ArgumentParser(description='This is a code to generate smiles list')
parser.add_argument('-i', '--input', type=str, help='The input smiles list, check the smiles is in the database or not')
opt = parser.parse_args()

if __name__ == '__main__':
    df = pd.read_csv(opt.input, sep='\s+', header=0)
    smiles_list = df.SMILES.unique().tolist()
    for i, smiles in enumerate(smiles_list):
        sys.stdout.write('\r %i / %i    ' % (i, len(smiles_list)))
        if len(get_stereo_isomer(smiles)) == 1:
            smiles = get_canonical_smiles(smiles)
            mols = session.query(Molecule).filter(Molecule.smiles == smiles)
            if mols.count() == 0 and get_heavy_atom_numbers(smiles) <= 15:
                print(smiles)
            elif mols.count() == 1:
                continue
            elif get_heavy_atom_numbers(smiles) <= 15:
                raise Exception('%s, count = %i, error' % (smiles, mols.count()))

