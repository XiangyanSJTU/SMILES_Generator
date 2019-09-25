import json, sys, pandas as pd
from app.models import *
from smiles import *

import argparse
parser = argparse.ArgumentParser(description='This is a code to delete')
opt = parser.parse_args()

if __name__ == '__main__':
    tasks = session.query(Task).filter(Task.id ==838)
    for task in tasks:
        #task.delete()
        task.delete_molecule()
