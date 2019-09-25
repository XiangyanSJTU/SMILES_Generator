from app.models import *
import json, sys
from rdkit.Chem import AllChem as Chem

def create_imidazolium(n_heavy):
    print('Create imidazolium SMILES\n')
    category = 'imidazolium-v1'
    n = n_heavy
    id = 0
    while id < n:
        tasks_all = session.query(Task)
        tasks = tasks_all.filter(Task.category==category).order_by(Task.n_heavy.desc())
        task = tasks.first()
        if tasks.count() == 0:
            new_info = Task()
            new_info.smiles_list = json.dumps([get_canonical_smiles('Cn1cc[n+](c1)C')])
            new_info.n_heavy = 7
            new_info.category = category
            new_info.molecule_number = 1
            session.add(new_info)
            session.commit()
            id =7
        elif task.n_heavy >= n:
            break
        elif tasks.count() + 6 == task.n_heavy:
            smiles_list = json.loads(task.smiles_list)
            new_smiles_list = []
            for smiles in smiles_list:
                rdk_mol = Chem.MolFromSmiles(smiles)
                for atom in rdk_mol.GetAtoms():
                    if atom.IsInRing():
                        continue
                    if atom.GetImplicitValence() == 0:
                        continue
                    new_mol = Chem.RWMol(rdk_mol)
                    new_id = new_mol.AddAtom(Chem.Atom(6))
                    new_mol.AddBond(atom.GetIdx(), new_id, Chem.BondType.SINGLE)
                    new_smiles = get_canonical_smiles(Chem.MolToSmiles(new_mol))
                    if new_smiles not in new_smiles_list:
                        new_smiles_list.append(new_smiles)
            new_info = Task()
            new_info.smiles_list = json.dumps(new_smiles_list)
            new_info.n_heavy = task.n_heavy + 1
            new_info.category = category
            new_info.molecule_number = len(new_smiles_list)
            if new_info.smiles_list != '[]':
                session.add(new_info)
            session.commit()
            id = new_info.n_heavy
        else:
            raise Exception('Bad infomation error\n')
        sys.stdout.write('\r present %i / %i' % (id, n))

def create_piperidinium(n_heavy):
    print('Create piperidinium SMILES\n')
    category = 'piperidinium-v1'
    n = n_heavy
    id = 0
    while id < n:
        tasks_all = session.query(Task)
        tasks = tasks_all.filter(Task.category==category).order_by(Task.n_heavy.desc())
        task = tasks.first()
        if tasks.count() == 0:
            new_info = Task()
            new_info.smiles_list = json.dumps([get_canonical_smiles('C1CCCC[N+](C)(C)1')])
            new_info.n_heavy = 8
            new_info.category = category
            new_info.molecule_number = 1
            session.add(new_info)
            session.commit()
            id = 8
        elif task.n_heavy >= n:
            break
        elif tasks.count() + 7 == task.n_heavy:
            smiles_list = json.loads(task.smiles_list)
            new_smiles_list = []
            for smiles in smiles_list:
                rdk_mol = Chem.MolFromSmiles(smiles)
                for atom in rdk_mol.GetAtoms():
                    if atom.IsInRing():
                        continue
                    if atom.GetImplicitValence() == 0:
                        continue
                    new_mol = Chem.RWMol(rdk_mol)
                    new_id = new_mol.AddAtom(Chem.Atom(6))
                    new_mol.AddBond(atom.GetIdx(), new_id, Chem.BondType.SINGLE)
                    new_smiles = get_canonical_smiles(Chem.MolToSmiles(new_mol))
                    if new_smiles not in new_smiles_list:
                        new_smiles_list.append(new_smiles)
            new_info = Task()
            new_info.smiles_list = json.dumps(new_smiles_list)
            new_info.n_heavy = task.n_heavy + 1
            new_info.category = category
            new_info.molecule_number = len(new_smiles_list)
            if new_info.smiles_list != '[]':
                session.add(new_info)
            session.commit()
            id = new_info.n_heavy
        else:
            raise Exception('Bad infomation error\n')
        sys.stdout.write('\r present %i / %i' % (id, n))

def create_pyrrolidinium(n_heavy):
    print('Create pyrrolidinium SMILES\n')
    category = 'pyrrolidinium-v1'
    n = n_heavy
    id = 0
    while id < n:
        tasks_all = session.query(Task)
        tasks = tasks_all.filter(Task.category==category).order_by(Task.n_heavy.desc())
        task = tasks.first()
        if tasks.count() == 0:
            new_info = Task()
            new_info.smiles_list = json.dumps([get_canonical_smiles('C1CCC[N+](C)(C)1')])
            new_info.n_heavy = 7
            new_info.category = category
            new_info.molecule_number = 1
            session.add(new_info)
            session.commit()
            id = 7
        elif task.n_heavy >= n:
            break
        elif tasks.count() + 6 == task.n_heavy:
            smiles_list = json.loads(task.smiles_list)
            new_smiles_list = []
            for smiles in smiles_list:
                rdk_mol = Chem.MolFromSmiles(smiles)
                for atom in rdk_mol.GetAtoms():
                    if atom.IsInRing():
                        continue
                    if atom.GetImplicitValence() == 0:
                        continue
                    new_mol = Chem.RWMol(rdk_mol)
                    new_id = new_mol.AddAtom(Chem.Atom(6))
                    new_mol.AddBond(atom.GetIdx(), new_id, Chem.BondType.SINGLE)
                    new_smiles = get_canonical_smiles(Chem.MolToSmiles(new_mol))
                    if new_smiles not in new_smiles_list:
                        new_smiles_list.append(new_smiles)
            new_info = Task()
            new_info.smiles_list = json.dumps(new_smiles_list)
            new_info.n_heavy = task.n_heavy + 1
            new_info.category = category
            new_info.molecule_number = len(new_smiles_list)
            if new_info.smiles_list != '[]':
                session.add(new_info)
            session.commit()
            id = new_info.n_heavy
        else:
            raise Exception('Bad infomation error\n')
        sys.stdout.write('\r present %i / %i' % (id, n))

def create_amide_0h(n_heavy):
    print('Create amide SMILES\n')
    category = 'amide-0h'
    n = n_heavy
    id = 0
    N = Chem.Atom(7)
    N.SetFormalCharge(1)
    while id < n:
        tasks_all = session.query(Task)
        tasks = tasks_all.filter(Task.category == category).filter(Task.n_heavy == id + 1)
        if tasks.count() == 1:
            id += 1
            continue
        elif tasks.count() == 0:
            tasks = tasks_all.filter(Task.category == 'alkane').filter(Task.n_heavy == id + 1)
            task = tasks.first()
            smiles_list = json.loads(task.smiles_list)
            new_smiles_list = []
            for smiles in smiles_list:
                rdk_mol = Chem.MolFromSmiles(smiles)
                for atom in rdk_mol.GetAtoms():
                    if atom.GetImplicitValence()==0:
                        new_mol = Chem.RWMol(rdk_mol)
                        new_mol.ReplaceAtom(atom.GetIdx(), N)
                        new_smiles = get_canonical_smiles(Chem.MolToSmiles(new_mol))
                        if new_smiles not in new_smiles_list:
                            new_smiles_list.append(new_smiles)
            new_info = Task()
            new_info.smiles_list = json.dumps(new_smiles_list)
            new_info.n_heavy = id + 1
            new_info.category = category
            new_info.molecule_number = len(new_smiles_list)
            if new_info.smiles_list != '[]':
                session.add(new_info)
            session.commit()
            id += 1
        else:
            raise Exception('Bad infomation error\n')
        sys.stdout.write('\r present %i / %i' % (id, n))

def create_amide_1h(n_heavy):
    print('Create amide SMILES\n')
    category = 'amide-1h'
    n = n_heavy
    id = 0
    N = Chem.Atom(7)
    N.SetFormalCharge(1)
    while id < n:
        tasks_all = session.query(Task)
        tasks = tasks_all.filter(Task.category == category).filter(Task.n_heavy == id + 1)
        if tasks.count() == 1:
            id += 1
            continue
        elif tasks.count() == 0:
            tasks = tasks_all.filter(Task.category == 'alkane').filter(Task.n_heavy == id + 1)
            task = tasks.first()
            smiles_list = json.loads(task.smiles_list)
            new_smiles_list = []
            for smiles in smiles_list:
                rdk_mol = Chem.MolFromSmiles(smiles)
                for atom in rdk_mol.GetAtoms():
                    if atom.GetImplicitValence()==1:
                        new_mol = Chem.RWMol(rdk_mol)
                        new_mol.ReplaceAtom(atom.GetIdx(), N)
                        new_smiles = get_canonical_smiles(Chem.MolToSmiles(new_mol))
                        if new_smiles not in new_smiles_list:
                            new_smiles_list.append(new_smiles)
            new_info = Task()
            new_info.smiles_list = json.dumps(new_smiles_list)
            new_info.n_heavy = id + 1
            new_info.category = category
            new_info.molecule_number = len(new_smiles_list)
            if new_info.smiles_list != '[]':
                session.add(new_info)
            session.commit()
            id += 1
        else:
            raise Exception('Bad infomation error\n')
        sys.stdout.write('\r present %i / %i' % (id, n))

def create_amide_2h(n_heavy):
    print('Create amide SMILES\n')
    category = 'amide-2h'
    n = n_heavy
    id = 0
    N = Chem.Atom(7)
    N.SetFormalCharge(1)
    while id < n:
        tasks_all = session.query(Task)
        tasks = tasks_all.filter(Task.category == category).filter(Task.n_heavy == id + 1)
        if tasks.count() == 1:
            id += 1
            continue
        elif tasks.count() == 0:
            tasks = tasks_all.filter(Task.category == 'alkane').filter(Task.n_heavy == id + 1)
            task = tasks.first()
            smiles_list = json.loads(task.smiles_list)
            new_smiles_list = []
            for smiles in smiles_list:
                rdk_mol = Chem.MolFromSmiles(smiles)
                for atom in rdk_mol.GetAtoms():
                    if atom.GetImplicitValence()==2:
                        new_mol = Chem.RWMol(rdk_mol)
                        new_mol.ReplaceAtom(atom.GetIdx(), N)
                        new_smiles = get_canonical_smiles(Chem.MolToSmiles(new_mol))
                        if new_smiles not in new_smiles_list:
                            new_smiles_list.append(new_smiles)
            new_info = Task()
            new_info.smiles_list = json.dumps(new_smiles_list)
            new_info.n_heavy = id + 1
            new_info.category = category
            new_info.molecule_number = len(new_smiles_list)
            if new_info.smiles_list != '[]':
                session.add(new_info)
            session.commit()
            id += 1
        else:
            raise Exception('Bad infomation error\n')
        sys.stdout.write('\r present %i / %i' % (id, n))

def create_pyridinium(n_heavy):
    print('Create pyridinium SMILES\n')
    category = 'pyridinium-v1'
    n = n_heavy
    id = 0
    while id < n:
        tasks_all = session.query(Task)
        tasks = tasks_all.filter(Task.category==category).order_by(Task.n_heavy.desc())
        task = tasks.first()
        if tasks.count() == 0:
            new_info = Task()
            new_info.smiles_list = json.dumps([get_canonical_smiles('C[n+]1ccccc1')])
            new_info.n_heavy = 7
            new_info.category = category
            new_info.molecule_number = 1
            session.add(new_info)
            session.commit()
            id = 7
        elif task.n_heavy >= n:
            break
        elif tasks.count() + 6 == task.n_heavy:
            smiles_list = json.loads(task.smiles_list)
            new_smiles_list = []
            for smiles in smiles_list:
                rdk_mol = Chem.MolFromSmiles(smiles)
                for atom in rdk_mol.GetAtoms():
                    if atom.IsInRing():
                        continue
                    if atom.GetImplicitValence() == 0:
                        continue
                    new_mol = Chem.RWMol(rdk_mol)
                    new_id = new_mol.AddAtom(Chem.Atom(6))
                    new_mol.AddBond(atom.GetIdx(), new_id, Chem.BondType.SINGLE)
                    new_smiles = get_canonical_smiles(Chem.MolToSmiles(new_mol))
                    if new_smiles not in new_smiles_list:
                        new_smiles_list.append(new_smiles)
            new_info = Task()
            new_info.smiles_list = json.dumps(new_smiles_list)
            new_info.n_heavy = task.n_heavy + 1
            new_info.category = category
            new_info.molecule_number = len(new_smiles_list)
            if new_info.smiles_list != '[]':
                session.add(new_info)
            session.commit()
            id = new_info.n_heavy
        else:
            raise Exception('Bad infomation error\n')
        sys.stdout.write('\r present %i / %i' % (id, n))

def create_phosphonium(n_heavy):
    print('Create phosphonium SMILES\n')
    category = 'phosphonium'
    n = n_heavy
    id = 0
    P = Chem.Atom(15)
    P.SetFormalCharge(1)
    while id < n:
        tasks_all = session.query(Task)
        tasks = tasks_all.filter(Task.category == category).filter(Task.n_heavy == id + 1)
        if tasks.count() == 1:
            id += 1
            continue
        elif tasks.count() == 0:
            tasks = tasks_all.filter(Task.category == 'alkane').filter(Task.n_heavy == id + 1)
            task = tasks.first()
            smiles_list = json.loads(task.smiles_list)
            new_smiles_list = []
            for smiles in smiles_list:
                rdk_mol = Chem.MolFromSmiles(smiles)
                for atom in rdk_mol.GetAtoms():
                    if atom.GetImplicitValence()==0:
                        new_mol = Chem.RWMol(rdk_mol)
                        new_mol.ReplaceAtom(atom.GetIdx(), P)
                        new_smiles = get_canonical_smiles(Chem.MolToSmiles(new_mol))
                        if new_smiles not in new_smiles_list:
                            new_smiles_list.append(new_smiles)
            new_info = Task()
            new_info.smiles_list = json.dumps(new_smiles_list)
            new_info.n_heavy = id + 1
            new_info.category = category
            new_info.molecule_number = len(new_smiles_list)
            if new_info.smiles_list != '[]':
                session.add(new_info)
            session.commit()
            id += 1
        else:
            raise Exception('Bad infomation error\n')
        sys.stdout.write('\r present %i / %i' % (id, n))