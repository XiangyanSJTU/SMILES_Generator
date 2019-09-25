from app.models import *
import json, sys
from rdkit.Chem import AllChem as Chem

def create_AlkylAcetate(n_heavy):
    print('Create AlkylAcetate SMILES\n')
    category = 'AlkylAcetate'
    n = n_heavy
    id = 0
    while id < n:
        tasks_all = session.query(Task)
        tasks = tasks_all.filter(Task.category==category).order_by(Task.n_heavy.desc())
        task = tasks.first()
        if tasks.count() == 0:
            new_info = Task()
            new_info.smiles_list = json.dumps([get_canonical_smiles('[O-]C(=O)C')])
            new_info.n_heavy = 4
            new_info.category = category
            new_info.molecule_number = 1
            session.add(new_info)
            session.commit()
            id = 4
        elif task.n_heavy >= n:
            break
        elif tasks.count() + 3 == task.n_heavy:
            smiles_list = json.loads(task.smiles_list)
            new_smiles_list = []
            for smiles in smiles_list:
                rdk_mol = Chem.MolFromSmiles(smiles)
                for atom in rdk_mol.GetAtoms():
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

def create_FluoroAlkylAcetate(n_heavy):
    print('Create FluoroAlkylAcetate SMILES\n')
    category = 'FluoroAlkylAcetate'
    n = n_heavy
    id = 0
    while id < n:
        tasks_all = session.query(Task)
        tasks = tasks_all.filter(Task.category == category).filter(Task.n_heavy == id + 1)
        if tasks.count() == 1:
            id += 1
            continue
        elif tasks.count() == 0:
            tasks = tasks_all.filter(Task.category == 'AlkylAcetate').filter(Task.n_heavy == id + 1)
            if tasks.count() == 0:
                id += 1
                continue
            task = tasks.first()
            smiles_list = json.loads(task.smiles_list)
            new_smiles_list = []
            for smiles in smiles_list:
                rdk_mol = Chem.MolFromSmiles(smiles)
                new_mol = Chem.RWMol(Chem.AddHs(rdk_mol))
                for atom in new_mol.GetAtoms():
                    if atom.GetAtomicNum() == 1:
                        new_mol.ReplaceAtom(atom.GetIdx(), Chem.Atom(9))
                new_smiles = get_canonical_smiles(Chem.MolToSmiles(new_mol))
                if new_smiles not in new_smiles_list:
                    new_smiles_list.append(new_smiles)
            new_info = Task()
            new_info.smiles_list = json.dumps(new_smiles_list)
            new_info.n_heavy = get_heavy_atom_numbers(new_smiles_list[0])
            new_info.category = category
            new_info.molecule_number = len(new_smiles_list)
            if new_info.smiles_list != '[]':
                session.add(new_info)
            session.commit()
            id += 1
        else:
            raise Exception('Bad infomation error\n')
        sys.stdout.write('\r present %i / %i' % (id, n))

def create_AlkylSulfonate(n_heavy):
    print('Create AlkylSulfonate SMILES\n')
    category = 'AlkylSulfonate'
    n = n_heavy
    id = 0
    while id < n:
        tasks_all = session.query(Task)
        tasks = tasks_all.filter(Task.category==category).order_by(Task.n_heavy.desc())
        task = tasks.first()
        if tasks.count() == 0:
            new_info = Task()
            new_info.smiles_list = json.dumps([get_canonical_smiles('CS(=O)(=O)[O-]')])
            new_info.n_heavy = 5
            new_info.category = category
            new_info.molecule_number = 1
            session.add(new_info)
            session.commit()
            id = 5
        elif task.n_heavy >= n:
            break
        elif tasks.count() + 4 == task.n_heavy:
            smiles_list = json.loads(task.smiles_list)
            new_smiles_list = []
            for smiles in smiles_list:
                rdk_mol = Chem.MolFromSmiles(smiles)
                for atom in rdk_mol.GetAtoms():
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

def create_FluoroAlkylSulfonate(n_heavy):
    print('Create FluoroAlkylSulfonate SMILES\n')
    category = 'FluoroAlkylSulfonate'
    n = n_heavy
    id = 0
    while id < n:
        tasks_all = session.query(Task)
        tasks = tasks_all.filter(Task.category == category).filter(Task.n_heavy == id + 1)
        if tasks.count() == 1:
            id += 1
            continue
        elif tasks.count() == 0:
            tasks = tasks_all.filter(Task.category == 'AlkylSulfonate').filter(Task.n_heavy == id + 1)
            if tasks.count() == 0:
                id += 1
                continue
            task = tasks.first()
            smiles_list = json.loads(task.smiles_list)
            new_smiles_list = []
            for smiles in smiles_list:
                rdk_mol = Chem.MolFromSmiles(smiles)
                new_mol = Chem.RWMol(Chem.AddHs(rdk_mol))
                for atom in new_mol.GetAtoms():
                    if atom.GetAtomicNum() == 1:
                        new_mol.ReplaceAtom(atom.GetIdx(), Chem.Atom(9))
                new_smiles = get_canonical_smiles(Chem.MolToSmiles(new_mol))
                if new_smiles not in new_smiles_list:
                    new_smiles_list.append(new_smiles)
            new_info = Task()
            new_info.smiles_list = json.dumps(new_smiles_list)
            new_info.n_heavy = get_heavy_atom_numbers(new_smiles_list[0])
            new_info.category = category
            new_info.molecule_number = len(new_smiles_list)
            if new_info.smiles_list != '[]':
                session.add(new_info)
            session.commit()
            id += 1
        else:
            raise Exception('Bad infomation error\n')
        sys.stdout.write('\r present %i / %i' % (id, n))

def create_AlkylSulfonylImide(n_heavy):
    print('Create AlkylSulfonylImide SMILES\n')
    category = 'AlkylSulfonylImide'
    n = n_heavy
    id = 0
    while id < n:
        tasks_all = session.query(Task)
        tasks = tasks_all.filter(Task.category==category).order_by(Task.n_heavy.desc())
        task = tasks.first()
        if tasks.count() == 0:
            new_info = Task()
            new_info.smiles_list = json.dumps([get_canonical_smiles('C(S(=O)(=O)[N-]S(=O)(=O)C)')])
            new_info.n_heavy = 9
            new_info.category = category
            new_info.molecule_number = 1
            session.add(new_info)
            session.commit()
            id = 9
        elif task.n_heavy >= n:
            break
        elif tasks.count() + 8 == task.n_heavy:
            smiles_list = json.loads(task.smiles_list)
            new_smiles_list = []
            for smiles in smiles_list:
                rdk_mol = Chem.MolFromSmiles(smiles)
                for atom in rdk_mol.GetAtoms():
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

def create_FluoroAlkylSulfonylImide(n_heavy):
    print('Create FluoroAlkylSulfonylImide SMILES\n')
    category = 'FluoroAlkylSulfonylImide'
    n = n_heavy
    id = 0
    while id < n:
        tasks_all = session.query(Task)
        tasks = tasks_all.filter(Task.category == category).filter(Task.n_heavy == id + 1)
        if tasks.count() == 1:
            id += 1
            continue
        elif tasks.count() == 0:
            tasks = tasks_all.filter(Task.category == 'AlkylSulfonylImide').filter(Task.n_heavy == id + 1)
            if tasks.count() == 0:
                id += 1
                continue
            task = tasks.first()
            smiles_list = json.loads(task.smiles_list)
            new_smiles_list = []
            for smiles in smiles_list:
                rdk_mol = Chem.MolFromSmiles(smiles)
                new_mol = Chem.RWMol(Chem.AddHs(rdk_mol))
                for atom in new_mol.GetAtoms():
                    if atom.GetAtomicNum()==1:
                        new_mol.ReplaceAtom(atom.GetIdx(), Chem.Atom(9))
                new_smiles = get_canonical_smiles(Chem.MolToSmiles(new_mol))
                if new_smiles not in new_smiles_list:
                    new_smiles_list.append(new_smiles)
            new_info = Task()
            new_info.smiles_list = json.dumps(new_smiles_list)
            new_info.n_heavy = get_heavy_atom_numbers(new_smiles_list[0])
            new_info.category = category
            new_info.molecule_number = len(new_smiles_list)
            if new_info.smiles_list != '[]':
                session.add(new_info)
            session.commit()
            id += 1
        else:
            raise Exception('Bad infomation error\n')
        sys.stdout.write('\r present %i / %i' % (id, n))