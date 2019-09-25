import numpy as np
import requests, pybel, json, sys, math
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker, relationship
from sqlalchemy import create_engine
from sqlalchemy import Column, Integer, Float, Text, Boolean, String, ForeignKey, UniqueConstraint
sys.path.append('..')
from config import *
sys.path.append(Config.MS_TOOLS_DIR)
from mstools.smiles.smiles import *

Base = declarative_base()
metadata = Base.metadata

db_file = 'sqlite:///../database/smiles.sqlite?check_same_thread=False'

engine = create_engine(db_file, echo=False)
Session = sessionmaker(engine)
session = Session()


class Task(Base):
    __tablename__ = 'task'
    id = Column(Integer, primary_key=True)
    smiles_list = Column(Text, unique=True)
    n_heavy = Column(Integer)
    category = Column(Text)
    molecule_number = Column(Integer)
    remark = Column(Text)
    charge = Column(Integer)
    molecules = relationship('Molecule', lazy='dynamic')

    def __repr__(self):
        return '<Task: %i %s %i>' % (self.id, self.category, self.n_heavy)

    def generate_Molecule_table(self, smiles_unique_check=False):
        print('generate_Molecule_table task_id = %i' % self.id)
        smiles_list = json.loads(self.smiles_list)
        if self.molecules.count() == len(smiles_list):
            return
        for smiles in smiles_list:
            if smiles_unique_check:
                mols = session.query(Molecule).filter(Molecule.smiles == smiles)
                if mols.count() == 1:
                    print('%i, %s' % (mols.first().id, mols.first().smiles))
                    continue
            mol = Molecule(task=self)
            mol.smiles = smiles
            if self.remark == 'unstable':
                mol.stability = False
            else:
                mol.stability = True
            smiles_stereo_isomer = get_stereo_isomer(smiles)
            if len(smiles_stereo_isomer) == 1:
                mol.has_StereoIsomer = False
            else:
                mol.has_StereoIsomer = True
            session.add(mol)
            mol.generate_Stereoisomer_table()
            session.commit()

    def delete(self):
        print('delete task_id = %i' % self.id)
        self.delete_molecule()
        session.delete(self)
        session.commit()

    def delete_molecule(self):
        print('delete molecules of task_id = %i' % self.id)
        for mol in self.molecules:
            session.delete(mol)
        session.commit()

    def get_mol_number_list(self, similarity=False):
        stereos_all = session.query(StereoIsomer)
        n = 10000
        N = math.ceil(stereos_all.count() / n)
        if similarity:
            count = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            for j in range(N):
                stereos = stereos_all.filter(StereoIsomer.id >= j * n).filter(StereoIsomer.id < (j + 1) * n)
                for i, stereo in enumerate(stereos):
                    if i % 100 == 0:
                        sys.stdout.write('\r %i / %i' % (i + j * n, N * n))
                    if stereo.molecule.task == self:
                        if stereo.similarity_010:
                            count[0] += 1
                        if stereo.similarity_020:
                            count[1] += 1
                        if stereo.similarity_030:
                            count[2] += 1
                        if stereo.similarity_050:
                            count[3] += 1
                        if stereo.similarity_080:
                            count[4] += 1
                        if stereo.similarity_100:
                            count[5] += 1
                        if stereo.similarity_150:
                            count[6] += 1
                        if stereo.similarity_200:
                            count[7] += 1
                        if stereo.training:
                            count[8] += 1
                        count[9] += 1
        else:
            count = [0, 0]
            for j in range(N):
                stereos = stereos_all.filter(StereoIsomer.id >= j * n).filter(StereoIsomer.id < (j + 1) * n)
                if stereos[0].molecule.task != self and stereos[-1].molecule.task != self and j != 0:
                    sys.stdout.write('\r %i / %i' % (j * n, N * n))
                    if count != [0, 0]:
                        return count
                    continue
                for i, stereo in enumerate(stereos):
                    if i % 100 == 0:
                        sys.stdout.write('\r %i / %i' % (i + j * n, N * n))
                    if stereo.molecule.task == self:
                        if stereo.training:
                            count[0] += 1
                        count[1] += 1
        return count


class Molecule(Base):
    __tablename__ = 'molecule'
    id = Column(Integer, primary_key=True)
    smiles = Column(Text, unique=True)
    task_id = Column(Integer, ForeignKey(Task.id))
    stability = Column(Boolean)
    has_StereoIsomer = Column(Boolean)
    remark = Column(Text)

    task = relationship('Task', foreign_keys='Molecule.task_id')
    stereoisomers = relationship('StereoIsomer', lazy='dynamic')

    def generate_Stereoisomer_table(self, topolfp=True, badfp=False):
        smiles_stereo_isomers = get_stereo_isomer(self.smiles)
        for smiles in smiles_stereo_isomers:
            isomer = StereoIsomer(molecule=self)
            isomer.smiles = smiles
            if Config.TOPOLOGICAL:
                isomer.generate_rdk_fp()
            if Config.MORGAN:
                isomer.generate_morgan_fp()
            if Config.PAIR:
                isomer.generate_pair_fp()
            if Config.TORSION:
                isomer.generate_torsion_fp()
            if Config.PATH:
                isomer.generate_path_fp()
            session.add(isomer)
        session.commit()

    def delete(self):
        for isomer in self.stereoisomers:
            session.delete(isomer)
        session.delete(self)
        session.commit()


class StereoIsomer(Base):
    __tablename__ = 'stereoisomer'
    id = Column(Integer, primary_key=True)
    smiles = Column(Text, unique=True)
    molecule_id = Column(Integer, ForeignKey(Molecule.id))
    simple_fp = Column(Text)
    rdk_fp = Column(Text)
    pair_fp = Column(Text)
    torsion_fp = Column(Text)
    path_fp = Column(Text)
    morgan_fp = Column(Text)
    final_fp = Column(Text)
    training = Column(Boolean, default=False)
    similarity_010 = Column(Boolean, default=False)
    similarity_020 = Column(Boolean, default=False)
    similarity_030 = Column(Boolean, default=False)
    similarity_050 = Column(Boolean, default=False)
    similarity_080 = Column(Boolean, default=False)
    similarity_100 = Column(Boolean, default=False)
    similarity_150 = Column(Boolean, default=False)
    similarity_200 = Column(Boolean, default=False)

    molecule = relationship('Molecule', foreign_keys='StereoIsomer.molecule_id')

    def generate_simple_fp(self):
        self.simple_fp = json.dumps(get_simple_fp(self.smiles))
        session.commit()

    def generate_rdk_fp(self, minPath=Config.TOPOLOGICAL_minPath, maxPath=Config.TOPOLOGICAL_maxPath):
        self.rdk_fp = json.dumps(get_fingerprint(self.smiles, type='rdk', count=True, minPath=minPath, maxPath=maxPath))
        session.commit()

    def generate_path_fp(self, maxPath=Config.PATH_maxPath + 1):
        fp_list = []
        for npath in range(1, maxPath):
            rdk_fp = get_fingerprint(self.smiles, type='rdk', count=True, minPath=npath, maxPath=npath)
            fp_list.append(sum(rdk_fp.values()))
        self.path_fp = json.dumps(fp_list)
        session.commit()

    def generate_pair_fp(self):
        self.pair_fp = json.dumps(get_fingerprint(self.smiles, type='pair', count=True))
        session.commit()

    def generate_torsion_fp(self):
        self.torsion_fp = json.dumps(get_fingerprint(self.smiles, type='torsion', count=True))
        session.commit()

    def generate_morgan_fp(self, radius=Config.MORGAN_radius):
        info = get_fingerprint(self.smiles, type='morgan', radius=radius, count=True)
        for name in info.keys():
            info[name] = len(info[name])
        self.morgan_fp = json.dumps(info)


def create_table():
    if not os.path.exists(os.path.join('..', 'database')):
        os.mkdir(os.path.join('..', 'database'))
    metadata.create_all(engine)

