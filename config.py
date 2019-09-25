import os

CWD = os.path.dirname(os.path.abspath(__file__))


class Config:
    classes = ['alkane', 'alcohol', 'ketone', 'ether', 'ammounia']

    TOPOLOGICAL = True
    TOPOLOGICAL_minPath = 1
    TOPOLOGICAL_maxPath = 7
    MORGAN = False
    MORGAN_radius = 1
    PAIR = False
    TORSION = False
    PATH = False
    PATH_maxPath = 10

    if os.path.exists(os.path.join(CWD, '..', 'AIMS_Tools')):
        MS_TOOLS_DIR = os.path.join(CWD, '..', 'AIMS_Tools')
    elif os.path.exists('/share/md1400/xiangyan/Github/AIMS_Tools'):
        MS_TOOLS_DIR = '/share/md1400/xiangyan/Github/AIMS_Tools'
    else:
        raise Exception('AIMS_Tools not found')
