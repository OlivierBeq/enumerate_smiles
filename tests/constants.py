# -*- coding: utf-8 -*-


"""Constant for chemopy unittesting."""

import os

from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom


SAMPLES_FOLDER = os.path.realpath(os.path.join(os.path.dirname(__file__), 'samples'))

MOL_MOLECULE = Chem.MolFromMolFile(os.path.join(SAMPLES_FOLDER, 'CID148124.mol'))
MOL2_MOLECULE = Chem.MolFromMol2File(os.path.join(SAMPLES_FOLDER, 'CID148124.mol2'))

with open(os.path.join(SAMPLES_FOLDER, 'testmol.smi')) as fh:
    SMILES_MOLECULES = [Chem.MolFromSmiles(line.strip()) for line in fh]

with open(os.path.join(SAMPLES_FOLDER, 'cas.txt')) as fh:
    CAS_NUMBERS = [line.strip() for line in fh.readlines()]

with open(os.path.join(SAMPLES_FOLDER, 'drugbank.txt')) as fh:
    DRUGBANK_NUMBERS = [line.strip() for line in fh.readlines()]

with open(os.path.join(SAMPLES_FOLDER, 'kegg.txt')) as fh:
    KEGG_NUMBERS = [line.strip() for line in fh.readlines()]

with open(os.path.join(SAMPLES_FOLDER, 'ncbi.txt')) as fh:
    NCBI_NUMBERS = [line.strip() for line in fh.readlines()]

with open(os.path.join(SAMPLES_FOLDER, 'ebi.txt')) as fh:
    EBI_NUMBERS = [line.strip() for line in fh.readlines()]

MOLECULES_2D = [MOL_MOLECULE, MOL2_MOLECULE] + SMILES_MOLECULES

MOLECULES_3D = [Chem.AddHs(mol) for mol in MOLECULES_2D]
_embedded_3d = map(AllChem.EmbedMolecule, MOLECULES_3D)
MOLECULES_3D = [mol for mol, success in zip(MOLECULES_3D, _embedded_3d) if success == 0]
