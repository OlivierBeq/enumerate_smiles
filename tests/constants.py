# -*- coding: utf-8 -*-


"""Constant for enumerate_siles unittesting."""

import os

from rdkit import Chem

from enumerate_smiles import MolSupplier

SAMPLES_FOLDER = os.path.realpath(os.path.join(os.path.dirname(__file__), 'samples'))

MOL_SUPPLIER = MolSupplier(os.path.join(SAMPLES_FOLDER, 'CID148124.mol'))
MOL2_SUPPLIER = MolSupplier(os.path.join(SAMPLES_FOLDER, 'CID148124.mol2'))
SMILES_SUPPLIER = MolSupplier(os.path.join(SAMPLES_FOLDER, 'testmol.smi'), titleLine=False)
MAE_SUPPLIER = MolSupplier(os.path.join(SAMPLES_FOLDER, 'CID18730.mae'))
SDF_SUPPLIER = MolSupplier(os.path.join(SAMPLES_FOLDER, 'drugs.sdf'))
