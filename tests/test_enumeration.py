# -*- coding: utf-8 -*-

"""Tests for molecular descriptors."""

import unittest

from enumerate_smiles import SmilesEnumerator

from tests.constants import *


class TestEnumeration(unittest.TestCase):
    """Tests enumeration."""

    def test_mol_enumeration(self):
        """Test enumeration from the MOL supplier."""
        se = SmilesEnumerator()
        smiles = se.enumerate(MOL_SUPPLIER)
        self.assertEqual(len(smiles), 1000)
        MOL_SUPPLIER.reset()

    def test_mol2_enumeration(self):
        """Test enumeration from the MOL2 supplier."""
        se = SmilesEnumerator()
        smiles = se.enumerate(MOL2_SUPPLIER)
        self.assertEqual(len(smiles), 1000)
        MOL2_SUPPLIER.reset()

    def test_mae_enumeration(self):
        """Test enumeration from the MAE supplier."""
        se = SmilesEnumerator()
        smiles = se.enumerate(MAE_SUPPLIER)
        self.assertEqual(len(smiles), 1000)
        MAE_SUPPLIER.reset()

    def test_smiles_enumeration(self):
        """Test enumeration from the SMILES supplier."""
        se = SmilesEnumerator(max_enum_smiles=60, max_out_smiles=10)
        smiles = se.enumerate(SMILES_SUPPLIER, njobs=-1)
        self.assertEqual(len(smiles), 400)
        SMILES_SUPPLIER.reset()

    def test_sdf_enumeration(self):
        """Test enumeration from the SDF supplier."""
        se = SmilesEnumerator(max_enum_smiles=50, max_out_smiles=10)
        smiles = se.enumerate(SDF_SUPPLIER, njobs=-1)
        self.assertEqual(len(smiles), 50)
        SDF_SUPPLIER.reset()

    def test_mol_withH_enumeration(self):
        """Test enumeration from the MOL supplier."""
        se = SmilesEnumerator(max_enum_smiles=60, max_out_smiles=10)
        smiles = se.enumerate((Chem.AddHs(mol) for mol in MOL_SUPPLIER))
        self.assertEqual(len(smiles), 10)
        MOL_SUPPLIER.reset()

    def test_mol2_withH_enumeration(self):
        """Test enumeration from the MOL2 supplier."""
        se = SmilesEnumerator(max_enum_smiles=100, max_out_smiles=10)
        smiles = se.enumerate((Chem.AddHs(mol) for mol in MOL2_SUPPLIER))
        self.assertEqual(len(smiles), 10)
        MOL2_SUPPLIER.reset()

    def test_mae_withH_enumeration(self):
        """Test enumeration from the MAE supplier."""
        se = SmilesEnumerator(max_enum_smiles=60, max_out_smiles=10)
        smiles = se.enumerate((Chem.AddHs(mol) for mol in MAE_SUPPLIER))
        self.assertEqual(len(smiles), 10)
        MAE_SUPPLIER.reset()

    def test_smiles_withH_enumeration(self):
        """Test enumeration from the SMILES supplier."""
        se = SmilesEnumerator(max_enum_smiles=60, max_out_smiles=10)
        smiles = se.enumerate((Chem.AddHs(mol) for mol in SMILES_SUPPLIER), njobs=-1)
        self.assertEqual(len(smiles), 400)
        SMILES_SUPPLIER.reset()

    def test_sdf_withH_enumeration(self):
        """Test enumeration from the SDF supplier."""
        se = SmilesEnumerator(max_enum_smiles=50, max_out_smiles=10)
        smiles = se.enumerate((Chem.AddHs(mol) for mol in SDF_SUPPLIER), njobs=-1)
        self.assertEqual(len(smiles), 50)
        SDF_SUPPLIER.reset()
