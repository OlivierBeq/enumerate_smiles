# -*- coding: utf-8 -*-

"""Tests for details about molecular descriptors/fingerprints."""

import unittest

from rdkit import Chem
from rdkit.Chem import AllChem
from chemopy import ChemoPy

from .constants import MOL_MOLECULE


class TestDetailsOfDescriptorsAndFingerprints(unittest.TestCase):
    """Tests molecular descriptors."""

    def setUp(self) -> None:
        """Create the molecular descriptor calculator."""
        self.cmp = ChemoPy(ignore_3D=False)
        mol = Chem.Mol(MOL_MOLECULE)
        _ = AllChem.EmbedMolecule(mol)
        self.desc_names = self.cmp.calculate([mol], show_banner=False).columns.tolist()

    def test_each_descriptor_details(self):
        """Test the details for each descriptor available."""
        for desc_name in self.desc_names:
            details = self.cmp.get_details(desc_name)
            self.assertIsNotNone(details)
            self.assertFalse(details.empty)
            self.assertEqual(details.shape[0], 1)
            self.assertEqual(details.Name.item(), desc_name)

    def test_all_descriptor_details(self):
        """Test the details of all descriptors/fingerprints available."""
        details = self.cmp.get_details()
        self.assertIsNotNone(details)
        self.assertFalse(details.empty)
        # Includes the 12 additional fingerprints
        self.assertEqual(details.shape[0], 1196)
