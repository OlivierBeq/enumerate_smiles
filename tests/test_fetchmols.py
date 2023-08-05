# -*- coding: utf-8 -*-

"""Tests for the fetching of molecules from online resources."""

import unittest

from chemopy import MolFrom

from tests.constants import *


class TestDescriptorsAndFingerprints(unittest.TestCase):
    """Tests for online molecule retrieval."""

    def test_cas(self):
        for casnum in CAS_NUMBERS:
            mol = MolFrom.CAS(casnum)
            self.assertIsNotNone(mol)
            self.assertNotEqual(Chem.MolToSmiles(mol), '')

    def test_ebi(self):
        for ebinum in EBI_NUMBERS:
            mol = MolFrom.EBI(ebinum)
            self.assertIsNotNone(mol)
            self.assertNotEqual(Chem.MolToSmiles(mol), '')

    def test_ncbi(self):
        for ncbinum in NCBI_NUMBERS:
            mol = MolFrom.NCBI(ncbinum)
            self.assertIsNotNone(mol)
            self.assertNotEqual(Chem.MolToSmiles(mol), '')

    def test_drugbank(self):
        for drugbankid in DRUGBANK_NUMBERS:
            mol = MolFrom.DrugBank(drugbankid)
            self.assertIsNotNone(mol)
            self.assertNotEqual(Chem.MolToSmiles(mol), '')

    def test_kegg(self):
        for keggid in KEGG_NUMBERS:
            mol = MolFrom.KEGG(keggid)
            self.assertIsNotNone(mol)
            self.assertNotEqual(Chem.MolToSmiles(mol), '')