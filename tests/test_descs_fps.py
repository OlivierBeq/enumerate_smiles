# -*- coding: utf-8 -*-

"""Tests for molecular descriptors."""

import unittest

from chemopy import ChemoPy, Fingerprint, Fingerprint3D

from tests.constants import *


class TestDescriptorsAndFingerprints(unittest.TestCase):
    """Tests molecular descriptors."""

    def setUp(self) -> None:
        """Create the molecular descriptor calculator."""
        self.cmp = ChemoPy(ignore_3D=True)
        self.cmp3d = ChemoPy(ignore_3D=False)
        self.cmpfp = ChemoPy(ignore_3D=True, include_fps=True)
        self.cmp3dfp = ChemoPy(ignore_3D=False, include_fps=True)

    def test_2D_descriptors(self):
        """Test the dimensions of the output descriptor dataframe."""
        values = self.cmp.calculate(MOLECULES_2D, show_banner=False)
        self.assertEqual(values.shape, (len(MOLECULES_2D), 632))
        self.assertEqual(len(values.columns.unique().tolist()), 632)

    def test_2D_descriptors_multithread(self):
        """Test the dimensions of the output descriptor dataframes calculated by different processes."""
        values = self.cmp.calculate(MOLECULES_2D, show_banner=False, njobs=-1, chunksize=1)
        self.assertEqual(values.shape, (len(MOLECULES_2D), 632))
        self.assertEqual(len(values.columns.unique().tolist()), 632)

    def test_3D_descriptors(self):
        """Test the dimensions of the output descriptor dataframe."""
        values = self.cmp3d.calculate(MOLECULES_3D[:2], show_banner=False)
        self.assertEqual(values.shape, (len(MOLECULES_3D[:2]), 1184))
        self.assertEqual(len(values.columns.unique().tolist()), 1184)

    def test_3D_descriptors_multithread(self):
        """Test the dimensions of the output descriptor dataframes calculated by different processes."""
        values = self.cmp3d.calculate(MOLECULES_3D[:2], show_banner=False, njobs=-1, chunksize=1)
        self.assertEqual(values.shape, (len(MOLECULES_3D[:2]), 1184))
        self.assertEqual(len(values.columns.unique().tolist()), 1184)

    def test_2D_fingerprints_2D_descriptors(self):
        """Test the dimensions of the output descriptor dataframe."""
        values = self.cmpfp.calculate(MOLECULES_2D, show_banner=False)
        self.assertEqual(values.shape, (len(MOLECULES_2D), 14551))
        self.assertEqual(len(values.columns.unique().tolist()), 14551)

    def test_2D_fingerprints_2D_descriptors_multithread(self):
        """Test the dimensions of the output descriptor dataframes calculated by different processes."""
        values = self.cmpfp.calculate(MOLECULES_2D, show_banner=False, njobs=-1, chunksize=1)
        self.assertEqual(values.shape, (len(MOLECULES_2D), 14551))
        self.assertEqual(len(values.columns.unique().tolist()), 14551)

    def test_2D_fingerprints(self):
        for mol in MOLECULES_2D:
            values = Fingerprint.get_all_fps(mol)
            self.assertEqual(len(values), 13919)

    def test_3D_fingerprints(self):
        for mol in MOLECULES_3D:
            values = Fingerprint3D.get_all_fps(mol)
            self.assertEqual(len(values), 2048)

    def test_3D_fingerprints_2D_descriptors(self):
        """Test the dimensions of the output descriptor dataframe."""
        values = self.cmp3dfp.calculate(MOLECULES_3D[:2], show_banner=False)
        self.assertEqual(values.shape, (len(MOLECULES_3D[:2]), 17151))
        self.assertEqual(len(values.columns.unique().tolist()), 17151)

    def test_3D_fingerprints_2D_descriptors_multithread(self):
        """Test the dimensions of the output descriptor dataframes calculated by different processes."""
        values = self.cmp3dfp.calculate(MOLECULES_3D[:2], show_banner=False, njobs=-1, chunksize=1)
        self.assertEqual(values.shape, (len(MOLECULES_3D[:2]), 17151))
        self.assertEqual(len(values.columns.unique().tolist()), 17151)
