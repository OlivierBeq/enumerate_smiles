# -*- coding: utf-8 -*-

import itertools
from typing import Iterator, List
from multiprocessing import cpu_count

import numpy as np
from rdkit import Chem

from rdkit.Chem import AllChem, MolStandardize
from rdkit.rdBase import BlockLogs
from rdkit.Chem.EnumerateHeterocycles import EnumerateHeterocycles
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
from tqdm.auto import tqdm
from bounded_pool_executor import BoundedProcessPoolExecutor

from .IO import MolSupplier

class SmilesEnumerator:

    def __init__(self,
                 enum_heterocycles: bool = True,
                 enum_tautomers: bool = True,
                 enum_resonanceforms: bool = True,
                 enum_stereoisomers: bool = True,
                 enum_smiles: bool = True,
                 max_num_smiles: int = 1000,
                 smiles_type: str = 'both',
                 random_state: int = 1234
                 ) -> None:
        """Create a SMILES enumerator.

        :param enum_heterocycles: enumerate all derived heterocycles
        :param enum_tautomers: enumerate all tautomeric forms
        :param enum_resonanceforms: enumerate resonance forms
        :param enum_stereoisomers: enumerate stereoisomers after having dropped the input stereochemistry
        :param enum_smiles: enumerate SMILES
        :param max_num_smiles: maximum number of SMILES to sample per molecule
        :param smiles_type: type of output SMILES: {'kekule', 'aromatic', 'both'}
        :param random_state: random seed for the selection of enumerated SMILES to be returned
        """
        if smiles_type not in ('both', 'kekule', 'aromatic'):
            raise ValueError("smiles_type can only be one of ('both', 'kekule', 'aromatic')")
        self.enum_heterocycles = enum_heterocycles
        self.enum_tautomers = enum_tautomers
        self.enum_resonanceforms = enum_resonanceforms
        self.enum_stereoisomers = enum_stereoisomers
        self.enum_smiles = enum_smiles
        self.max_num_smiles = max_num_smiles
        self.smiles_type = smiles_type
        # Generators
        self.tauto_enumerator = Chem.MolStandardize.rdMolStandardize.TautomerEnumerator()
        self.random = np.random.default_rng(random_state)
        self.block = BlockLogs()

    def __del__(self):
        del self.block

    def enumerate(self, mols: List[Chem.Mol], njobs: int = -1) -> List[str]:
        """Enumerate SMILES.

        :param mols: molecules to enumerate the SMILES of.
        :param njobs: number of parallel jobs
        :return: the list of enumerated SMILES
        """
        # Parallelize should need be
        if njobs != 1:
            if njobs < 0:
                njobs += cpu_count()
            with BoundedProcessPoolExecutor(max_workers=njobs) as worker:
                futures = [worker.submit(self._enumerator(mol))
                           for mol in mols]
            # Collect results
            return list(itertools.chain.from_iterable(future.result() for future in futures))
        else:
            # Single process
            result = list(itertools.chain.from_iterable(self._enumerator(mol) for mol in mols))
        return result

    def _sanitize(self, mol: Chem.Mol):
        """Sanitize and return the input molecule or None upon failure."""
        if Chem.SanitizeMol(mol, catchErrors=True) != 0:
            return None
        return mol

    def _copy_props(self, from_mol: Chem.Mol, to_mol: Chem.Mol):
        """Copy all properties of a molecule to another."""
        props = from_mol.GetPropsAsDict(includePrivate=True)
        for key, value in props.items():
            if isinstance(value, int):
                to_mol.SetIntProp(key, value)
            elif isinstance(value, float):
                to_mol.SetDoubleProp(key, value)
            elif isinstance(value, bool):
                to_mol.SetBoolProp(key, value)
            else:
                to_mol.SetProp(key, str(value))
        return to_mol

    def _enumerator(self, mol: Chem.Mol) -> List[str]:
        """Obtain an iterator of SMILES for one single molecule.

        :param mol: Molecule to iterate the SMILES of.
        :return: A generator of SMILES
        """
        mols = [mol]
        if self.enum_heterocycles:
            mols = list(itertools.chain(mols,
                                        [self._copy_props(mol, self._sanitize(het))
                                         for mol in mols
                                         for het in EnumerateHeterocycles(mol)]))
            mols = self.random.choice(mols, self.max_num_smiles)
        if self.enum_tautomers:
            mols = list(itertools.chain(mols,
                                        [self._copy_props(mol, self._sanitize(tauto))
                                         for mol in mols
                                         for tauto in self.tauto_enumerator.Enumerate(mol)
                                         if mol is not None]))
            mols = self.random.choice(mols, self.max_num_smiles)
        if self.enum_resonanceforms:
            mols = list(itertools.chain(mols,
                                        [self._copy_props(mol, self._sanitize(res))
                                         for mol in mols
                                         for res in Chem.ResonanceMolSupplier(mol)
                                         if mol is not None]))
            mols = self.random.choice(mols, self.max_num_smiles)
        if self.enum_stereoisomers:
            mols = list(itertools.chain(mols,
                                        [self._copy_props(mol, self._sanitize(stereoi))
                                         for mol in mols
                                         for stereoi in EnumerateStereoisomers(mol)
                                         if mol is not None]))
            mols = self.random.choice(mols, self.max_num_smiles)
        if self.enum_smiles:
            if self.smiles_type == 'both':
                mols1 = [(smiles, mol.GetPropsAsDict(includePrivate=True).get('_Name', ''))
                        for mol in mols
                        for smiles in list(set(Chem.MolToSmiles(mol, doRandom=True, canonical=False, kekuleSmiles=True)
                                               for _ in range(self.max_num_smiles * 10)
                                               ))
                        if mol is not None
                        ]
                mols1 = self.random.choice(mols1, int(self.max_num_smiles / 2))
                mols2 = [(smiles, mol.GetPropsAsDict(includePrivate=True).get('_Name', ''))
                         for mol in mols
                         for smiles in list(set(Chem.MolToSmiles(mol, doRandom=True, canonical=False, kekuleSmiles=False)
                                                for _ in range(self.max_num_smiles * 10)
                                                ))
                         if mol is not None
                        ]
                mols2 = self.random.choice(mols2, int(self.max_num_smiles / 2))
                mols = itertools.chain(mols1, mols2)
            elif self.smiles_type == 'kekule':
                mols = [(smiles, mol.GetPropsAsDict(includePrivate=True).get('_Name', ''))
                        for mol in mols
                        for smiles in list(set(Chem.MolToSmiles(mol, doRandom=True, canonical=False, kekuleSmiles=True)
                                               for _ in range(self.max_num_smiles * 10)
                                               ))
                        if mol is not None
                        ]
                mols = self.random.choice(mols, self.max_num_smiles)
            else: # aromatic
                mols = [(smiles, mol.GetPropsAsDict(includePrivate=True).get('_Name', ''))
                        for mol in mols
                        for smiles in list(set(Chem.MolToSmiles(mol, doRandom=True, canonical=False, kekuleSmiles=False)
                                               for _ in range(self.max_num_smiles * 10)
                                               ))
                        if mol is not None
                        ]
                mols = self.random.choice(mols, self.max_num_smiles)
        else:
            mols = [(Chem.MolToSmiles(mol), mol.GetPropsAsDict(includePrivate=True).get('_Name', ''))
                    for mol in mols
                    if mol is not None
                    ]
        return mols

    def from_file_to_file(self,
                          molsupplier: MolSupplier,
                          out_file: str,
                          max_in_mols: int = None,
                          njobs: int = -1,
                          progress: bool = True
                          ) -> None:
        """Enumerate all derivatives of SMILES contained in the input file and write them in the output file.

        :param in_file: MolSupplier of the input file containing molecules to enumerate SMILES from.
        :param out_file: output SMILES file to write results to.
        :param max_in_mols: maximum number of molecules of the input file to be enumerated
        :param njobs: number of parallel processes
        :param progress: should progress be displayed
        """
        # Determine the number of molecules
        len_supplier = len(molsupplier)
        if max_in_mols is not None:
            molsupplier = itertools.islice(molsupplier, max_in_mols)
            len_supplier = max_in_mols
        # Get molecules
        mols = list(molsupplier)
        # Enumerate SMILES
        if njobs != 1:
            if njobs < 0:
                njobs += cpu_count() + 1
            with BoundedProcessPoolExecutor(max_workers=njobs) as worker:
                futures = [worker.submit(self._enumerator(mol))
                           for mol in mols]
            # Collect results
            smiles = [future.result() for future in futures]
            # smiles = itertools.chain.from_iterable(future.result() for future in futures)
        else:
            # Single process
            smiles = itertools.chain.from_iterable(self._enumerator(mol) for mol in mols)
        # Write results to output
        with open(out_file, 'w') as oh:
            pbar = tqdm(smiles, total=len_supplier * self.max_num_smiles) if progress else smiles
            for smiles, mol_id in pbar:
                oh.write(f'{smiles}\t{mol_id}\n')
