# -*- coding: utf-8 -*-

import itertools
import multiprocessing
from functools import partial
from typing import List, Callable, Union

import numpy as np
from bounded_pool_executor import BoundedProcessPoolExecutor
from more_itertools import chunked
from rdkit import Chem
from rdkit.Chem import MolStandardize
from rdkit.Chem.EnumerateHeterocycles import EnumerateHeterocycles
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from rdkit.rdBase import BlockLogs
from tqdm.auto import tqdm

from .IO import MolSupplier


class SmilesEnumerator:

    def __init__(self,
                 enum_heterocycles: bool = True,
                 enum_tautomers: bool = True,
                 enum_resonanceforms: bool = True,
                 enum_stereoisomers: bool = True,
                 enum_smiles: bool = True,
                 max_enum_heterocycles: int = 100,
                 max_enum_tautomers: int = 100,
                 max_enum_resonanceforms: int = 100,
                 max_enum_stereoisomers: int = 100,
                 max_enum_smiles: int = 100,
                 random_choice_heterocycles: bool = True,
                 random_choice_tautomers: bool = True,
                 random_choice_resonanceforms: bool = True,
                 random_choice_stereoisomers: bool = True,
                 max_out_heterocycles: int = 70,
                 max_out_tautomers: int = 125,
                 max_out_resonanceforms: int = 250,
                 max_out_stereoisomers: int = 500,
                 max_out_smiles: int = 1000,
                 smiles_type: str = 'both',
                 random_state: int = 1234
                 ) -> None:
        """Create a SMILES enumerator.

        The process follows the following workflow:
            1) a - enumerate heterocycles from input molecules
               b - combine input molecules and outputs of step 1a
               c - obtain a random subsample
            2) a - enumerate tautomers from outputs of step 1c
               b - combine outputs of step 1c and 2a
               c - obtain a random subsample
            3) a - enumerate resonance forms from outputs of step 2c
               b - combine outputs of step 2c and 3a
               c - obtain a random subsample
            4) a - enumerate stereoisomers from outputs of step 3c
               b - combine outputs of step 3c and 4a
               c - obtain a random subsample
            5) a - enumerate Kekule and/or aromatic SMILES from outputs of step 4c

        :param enum_heterocycles: enumerate all derived heterocycles
        :param enum_tautomers: enumerate all tautomeric forms
        :param enum_resonanceforms: enumerate resonance forms
        :param enum_stereoisomers: enumerate stereoisomers after having dropped the input stereochemistry
        :param enum_smiles: enumerate SMILES
        :param max_enum_heterocycles: maximum number of heterocycles to sample per molecule
        :param max_enum_tautomers: maximum number of tautomers to sample per molecule
        :param max_enum_resonanceforms: maximum number of resonance forms to sample per molecule
        :param max_enum_stereoisomers: maximum number of stereoisomers to sample per molecule
        :param max_enum_smiles: maximum number of enumerated SMILES to sample per molecule
        :param random_choice_heterocycles: number of heterocycles randomly sampled from the obtained `max_enum_heterocycles`
        :param random_choice_tautomers: number of tautomers randomly sampled from the obtained `max_enum_tautomers`
        :param random_choice_resonanceforms: number of resonance forms randomly sampled from the obtained `max_enum_resonanceforms`
        :param random_choice_stereoisomers: number of stereoisomers randomly sampled from the obtained `max_enum_stereoisomers`
        :param max_out_heterocycles: maximum number of heterocycles to be sampled per molecule
        :param max_out_tautomers: maximum number of tautomers to be sampled per molecule
        :param max_out_resonanceforms: maximum number of resonance forms to be sampled per molecule
        :param max_out_stereoisomers: maximum number of stereoisomers to be sampled per molecule
        :param max_out_smiles: maximum number of SMILES to sampled per molecule
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
        self.max_enum_heterocycles = max_enum_heterocycles
        self.max_enum_tautomers = max_enum_tautomers
        self.max_enum_resonanceforms = max_enum_resonanceforms
        self.max_enum_stereoisomers = max_enum_stereoisomers if max_enum_stereoisomers is not None else 2**32
        self.max_enum_smiles = max_enum_smiles
        self.random_choice_heterocycles = random_choice_heterocycles
        self.random_choice_tautomers = random_choice_tautomers
        self.random_choice_resonanceforms = random_choice_resonanceforms
        self.random_choice_stereoisomers = random_choice_stereoisomers
        self.max_out_heterocycles = max_out_heterocycles
        self.max_out_tautomers = max_out_tautomers
        self.max_out_resonanceforms = max_out_resonanceforms
        self.max_out_stereoisomers = max_out_stereoisomers
        self.max_out_smiles = max_out_smiles
        self.smiles_type = smiles_type
        # Generator
        self.random_state = random_state

    def enumerate(self, mols: List[Chem.Mol], njobs: int = 1, chunksize: int = 1) -> List[str]:
        """Enumerate SMILES.

        :param mols: molecules to enumerate the SMILES of.
        :param njobs: number of parallel jobs
        :param chunksize: number of molecules given to each subprocess to enumerate
        :return: the list of enumerated SMILES
        """
        # Parallelize should need be
        if njobs != 1:
            if njobs < 0:
                njobs += multiprocessing.cpu_count() + 1
            with BoundedProcessPoolExecutor(max_workers=njobs) as worker:
                Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)
                futures = [worker.submit(self._enumerator, mol_list)
                           for mol_list in chunked(mols, chunksize)]
            # Collect results
            result = itertools.chain.from_iterable(future.result() for future in futures)
        else:
            # Single process
            result = itertools.chain.from_iterable(self._enumerator(mol) for mol in mols)
        return list(result)

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

    def _enumerator(self, mol: Chem.Mol | List[Chem.Mol]) -> List[str]:
        """Obtain an iterator of SMILES for one single molecule.

        :param mol: Molecule to iterate the SMILES of.
        :return: A generator of SMILES
        """
        block = BlockLogs()
        # Generators
        tauto_enumerator = Chem.MolStandardize.rdMolStandardize.TautomerEnumerator()
        random = np.random.default_rng(self.random_state)
        # Transform input
        mols = mol if isinstance(mol, list) else [mol]
        # Start enumerating
        if self.enum_heterocycles:
            # Obtain heterocycles
            mols = list(itertools.chain(mols,
                                        [self._copy_props(mol, self._sanitize(het))
                                         for mol in mols
                                         for het in map(partial(trycatch, lambda x: x),
                                                        itertools.islice(EnumerateHeterocycles(mol),
                                                                         self.max_enum_heterocycles)
                                                        )
                                         if mol is not None and het is not None]))
            # Subsample heterocycles
            if self.random_choice_heterocycles:
                # Randomly
                n = min(len(mols), self.max_out_heterocycles) if self.max_out_heterocycles is not None else len(mols)
                mols = random.choice(mols, n, replace=False).tolist()
            elif self.max_out_heterocycles is not None:
                # First n
                mols = mols[:self.max_out_heterocycles]
        if self.enum_tautomers:
            # Obtain tautomers
            mols = list(itertools.chain(mols,
                                        [self._copy_props(mol, self._sanitize(tauto))
                                         for mol in mols
                                         for tauto in map(partial(trycatch, lambda x: x),
                                                          itertools.islice(tauto_enumerator.Enumerate(mol),
                                                                           self.max_enum_tautomers)
                                                          )
                                         if mol is not None and tauto is not None and self._sanitize(tauto) is not None]))
            # Subsample tautomers
            if self.random_choice_tautomers:
                # Randomly
                n = min(len(mols), self.max_out_tautomers) if self.max_out_tautomers is not None else len(mols)
                mols = random.choice(mols,  n, replace=False).tolist()
            elif self.max_out_tautomers is not None:
                # First n
                mols = mols[:self.max_out_tautomers]
        if self.enum_resonanceforms:
            # Obtain resonance forms
            mols = list(itertools.chain(mols,
                                        [self._copy_props(mol, self._sanitize(res))
                                         for mol in mols
                                         for res in map(partial(trycatch, lambda x: x),
                                                        itertools.islice(Chem.ResonanceMolSupplier(mol),
                                                                     self.max_enum_resonanceforms)
                                                        )
                                         if mol is not None and res is not None]))
            # Subsample resonance forms
            if self.random_choice_resonanceforms:
                # Randomly
                n = min(len(mols), self.max_out_resonanceforms) if self.max_out_resonanceforms is not None else len(mols)
                mols = random.choice(mols,  n, replace=False).tolist()
            elif self.max_out_resonanceforms is not None:
                # First n
                mols = mols[:self.max_out_resonanceforms]
        if self.enum_stereoisomers:
            # Obtain stereoisomers
            stereo_options= StereoEnumerationOptions(onlyUnassigned=False, maxIsomers=self.max_enum_stereoisomers)
            mols = list(itertools.chain(mols,
                                        [self._copy_props(mol, self._sanitize(stereoi))
                                         for mol in mols
                                         for stereoi in map(partial(trycatch, lambda x: x),
                                                            EnumerateStereoisomers(mol, stereo_options)
                                                            )
                                         if mol is not None and stereoi is not None]))
            # Subsample stereoisomers
            if self.random_choice_stereoisomers:
                # Randomly
                n = min(len(mols), self.max_out_stereoisomers) if self.max_out_stereoisomers is not None else len(mols)
                mols = random.choice(mols,  n, replace=False).tolist()
            elif self.max_out_stereoisomers is not None:
                mols = mols[:self.max_out_stereoisomers]
        # Obtain enumerated SMILES
        if self.enum_smiles:
            # Obtain both Kekule and aromatic SMILES
            if self.smiles_type == 'both':
                mols1 = set([(smiles, mol.GetPropsAsDict(includePrivate=True).get('_Name', ''))
                             for mol in mols
                             for smiles in list(set(trycatch(Chem.MolToSmiles, mol,
                                                             doRandom=True,
                                                             canonical=False,
                                                             kekuleSmiles=True)
                                                    for _ in range(self.max_enum_smiles)
                                                    ))
                             if mol is not None and smiles is not None
                             ])
                # Subsample Kekule SMILES
                if self.max_out_smiles is not None:
                    mols1 = random.choice(list(mols1),  min(len(mols), int(self.max_out_smiles / 2)), replace=False).tolist()
                mols2 = set([(smiles, mol.GetPropsAsDict(includePrivate=True).get('_Name', ''))
                             for mol in mols
                             for smiles in list(set(trycatch(Chem.MolToSmiles, mol,
                                                             doRandom=True,
                                                             canonical=False,
                                                             kekuleSmiles=False)
                                                    for _ in range(self.max_enum_smiles)
                                                    ))
                             if mol is not None and smiles is not None
                             ])
                # Subsample aromatic SMILES
                if self.max_out_smiles is not None:
                    mols2 = random.choice(list(mols2),  min(len(mols), int(self.max_out_smiles / 2)), replace=False).tolist()
                mols = list(itertools.chain(mols1, mols2))
            # Obtain only Kekule SMILES
            elif self.smiles_type == 'kekule':
                mols = set([(smiles, mol.GetPropsAsDict(includePrivate=True).get('_Name', ''))
                            for mol in mols
                            for smiles in list(set(trycatch(Chem.MolToSmiles, mol,
                                                            doRandom=True,
                                                            canonical=False,
                                                            kekuleSmiles=True)
                                                   for _ in range(self.max_enum_smiles)
                                                   ))
                            if mol is not None and smiles is not None
                            ])
                # Subsample Kekule SMILES
                if self.max_out_smiles is not None:
                    mols = random.choice(list(mols),  min(len(mols), self.max_out_smiles), replace=False).tolist()
            # Obtain only aromatic SMILES
            else:
                mols = set([(smiles, mol.GetPropsAsDict(includePrivate=True).get('_Name', ''))
                            for mol in mols
                            for smiles in list(set(trycatch(Chem.MolToSmiles, mol,
                                                            doRandom=True,
                                                            canonical=False,
                                                            kekuleSmiles=False)
                                                   for _ in range(self.max_enum_smiles)
                                                   ))
                            if mol is not None and smiles is not None
                            ])
                # Subsample aromatic SMILES
                if self.max_out_smiles is not None:
                    mols = random.choice(list(mols),  min(len(mols), self.max_out_smiles), replace=False).tolist()
        else:
            mols = [(Chem.MolToSmiles(mol), mol.GetPropsAsDict(includePrivate=True).get('_Name', ''))
                    for mol in mols
                    ]
        del block
        return mols

    def _enumeration_job(self, chunk, queue: multiprocessing.Queue):
        results = self._enumerator(chunk)
        queue.put(results)
        return True

    def _write_results(self, queue: multiprocessing.Queue, out_path: str, progress: bool, total: int):
       if progress:
           pbar = tqdm(total=total, desc='Enumerating molecules', ncols=90, smoothing=0.0)
       with open(out_path, 'w') as oh:
            while True:
                smiles_list = queue.get()
                if smiles_list == (None, None):
                    return True
                for smiles, mol_id in smiles_list:
                    oh.write(f'{smiles}\t{mol_id}\n')
                if progress:
                    pbar.update()

    def enumerate_to_file(self,
                          molsupplier: Union[MolSupplier, List[Chem.Mol]],
                          out_file: str,
                          max_in_mols: int = None,
                          njobs: int = -1,
                          chunksize: int = 1,
                          progress: bool = True
                          ) -> None:
        """Enumerate all derivatives of SMILES contained in the input file and write them in the output file.

        :param in_file: MolSupplier of the input file containing molecules to enumerate SMILES from.
        :param out_file: output SMILES file to write results to.
        :param max_in_mols: maximum number of molecules of the input file to be enumerated
        :param njobs: number of parallel processes
        :param chunksize: number of molecules given to each subprocess to enumerate
        :param progress: should progress be displayed
        """
        if max_in_mols is not None:
            molsupplier = itertools.islice(molsupplier, max_in_mols)
        if njobs == 1:
            result = itertools.chain.from_iterable(self._enumerator(mol) for mol in molsupplier)
            # Write results to output
            with open(out_file, 'w') as oh:
                pbar = tqdm(result, desc='Writing to disk', ncols=90, smoothing=0.0) if progress else result
                for smiles, mol_id in pbar:
                    oh.write(f'{smiles}\t{mol_id}\n')
        else:
            # Create queue of `njobs` max items (lists of SMILES)
            manager = multiprocessing.Manager()
            queue = manager.Queue(njobs)
            # Start the process writing to disk
            writer = multiprocessing.Process(target=self._write_results, args=(queue,
                                                                               out_file,
                                                                               progress,
                                                                               len(molsupplier)),
                                             daemon=False)
            writer.start()
            # Start the submitting job
            if njobs < 0:
                njobs += multiprocessing.cpu_count()  # 1 process needs to be the writer
            with BoundedProcessPoolExecutor(max_workers=njobs) as worker:
                Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)
                futures = [worker.submit(self._enumeration_job, chunk, queue)
                           for chunk in chunked(molsupplier, chunksize)]
            # Ensure all enumerations also are
            assert all([future.result() for future in futures])
            # Signal writer to stop
            for _ in range(njobs):
                queue.put((None, None))
            writer.join()
            writer.close()
            # Free queue
            #queue._close()


def trycatch(func: Callable, *args, handle: Callable = lambda e : e, **kwargs):
    """Convenience function catching exceptions in list comprehensions.

    From @bryan-head & @starball's StackOverflow answer
    https://stackoverflow.com/questions/1528237/how-to-handle-exceptions-in-a-list-comprehensions
    """
    try:
        return func(*args, **kwargs)
    except Exception as e:
        return None
