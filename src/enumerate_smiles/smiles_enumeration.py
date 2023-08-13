# -*- coding: utf-8 -*-

import itertools
import multiprocessing
from typing import List, Callable

import numpy as np
from bounded_pool_executor import BoundedProcessPoolExecutor
from more_itertools import chunked
from rdkit import Chem
from rdkit.Chem import MolStandardize
from rdkit.Chem.EnumerateHeterocycles import EnumerateHeterocycles
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
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
                 max_num_smiles: int = 10,
                 max_enum_thresh: int = 1000,
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
        :param max_enum_thresh: maximum number of variants enumerated that `max_num_smiles` are sampled from
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
        self.max_enum_thresh = max_enum_thresh
        self.smiles_type = smiles_type
        # Generators
        self.random_state = random_state

    def enumerate(self, mols: List[Chem.Mol], njobs: int = -1, chunksize: int = 1) -> List[str]:
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
        # Trnasform input
        mols = mol if isinstance(mol, list) else [mol]
        # Start enumerating
        if self.enum_heterocycles:
            mols = list(itertools.chain(mols,
                                        [self._copy_props(mol, self._sanitize(het))
                                         for mol in mols
                                         for het in itertools.islice(EnumerateHeterocycles(mol),
                                                                     self.max_enum_thresh)
                                         ]))
            mols = random.choice(mols, self.max_num_smiles).tolist()
            # print(f'Enumeration of heterocycles finished ({len(mols)} variants)')
        if self.enum_tautomers:
            mols = list(itertools.chain(mols,
                                        [self._copy_props(mol, self._sanitize(tauto))
                                         for mol in mols
                                         for tauto in itertools.islice(tauto_enumerator.Enumerate(mol),
                                                                       self.max_enum_thresh)
                                         if mol is not None]))
            mols = random.choice(mols, self.max_num_smiles).tolist()
            # print(f'Enumeration of tautomers finished ({len(mols)} variants)')
        if self.enum_resonanceforms:
            mols = list(itertools.chain(mols,
                                        [self._copy_props(mol, self._sanitize(res))
                                         for mol in mols
                                         for res in itertools.islice(Chem.ResonanceMolSupplier(mol),
                                                                     self.max_enum_thresh)
                                         if mol is not None]))
            mols = random.choice(mols, self.max_num_smiles).tolist()
            # print(f'Enumeration of resonance forms finished ({len(mols)} variants)')
        if self.enum_stereoisomers:
            mols = list(itertools.chain(mols,
                                        [self._copy_props(mol, self._sanitize(stereoi))
                                         for mol in mols
                                         for stereoi in itertools.islice(EnumerateStereoisomers(mol),
                                                                         self.max_enum_thresh)
                                         if mol is not None]))
            mols = random.choice(mols, self.max_num_smiles).tolist()
            # print(f'Enumeration of stereoisomers finished ({len(mols)} variants)')
        if self.enum_smiles:
            if self.smiles_type == 'both':
                mols1 = [(smiles, mol.GetPropsAsDict(includePrivate=True).get('_Name', ''))
                        for mol in mols
                        for smiles in list(set(trycatch(Chem.MolToSmiles, mol, doRandom=True, canonical=False, kekuleSmiles=True)
                                               for _ in range(self.max_enum_thresh)
                                               ))
                        if mol is not None and not isinstance(smiles, Exception)
                        ]
                mols1 = random.choice(mols1, int(self.max_num_smiles / 2)).tolist()
                mols2 = [(smiles, mol.GetPropsAsDict(includePrivate=True).get('_Name', ''))
                         for mol in mols
                         for smiles in list(set(trycatch(Chem.MolToSmiles, mol, doRandom=True, canonical=False, kekuleSmiles=False)
                                                for _ in range(self.max_enum_thresh)
                                                ))
                         if mol is not None and not isinstance(smiles, Exception)
                        ]
                mols2 = random.choice(mols2, int(self.max_num_smiles / 2)).tolist()
                mols = list(itertools.chain(mols1, mols2))
                # print(f'Enumeration of SMILES finished ({len(mols)} variants)')
            elif self.smiles_type == 'kekule':
                mols = [(smiles, mol.GetPropsAsDict(includePrivate=True).get('_Name', ''))
                        for mol in mols
                        for smiles in list(set(trycatch(Chem.MolToSmiles, mol, doRandom=True, canonical=False, kekuleSmiles=True)
                                               for _ in range(self.max_enum_thresh)
                                               ))
                        if mol is not None and not isinstance(smiles, Exception)
                        ]
                mols = random.choice(mols, self.max_num_smiles).tolist()
                # print(f'Enumeration of SMILES finished ({len(mols)} variants)')
            else: # aromatic
                mols = [(smiles, mol.GetPropsAsDict(includePrivate=True).get('_Name', ''))
                        for mol in mols
                        for smiles in list(set(trycatch(Chem.MolToSmiles, mol, doRandom=True, canonical=False, kekuleSmiles=False)
                                               for _ in range(self.max_enum_thresh)
                                               ))
                        if mol is not None and not isinstance(smiles, Exception)
                        ]
                mols = random.choice(mols, self.max_num_smiles).tolist()
                # print(f'Enumeration of SMILES finished ({len(mols)} variants)')
        else:
            mols = [(Chem.MolToSmiles(mol), mol.GetPropsAsDict(includePrivate=True).get('_Name', ''))
                    for mol in mols
                    ]
        del block
        return mols

    def _enumeration_job(self, chunk, queue: multiprocessing.Queue):
        print('started')
        results = self._enumerator(chunk)
        print('populating')
        for result in results:
            queue.put(result)
        print('finished')
        return True

    def _write_results(self, queue: multiprocessing.Queue, out_path: str, total: int, progress: bool):
        with open(out_path, 'w') as oh:
            if progress:
                pbar = tqdm(total=total)
            while True:
                smiles, mol_id = queue.get()
                if (smiles, mol_id) == (None, None):
                    if progress:
                        pbar.close()
                    return True
                oh.write(f'{smiles}\t{mol_id}\n')
                if progress:
                    pbar.update()

    def enumerate_to_file(self,
                          molsupplier: MolSupplier,
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
        # Determine the number of molecules
        len_supplier = len(molsupplier)
        if max_in_mols is not None:
            molsupplier = itertools.islice(molsupplier, max_in_mols)
            len_supplier = max_in_mols
        if njobs == 1:
            result = itertools.chain.from_iterable(self._enumerator(mol) for mol in molsupplier)
            # Write results to output
            with open(out_file, 'w') as oh:
                pbar = tqdm(result, total=len_supplier * self.max_num_smiles) if progress else result
                for smiles, mol_id in pbar:
                    oh.write(f'{smiles}\t{mol_id}\n')
        else:
            # Create queue of 1M items
            queue = multiprocessing.Queue(2**30)
            # Start the process writing to disk
            writer = multiprocessing.Process(target=self._write_results, args=(queue,
                                                                               out_file,
                                                                               len_supplier * self.max_num_smiles,
                                                                               progress))
            writer.start()
            # Start the submitting job
            if njobs < 0:
                njobs += multiprocessing.cpu_count()  # 1 process needs to be the writer
            with BoundedProcessPoolExecutor(max_workers=njobs) as worker:
                Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)
                futures = [worker.submit(self._enumeration_job, chunk, queue)
                           for chunk in chunked(molsupplier, chunksize)]
            # Ensure all enumerations also are
            assert all([future.result() for future in futures]) is True
            # Signal writer to stop
            queue.put((None, None))
            writer.join()
            writer.close()
            # Free queue
            queue.close()


def trycatch(func: Callable, *args, handle: Callable = lambda e : e, **kwargs):
    """Convenience function catching exceptions in list comprehensions.

    From @bryan-head & @starball's StackOverflow answer
    https://stackoverflow.com/questions/1528237/how-to-handle-exceptions-in-a-list-comprehensions
    """
    try:
        return func(*args, **kwargs)
    except Exception as e:
        return handle(e)
