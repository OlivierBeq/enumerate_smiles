# -*- coding: utf-8 -*-

import bz2
import gzip
import io
import lzma
import re
from typing import Iterable, Optional, Tuple, Callable, Union

from rdkit import Chem
from rdkit.rdBase import BlockLogs
from rdkit.Chem import (ForwardSDMolSupplier, MaeMolSupplier, MolFromMol2Block, SmilesMolSupplierFromText,
                        SmilesMolSupplier)
from tqdm.auto import tqdm


class ForwardMol2MolSupplier:
    def __init__(self, fileobj: Union[str, io.TextIOBase],
                 sanitize: bool = True,
                 removeHs: bool = True,
                 cleanupSubstructures: bool = True):
        self.sanitize = sanitize
        self.removeHs = removeHs
        self.cleanupSubstructures = cleanupSubstructures
        self._buffer_size = 32768  # 32kB
        self._buffer = b''
        self._mol_delimiter = b'@<TRIPOS>MOLECULE'
        if isinstance(fileobj, str):
            self._open_supplier = True
            self.supplier = open(fileobj, 'rb')
        else:
            self._open_supplier = False
            self.supplier = fileobj

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def _iterate(self):
        self._buffer = self.supplier.read(self._buffer_size)
        while True:
            i_seps = [x.start() for x in re.finditer(self._mol_delimiter, self._buffer) if x.start() != 0]
            if not i_seps:
                new_buffer = self.supplier.read(self._buffer_size)
                if len(new_buffer):
                    self._buffer += new_buffer
                else:
                    mol = MolFromMol2Block(self._buffer,
                                           self.sanitize,
                                           self.removeHs,
                                           self.cleanupSubstructures)
                    yield mol
                    break
            else:
                mol = MolFromMol2Block(self._buffer[:i_seps[0]])
                yield mol
                self._buffer = self._buffer[i_seps[0]:]
                del i_seps[0]

    def __iter__(self):
        if not hasattr(self, '_iterator'):
            self._iterator = self._iterate()
        for values in self._iterator:
            yield values

    def __next__(self):
        if not hasattr(self, '_iterator'):
            self._iterator = self._iterate()
        return next(self._iterator)

    def close(self):
        if self._open_supplier:
            self.supplier.close()


class ForwardSmilesMolSupplier:
    def __init__(self, fileobj: Union[str, io.TextIOBase],
                 delimiter: str = '\t',
                 smilesColumn: int = 0,
                 nameColumn: int = 1,
                 titleLine: bool = True,
                 sanitize: bool = True):
        self.delimiter = delimiter
        self.smilesColumn = smilesColumn
        self.nameColumn = nameColumn
        self.titleLine = titleLine
        self.sanitize = sanitize
        self._buffer_size = 32768  # 32kB
        self._buffer = b''
        self._mol_delimiter = b'\n'
        if isinstance(fileobj, str):
            self._open_supplier = True
            self.supplier = None
            self._iterator = SmilesMolSupplier(fileobj, self.delimiter, self.smilesColumn, self.nameColumn,
                                               self.titleLine, self.sanitize)
        else:
            self._open_supplier = False
            self.supplier = fileobj

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def _iterate(self):
        if self.titleLine:
            self.supplier.readline()
        self._buffer = self.supplier.read(self._buffer_size)
        while True:
            i_seps = [x.start() for x in re.finditer(self._mol_delimiter, self._buffer)]
            if not i_seps:
                new_buffer = self.supplier.read(self._buffer_size)
                if len(new_buffer):
                    self._buffer += new_buffer
                else:
                    if len(self._buffer):
                        block = BlockLogs() # Disable logger if no name column
                        mol = next(SmilesMolSupplierFromText(self._buffer, self.delimiter, self.smilesColumn,
                                                             self.nameColumn, False, self.sanitize))
                        del block  # Disable logger if no name column
                        yield mol
                    break
            else:
                block = BlockLogs() # Disable logger if no name column
                mol = next(
                    SmilesMolSupplierFromText(self._buffer[:i_seps[0] + len(self._mol_delimiter)], self.delimiter,
                                              self.smilesColumn, self.nameColumn, False, self.sanitize))
                del block  # Disable logger if no name column
                yield mol
                self._buffer = self._buffer[i_seps[0] + len(self._mol_delimiter):]
                del i_seps[0]

    def __iter__(self):
        if not hasattr(self, '_iterator'):
            self._iterator = self._iterate()
        for values in self._iterator:
            yield values

    def __next__(self):
        if not hasattr(self, '_iterator'):
            self._iterator = self._iterate()
        return next(self._iterator)

    def close(self):
        if self._open_supplier:
            del self._iterator
            self._iterator = None


class MolSupplier:
    # class properties
    valid_formats = ('smi', 'mae', 'sd', 'sdf', 'mol', 'mol2')
    valid_compression = ('lzma', 'zlib', 'bz2')

    def __init__(self, source: Union[str, io.TextIOBase, io.BufferedIOBase],
                 format: str = None,
                 compression: str = None, **kwargs):
        f"""Molecular supplier handling format and compression.
        :param source: filename or file-like object;
                       when using a context manager, file-like objects
                       are not closed upon exit
        :param format: data format {self.valid_formats}
                       can be detected if source is a file name,
                       must be provided if source is a not file name,
                       ignored if supplier is not None
        :param compression: compression type {self.valid_compression}
                            can be detected if source is a file name,
                            ignored otherwise
        :param kwargs: keyworded arguments to be passed to the underlying supplier,
                       ignored if source is a supplier
                       can also hold values for 'total' and 'show_progress'
                       to be considered when used as an iterable
        """
        # source is None
        if source is None:
            raise ValueError('source must be supplied')
        # Default attributes
        self._open_substream = False  # should a file be opened
        self.filename = None  # name of file to be opened
        self.open_fn = None  # function opening file and handling compression
        self._handle = None  # handle to opened file
        self.supplier = None  # molecule supplier
        self.compression = None
        self.format = None
        self.kwargs = kwargs  # additional parameters for suppliers
        self._iter_total = self.kwargs.pop('total', None)
        self._iter_progress = self.kwargs.pop('show_progress', None)
        # source is a file name
        if isinstance(source, str):
            self.filename = source
            self._open_substream = True
            # Handle compressions
            if compression is not None:
                if compression not in self.valid_compression:
                    raise ValueError(f'compression must be one of {self.valid_compression}')
                self.compression = compression
            else:
                self.compression, self._trunc_filename = self._get_compression(self.filename)
            self.open_fn = self._get_compression_handler(self.compression)
            # Handle file types
            if format is not None:
                if format not in self.valid_formats:
                    raise ValueError(f'format must be one of {self.valid_formats}')
                self.format = format
            else:
                self.format = self._get_format(self._trunc_filename)
        # source is file-like object
        elif isinstance(source, (io.TextIOBase, io.BufferedIOBase)):
            if format is None:
                raise ValueError('format must be specified with text or binary readers')
            self._handle = source
            self.format = format
        else:
            raise ValueError('source must either be filename or file-like object')
        # Create rdkit suppliers
        if self._open_substream:
            self._handle = self.open_fn(self.filename, 'rb')
        # Open a supplier
        if self.format == 'smi':
            self.supplier = ForwardSmilesMolSupplier(self._handle, **self.kwargs)
        elif self.format == 'mae':
            self.supplier = MaeMolSupplier(self._handle, **self.kwargs)
        elif self.format in ['sd', 'sdf', 'mol']:
            self.supplier = ForwardSDMolSupplier(self._handle, **self.kwargs)
        elif self.format == 'mol2':
            self.supplier = ForwardMol2MolSupplier(self._handle, **self.kwargs)

    def set_start_progress_total(self, start: int=0, progress: bool=True, total: Optional[int]=None):
        """Set the start, progress and total for iterating through the supplier.
        :param start: starting value for generated identifiers while enumerating molecules
        :param progress: whether a progress bar should be displayed
        :param total: total number of molecules in the supplier
        """
        self._iter_start = start
        self._iter_total = total
        self._iter_progress = progress

    @staticmethod
    def _get_compression(filename: str) -> Tuple[Optional[str], str]:
        """Get compression type and stripped filename."""
        if filename.endswith('.xz'):
            return 'lzma', filename.rstrip('.xz')
        elif filename.endswith('.gz'):
            return 'zlib', filename.rstrip('.gz')
        elif filename.endswith('.bz2'):
            return 'bz2', filename.rstrip('.bz2')
        else:
            return None, filename

    @staticmethod
    def _get_compression_handler(compression_type) -> Callable:
        """Get function to deal with the compression."""
        if compression_type == 'lzma':
            return lzma.open
        elif compression_type == 'zlib':
            return gzip.open
        elif compression_type == 'bz2':
            return bz2.open
        elif compression_type is None:
            return open
        else:
            raise ValueError(f'type compression not handled: {compression_type}')

    @staticmethod
    def _get_format(filename) -> str:
        """Get file format from filename."""
        if filename.endswith('.smi'):
            return 'smi'
        elif filename.endswith('.mae'):
            return 'mae'
        elif filename.endswith(('.sd', '.sdf')):
            return 'sd'
        elif filename.endswith('.mol2'):
            return 'mol2'
        elif filename.endswith('.mol'):
            return 'mol'

    def _processed_mol_supplier(self) -> Iterable[Tuple[int, Chem.Mol]]:
        """Generator function that reads from a rdkit molecule supplier."""
        # handle showing progress
        if self._iter_progress:
            pbar = tqdm(self.supplier, total=self._iter_total, ncols=100)
        else:
            pbar = self.supplier
        for rdmol in pbar:
            if rdmol:
                yield rdmol
            else:
                continue

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def __iter__(self):
        if not hasattr(self, '_iterator'):
            self._iterator = self._processed_mol_supplier()
        for values in self._iterator:
            yield values

    def __next__(self):
        if not hasattr(self, '_iterator'):
            self._iterator = self._processed_mol_supplier()
        return next(self._iterator)

    def close(self):
        del self.supplier
        self.supplier = None
        if self._open_substream:
            self._handle.close()

    def reset(self):
        """Reset the supplier's cursor to the beginning of the file."""
        self._handle.seek(0, 0)
        if self.format == 'smi':
            self.supplier = ForwardSmilesMolSupplier(self._handle, **self.kwargs)
        elif self.format == 'mae':
            self.supplier = MaeMolSupplier(self._handle, **self.kwargs)
        elif self.format in ['sd', 'sdf', 'mol']:
            self.supplier = ForwardSDMolSupplier(self._handle, **self.kwargs)
        elif self.format == 'mol2':
            self.supplier = ForwardMol2MolSupplier(self._handle, **self.kwargs)
        # Reset iterator
        del self._iterator

    def __len__(self):
        return self.size

    @property
    def size(self):
        return self.get_size()

    def get_size(self):
        """Obtain the number of molecules."""
        if hasattr(self, '_size'):
            return self._size
        self._size = sum(1 for mol in self if mol is not None)
        self.reset()
        return self._size
