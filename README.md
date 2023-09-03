enumerate-smiles
================

A Python library enumerating heterocycles, stereoisomers, tautomers and SMILES.


## 💪 Getting Started

```python
from enumerate_smiles import SmilesEnumerator
from rdkit import Chem

smiles_list = [
  # erlotinib
  "n1cnc(c2cc(c(cc12)OCCOC)OCCOC)Nc1cc(ccc1)C#C",
  # midecamycin
  "CCC(=O)O[C@@H]1CC(=O)O[C@@H](C/C=C/C=C/[C@@H]([C@@H](C[C@@H]([C@@H]([C@H]1OC)O[C@H]2[C@@H]([C@H]([C@@H]([C@H](O2)C)O[C@H]3C[C@@]([C@H]([C@@H](O3)C)OC(=O)CC)(C)O)N(C)C)O)CC=O)C)O)C",
  # selenofolate
  "C1=CC(=CC=C1C(=O)NC(CCC(=O)OCC[Se]C#N)C(=O)O)NCC2=CN=C3C(=N2)C(=O)NC(=N3)N",
]
# Ensure hydrogens are explicit
mols = [Chem.AddHs(Chem.MolFromSmiles(smiles)) for smiles in smiles_list]

se = SmilesEnumerator()
print(se.enumerate(mols))
```

The above enumerates heterocycles, tautomers, resonance forms, stereoismers and SMILES.<br/>

Multiprocessing can be enabled by setting the `njobs` argument. 

One may output the enumerated SMILES to a file as follows:


```python
from enumerate_smiles import SmilesEnumerator
from io import StringIO

output = StringIO()

se = SmilesEnumerator()
se.enumerate_to_file(mols, output, njobs=5)
```

For convenience, users may use a `MolSupplier` instance instead of a list of molecules.

```python
from enumerate_smiles import MolSupplier

input = StringIO("""n1cnc(c2cc(c(cc12)OCCOC)OCCOC)Nc1cc(ccc1)C#C   erlotinib
CCC(=O)O[C@@H]1CC(=O)O[C@@H](C/C=C/C=C/[C@@H]([C@@H](C[C@@H]([C@@H]([C@H]1OC)O[C@H]2[C@@H]([C@H]([C@@H]([C@H](O2)C)O[C@H]3C[C@@]([C@H]([C@@H](O3)C)OC(=O)CC)(C)O)N(C)C)O)CC=O)C)O)C midecamycin
C1=CC(=CC=C1C(=O)NC(CCC(=O)OCC[Se]C#N)C(=O)O)NCC2=CN=C3C(=N2)C(=O)NC(=N3)N  selenofolate
""")
supplier = MolSupplier(input, format='smi')

output = StringIO()

se.enumerate_to_file(supplier, output, njobs=5)
```

## 🚀 Installation

The most recent release can be installed from
[PyPI](https://pypi.org/project/enumerate-smiles/) with:

```shell
$ pip install enumerate-smiles
```

The most recent code and data can be installed directly from GitHub with:

```bash
$ pip install git+https://github.com/OlivierBeq/enumerate-smiles.git
```

## 👐 Contributing

Contributions, whether filing an issue, making a pull request, or forking, are appreciated. See
[CONTRIBUTING.md](https://github.com/OlivierBeq/enumerate-smiles/blob/master/.github/CONTRIBUTING.md) for more information on getting involved.

## 👋 Attribution

Please cite the current repository should you communicate results derived from it.  

### 🍪 Cookiecutter

This package was created with [@audreyfeldroy](https://github.com/audreyfeldroy)'s
[cookiecutter](https://github.com/cookiecutter/cookiecutter) package using [@cthoyt](https://github.com/cthoyt)'s
[cookiecutter-snekpack](https://github.com/cthoyt/cookiecutter-snekpack) template.

