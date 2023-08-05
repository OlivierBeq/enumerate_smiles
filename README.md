enumerate-smiles
================

A Python library enumerating heterocycles, stereoisomers, tautomers and SMILES.


## 💪 Getting Started

```python
from chemopy import ChemoPy
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

cmp = ChemoPy()
print(cmp.calculate(mols))
```

The above calculates 632 two-dimensional (2D) molecular descriptors.<br/>

The additional 552 three-dimensional (3D) molecular descriptors can be computed as follows:<br/>
:warning: Molecules are required to have conformers for descriptors to be calculated.

```python
from rdkit.Chem import AllChem

# Ensure molecules have 3D conformations
for mol in mols:
    _ = AllChem.EmbedMolecule(mol)

cmp = ChemoPy(ignore_3D=False)
print(cmp.calculate(mols))
```

To obtain 11 2D molecular fingerprints with default folding size, one can use the following:

```python
from chemopy import Fingerprint

for mol in mols:
    print(Fingerprint.get_all_fps(mol))
```

Other methods of the `Fingerprint` submodule allow to change folding size and depth of the compatible fingerprints.

Currently, only one 3D fingerprint may be obtained with the following:

```python
from chemopy import Fingerprint3D

for mol in mols:
    print(Fingerprint3D.get_all_fps(mol))
```

Conveniently, the calculation of molecular descriptors and fingerprints can be coupled:
```python
# 2D molecular descriptors and fingerprints
cmp =  ChemoPy(include_fps=True)
print(cmp.calculate(mols))

# 2D and 3D molecular descriptors and fingerprints
cmp =  ChemoPy(ignore_3D=False, include_fps=True)
print(cmp.calculate(mols))
```

Finally, details on either one or all molecular descriptors can be obtained like so:

```python
# Obtain details about Thara
print(cmp.get_details('Thara'))

# Obtain details for all descriptors and fingerprints
print(cmp.get_details())
```

## 🚀 Installation

**ChemoPy requires OpenMOPAC and OpenBabel to be installed.<br/>**
```
conda install openbabel mopac -c conda-forge
```

The most recent release can be installed from
[PyPI](https://pypi.org/project/chemopy2/) with:

```shell
$ pip install chemopy2
```

The most recent code and data can be installed directly from GitHub with:

```bash
$ pip install git+https://github.com/OlivierBeq/chemopy.git
```

## 👐 Contributing

Contributions, whether filing an issue, making a pull request, or forking, are appreciated. See
[CONTRIBUTING.md](https://github.com/OlivierBeq/chemopy/blob/master/.github/CONTRIBUTING.md) for more information on getting involved.

## 👋 Attribution

### 📖 Citation

1. Cao et al., ChemoPy: freely available python package for computational biology and chemoinformatics. Bioinformatics 2013; 29(8), 1092–1094. doi:10.1093/bioinformatics/btt105

### 🍪 Cookiecutter

This package was created with [@audreyfeldroy](https://github.com/audreyfeldroy)'s
[cookiecutter](https://github.com/cookiecutter/cookiecutter) package using [@cthoyt](https://github.com/cthoyt)'s
[cookiecutter-snekpack](https://github.com/cthoyt/cookiecutter-snekpack) template.

