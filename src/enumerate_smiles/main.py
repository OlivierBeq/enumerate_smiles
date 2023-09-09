from itertools import islice
from IO import MolSupplier
from smiles_enumeration import SmilesEnumerator
from time import perf_counter

if __name__ == '__main__':
    supplier = MolSupplier(r'C:\Users\ojbeq\.data\jax-cddd\data\Combined_dataset-6.smi')
    se = SmilesEnumerator(max_out_smiles=100,
                          max_out_heterocycles=100,
                          max_out_tautomers=100,
                          max_out_resonanceforms=100,
                          max_out_stereoisomers=100,
                          max_enum_heterocycles=10,
                          max_enum_tautomers=10,
                          max_enum_resonanceforms=10,
                          max_enum_stereoisomers=10,
                          max_enum_smiles=10)
    # se.enumerate_to_file(supplier, 'C:/Users/ojbeq/Downloads/test_smilesenum.smi')
    x = perf_counter()
    se.enumerate(islice(supplier,  10_000))
    y = perf_counter()
    print(y - x)
