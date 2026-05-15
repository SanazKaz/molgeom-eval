import os
from rdkit import Chem

# use env with rdkit and posebusters! not the csd env!
def sdf_to_mol(directory_path, max_files=100):
    """
    Load RDKit molecules from SDF files in a flat directory.

    Parameters
    ----------
    directory_path : str
    max_files : int
        Cap on number of SDF files to read.

    Returns
    -------
    list of rdkit.Chem.Mol
    """
    sdf_files = sorted([
        os.path.join(directory_path, f)
        for f in os.listdir(directory_path)
        if f.endswith('.sdf')
    ])[:max_files]

    mols = []
    for sdf_file in sdf_files:
        suppl = Chem.SDMolSupplier(sdf_file)
        mols.extend([mol for mol in suppl if mol is not None])
    return mols