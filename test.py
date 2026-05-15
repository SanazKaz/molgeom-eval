from ccdc.io import MoleculeReader
from ccdc.conformer import GeometryAnalyser

def build_analyser():
    analyser = GeometryAnalyser()
    analyser.settings.generalisation = False       # faster — skips generalised CSD searches
    analyser.settings.organometallic_filter = 'Organic'  # drug-like molecules only
    return analyser

# --- try it on a single molecule ---
analyser = build_analyser()

# load the first molecule from an SDF file
mol = MoleculeReader('data/drugflow_extracted/1a2g-A-rec-4jmv-1ly-lig-tt-min-0-pocket10_1a2g-A-rec-4jmv-1ly-lig-tt-min-0.sdf')[0]

# standardise bond types and add hydrogens — recommended before any Mogul analysis
# this ensures the fragment matching against the CSD is as accurate as possible
mol.assign_bond_types(which='unknown')      # assign any untyped bonds
mol.standardise_aromatic_bonds()            # normalise aromatic bond representations
mol.standardise_delocalised_bonds()         # normalise delocalised systems e.g. carboxylates
mol.add_hydrogens()                         # Mogul needs explicit hydrogens for fragment matching

# run the full geometry analysis — returns the molecule with analysis attributes attached
analysed = analyser.analyse_molecule(mol)

# for each geometry feature, filter to those where Mogul found enough CSD hits
# to make a confident judgement (default threshold: >= 15 hits)
# then compute the fraction flagged as unusual
for feature_name, features in [
    ('Bond lengths',  analysed.analysed_bonds),
    ('Bond angles',   analysed.analysed_angles),
    ('Torsions',      analysed.analysed_torsions),
    ('Rings',         analysed.analysed_rings),
]:
    # exclude features without enough CSD data — not enough hits to call unusual
    valid = [f for f in features if f.enough_hits]

    if not valid:
        print(f"{feature_name}: no features with enough CSD hits")
        continue

    # fraction of valid features flagged as geometrically unusual by Mogul
    unusual_fraction = sum(1 for f in valid if f.unusual) / len(valid)
    print(f"{feature_name}: {unusual_fraction:.2f} unusual ({len(valid)} features checked)")