import os
import pandas as pd
from tqdm import tqdm
from ccdc.io import MoleculeReader
from ccdc.conformer import GeometryAnalyser

# Self-contained script — ccdc env does not have rdkit/posebusters
# use ccdc env to run this, and run_pb.py in rdkit env


def get_sdf_files(directory_path, max_files=100):
    """Recursively collect SDF file paths from a directory, up to max_files."""
    sdf_files = []
    for root, dirs, files in os.walk(directory_path):
        for file in files:
            if file.endswith('.sdf'):
                sdf_files.append(os.path.join(root, file))
                if len(sdf_files) >= max_files:
                    return sdf_files
    return sdf_files


def build_analyser():
    """Configure GeometryAnalyser for speed: organic-only, no generalisation."""
    analyser = GeometryAnalyser()
    analyser.settings.generalisation = False
    analyser.settings.organometallic_filter = 'Organic'
    return analyser


def compute_unusual_fractions(mol, analyser):
    """
    Compute fraction of unusual bond lengths, angles, torsions, and ring conformations.
    Only includes features with enough CSD hits. Returns None on failure.
    """
    try:
        mol.assign_bond_types(which='unknown')
        mol.standardise_aromatic_bonds()
        mol.standardise_delocalised_bonds()
        mol.add_hydrogens()
        analysed = analyser.analyse_molecule(mol)
    except Exception:
        return None

    results = {}
    for feature_name, features in [
        ('bond_lengths', analysed.analysed_bonds),
        ('bond_angles', analysed.analysed_angles),
        ('torsions', analysed.analysed_torsions),
        ('rings', analysed.analysed_rings),
    ]:
        valid = [f for f in features if f.enough_hits]
        results[feature_name] = sum(1 for f in valid if f.unusual) / len(valid) if valid else None

    return results


def analyse_directory(directory_path, analyser, max_files=100):
    """Run Mogul analysis on SDF files in a directory. Returns DataFrame of fractions."""
    records = []
    sdf_files = get_sdf_files(directory_path, max_files)
    for sdf_path in tqdm(sdf_files, desc=f'{os.path.basename(directory_path)}', unit='file'):
        for mol in MoleculeReader(sdf_path):
            result = compute_unusual_fractions(mol, analyser)
            if result is not None:
                result['mol_name'] = mol.identifier
                result['source_file'] = os.path.basename(sdf_path)
                records.append(result)
    return pd.DataFrame(records)


def summarise_mogul(df, method_name):
    """Return a summary dict of mean unusual fractions for a method."""
    means = df.mean(numeric_only=True)
    return {
        'Method': method_name,
        'Bond Lengths (mean unusual)': round(means.get('bond_lengths', float('nan')), 4),
        'Bond Angles (mean unusual)': round(means.get('bond_angles', float('nan')), 4),
        'Torsions (mean unusual)': round(means.get('torsions', float('nan')), 4),
        'Rings (mean unusual)': round(means.get('rings', float('nan')), 4),
    }


if __name__ == "__main__":
    directory_paths = ['drugflow_extracted', 'pocket2mol_extracted', 'targetdiff_extracted']
    analyser = build_analyser()

    os.makedirs('results', exist_ok=True)

    summary_rows = []
    for directory_path in tqdm(directory_paths, desc='Methods', unit='method'):
        method_name = directory_path.split('_')[0]
        df = analyse_directory(directory_path, analyser, max_files=100)
        df.to_csv(f'results/{method_name}_mogul_results.csv', index=False)
        summary_rows.append(summarise_mogul(df, method_name))

    output_overall = pd.DataFrame(summary_rows)
    output_overall.to_csv('results/all_methods_mogul_summary_generalisation.csv', index=False)

    print("Summary of Mogul Geometry Analysis:")
    for _, row in output_overall.iterrows():
        print(
            f"{row['Method']}: "
            f"bonds={row['Bond Lengths (mean unusual)']:.3f}, "
            f"angles={row['Bond Angles (mean unusual)']:.3f}, "
            f"torsions={row['Torsions (mean unusual)']:.3f}, "
            f"rings={row['Rings (mean unusual)']:.3f}"
        )