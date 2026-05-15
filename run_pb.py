import os
from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
from posebusters import PoseBusters

from utils import sdf_to_mol

def run_pose_busters(method_name, set_of_mols):
    """Run PoseBusters on a set of mols"""
    buster = PoseBusters(config='mol_fast') # doesnt use energy ratio.
    results = buster.bust(mol_pred=set_of_mols, full_report=False)
    os.makedirs(os.path.join(os.getcwd(), 'results'), exist_ok=True)
    results.to_csv(f'results/{method_name}_posebusters_results.csv', index=False)
    return results

def count_all_pass(results):
    """Number of rows where all columns are True."""
    return int(results.all(axis=1).sum())

if __name__ == "__main__":
    all_results = {}
    directory_paths = ['data/drugflow_extracted', 'data/pocket2mol_extracted', 'data/targetdiff_extracted']

    for directory_path in directory_paths:
        method_name = directory_path.split('/')[-1].split('_')[0]
        mols = sdf_to_mol(directory_path)
        results = run_pose_busters(method_name, mols)
        all_results[method_name] = {
            'passed': count_all_pass(results),
            'total': len(mols)
        }

    summary_rows = [
        {
            'Method': method,
            'Passed': data['passed'],
            'Total': data['total'],
            'Pass Rate (%)': round(data['passed'] / data['total'] * 100, 2) if data['total'] > 0 else 0
        }
        for method, data in all_results.items()
    ]
    output_overall = pd.DataFrame(summary_rows)
    output_overall.to_csv('results/all_methods_posebusters_summary.csv', index=False)

    print("Summary of PoseBusters Results:")
    for _, row in output_overall.iterrows():
        print(f"{row['Method']}: {int(row['Passed'])} out of {int(row['Total'])} passed all checks")