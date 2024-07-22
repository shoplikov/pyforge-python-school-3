from rdkit import Chem

def substructure_search(mols, mol):
    subs = []
    molecules = [Chem.MolFromSmiles(smile) for smile in mols]
    converted_mol = Chem.MolFromSmiles(mol)
    for molecule in molecules:
        if molecule.HasSubstructMatch(converted_mol):
            subs.append(Chem.MolToSmiles(molecule))
    return subs

print(substructure_search(["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"], "c1ccccc1"))