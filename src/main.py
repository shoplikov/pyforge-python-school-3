from rdkit import Chem

def substructure_search(mols, mol):
    subs = []
    mols_from_smiles = [Chem.MolFromSmiles(mol) for mol in mols]
    for single_mol in mols_from_smiles:
        if single_mol.HasSubstructMatch(Chem.MolFromSmiles(mol)):
            subs.append(Chem.MolToSmiles(single_mol))
    print(subs)

substructure_search(["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"], "c1ccccc1")