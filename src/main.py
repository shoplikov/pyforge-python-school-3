from rdkit import Chem
from typing import List, Dict, Optional
from fastapi import FastAPI, HTTPException
from .models import Molecule

app = FastAPI()

molecule_db: Dict[str, Molecule] = {}

@app.post('/molecule/add')
async def add_molecule(molecule: Molecule):
    if molecule.identifier in molecule_db:
        raise HTTPException(status_code=400, detail="identifier already exists")
    molecule_db[molecule.identifier] = molecule
    return {"message": "molecule was added"}

@app.get('/molecule/{identifier}')
async def get_molecule(identifier: str):
    molecule = molecule_db.get(identifier)
    if molecule is None:
        raise HTTPException(status_code=404, detail="there's no such molecule")
    return molecule

@app.put('/molecule/{identifier}')
async def update_molecule(identifier: str, molecule: Molecule):
    if identifier not in molecule_db:
        raise HTTPException(status_code=404, detail="there's no such molecule")
    molecule_db[identifier] = molecule
    return {"message": "molecule was updated"}

@app.delete('/molecule/{identifier}')
async def delete_molecule(identifier: str):
    molecule = molecule_db.pop(identifier, None)
    if molecule is None:
        raise HTTPException(status_code=404, detail="there's no such molecule")
    return {"message": "molecule was deleted"}

@app.get('/molecules')
async def list_molecules():
    return {"molecules": list(molecule_db.values())}

@app.post('/search')
async def substructure_search(mol: str) -> List[Molecule]:
    subs = []
    converted_mol = Chem.MolFromSmiles(mol)
    if converted_mol is None:
        raise HTTPException(status_code=400, detail="wrong smile format")
    
    for molecule in molecule_db.values():
        molecules = Chem.MolFromSmiles(molecule.smiles)
        if molecules and molecules.HasSubstructMatch(converted_mol):
            subs.append(molecule)
    return subs