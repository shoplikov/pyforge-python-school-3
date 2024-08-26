from fastapi import FastAPI, HTTPException
from typing import List, Dict
from os import getenv
import os
from rdkit import Chem
from src.models import Molecule
from src.models import Molecule as MoleculeModel
from fastapi_sqlalchemy import DBSessionMiddleware, db
from dotenv import load_dotenv
from src.schema import Molecule as SchemaMolecule

load_dotenv('.env')

app = FastAPI()

app.add_middleware(DBSessionMiddleware, db_url=os.environ["DATABASE_URL"])

@app.get("/")
async def get_server():
    return {"server_id": getenv("SERVER_ID", "1")}


@app.get("/molecules/", response_model=List[SchemaMolecule])
async def list_molecules():
    molecules = db.session.query(MoleculeModel).all()
    return molecules


@app.post("/molecules/", response_model=SchemaMolecule)
async def add_molecule(molecule: SchemaMolecule):
    if not Chem.MolFromSmiles(molecule.smiles):
        raise HTTPException(status_code=400, detail="Invalid SMILES")

    db_molecule = MoleculeModel(identifier=molecule.identifier, smiles=molecule.smiles)
    db.session.add(db_molecule)
    db.session.commit()

    return db_molecule
# $###############################

@app.get("/molecules/{identifier}", response_model=SchemaMolecule)
async def get_molecule(identifier: str):
    molecule = db.session.query(MoleculeModel).filter_by(identifier=identifier).first()
    if not molecule:
        raise HTTPException(status_code=404, detail="Molecule not found")
    return molecule


@app.put("/molecules/{identifier}", response_model=SchemaMolecule)
async def update_molecule(identifier: str, molecule: SchemaMolecule):
    db_molecule = db.session.query(MoleculeModel).filter_by(identifier=identifier).first()
    if not db_molecule:
        raise HTTPException(status_code=404, detail="Molecule not found")

    if not Chem.MolFromSmiles(molecule.smiles):
        raise HTTPException(status_code=400, detail="Invalid SMILES")

    db_molecule.smiles = molecule.smiles
    db.session.commit()
    return db_molecule


@app.delete("/molecules/{identifier}")
async def delete_molecule(identifier: str):
    db_molecule = db.session.query(MoleculeModel).filter_by(identifier=identifier).first()
    if not db_molecule:
        raise HTTPException(status_code=404, detail="Molecule not found")

    db.session.delete(db_molecule)
    db.session.commit()
    return {"detail": "Molecule deleted"}


@app.post("/molecules/substructure/")
async def substructure_search(smiles: str):
    query_mol = Chem.MolFromSmiles(smiles)
    if not query_mol:
        raise HTTPException(status_code=400, detail="Invalid SMILES for substructure query")

    molecules = db.session.query(MoleculeModel).all()
    matching_molecules = []

    for mol in molecules:
        db_mol = Chem.MolFromSmiles(mol.smiles)
        if db_mol and db_mol.HasSubstructMatch(query_mol):
            matching_molecules.append(mol)

    return matching_molecules
