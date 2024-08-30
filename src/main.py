from fastapi import FastAPI, HTTPException, Query
from typing import List
from os import getenv
import os
from rdkit import Chem
from src.models import Molecule as MoleculeModel
from fastapi_sqlalchemy import DBSessionMiddleware, db
from dotenv import load_dotenv
from src.schema import Molecule as SchemaMolecule
from src.logger import logger

load_dotenv('.env')

app = FastAPI()

app.add_middleware(DBSessionMiddleware, db_url=os.environ["DATABASE_URL"])

@app.get("/")
async def get_server():
    return {"server_id": getenv("SERVER_ID", "1")}

@app.get("/molecules/", response_model=List[SchemaMolecule])
async def list_molecules(skip: int = Query(0, alias="page", ge=0), limit: int = Query(10, le=100)):
    logger.info(f"Listing molecules with pagination: skip={skip}, limit={limit}")
    molecules = db.session.query(MoleculeModel).offset(skip).limit(limit).all()
    return molecules

@app.post("/molecules/", response_model=SchemaMolecule)
async def add_molecule(molecule: SchemaMolecule):
    logger.info(f"Adding molecule with identifier: {molecule.identifier}")
    if not Chem.MolFromSmiles(molecule.smiles):
        raise HTTPException(status_code=400, detail="Invalid SMILES")

    db_molecule = MoleculeModel(identifier=molecule.identifier, smiles=molecule.smiles)
    db.session.add(db_molecule)
    db.session.commit()

    logger.info(f"Molecule added with identifier: {db_molecule.identifier}")
    return db_molecule

@app.get("/molecules/{identifier}", response_model=SchemaMolecule)
async def get_molecule(identifier: str):
    logger.info(f"Retrieving molecule with identifier: {identifier}")
    molecule = db.session.query(MoleculeModel).filter_by(identifier=identifier).first()
    if not molecule:
        raise HTTPException(status_code=404, detail="Molecule not found")
    return molecule

@app.put("/molecules/{identifier}", response_model=SchemaMolecule)
async def update_molecule(identifier: str, molecule: SchemaMolecule):
    logger.info(f"Updating molecule with identifier: {identifier}")
    db_molecule = db.session.query(MoleculeModel).filter_by(identifier=identifier).first()
    if not db_molecule:
        raise HTTPException(status_code=404, detail="Molecule not found")

    if not Chem.MolFromSmiles(molecule.smiles):
        raise HTTPException(status_code=400, detail="Invalid SMILES")

    db_molecule.smiles = molecule.smiles
    db_molecule.identifier = molecule.identifier
    db.session.commit()
    logger.info(f"Molecule updated with identifier: {db_molecule.identifier}")
    return db_molecule

@app.delete("/molecules/{identifier}")
async def delete_molecule(identifier: str):
    logger.info(f"Deleting molecule with identifier: {identifier}")
    db_molecule = db.session.query(MoleculeModel).filter_by(identifier=identifier).first()
    if not db_molecule:
        raise HTTPException(status_code=404, detail="Molecule not found")

    db.session.delete(db_molecule)
    db.session.commit()
    logger.info(f"Molecule deleted with identifier: {identifier}")
    return {"detail": "Molecule deleted"}

@app.post("/molecules/substructure/")
async def substructure_search(smiles: str):
    logger.info(f"Searching for substructures with SMILES: {smiles}")
    query_mol = Chem.MolFromSmiles(smiles)
    if not query_mol:
        raise HTTPException(status_code=400, detail="Invalid SMILES for substructure query")

    molecules = db.session.query(MoleculeModel).all()
    matching_molecules = []

    for mol in molecules:
        db_mol = Chem.MolFromSmiles(mol.smiles)
        if db_mol and db_mol.HasSubstructMatch(query_mol):
            matching_molecules.append(mol)

    logger.info(f"Found {len(matching_molecules)} matching molecules")
    return matching_molecules
