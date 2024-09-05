import json
import time
import pytest
from rdkit import Chem
from src.schema import Molecule as SchemaMolecule
from src.models import Molecule as MoleculeModel

@pytest.fixture
def create_molecule(db_session):
    def _create(identifier, smiles):
        molecule = MoleculeModel(identifier=identifier, smiles=smiles)
        db_session.add(molecule)
        db_session.commit()
        return molecule
    return _create

def test_add_molecule(client, create_molecule):
    response = client.post("/molecules/", json={"identifier": "mol1", "smiles": "C1=CC=CC=C1"})
    assert response.status_code == 200
    assert response.json()["identifier"] == "mol1"

def test_add_duplicate_molecule(client, create_molecule):
    create_molecule("mol1", "C1=CC=CC=C1")
    response = client.post("/molecules/", json={"identifier": "mol1", "smiles": "C1=CC=CC=C1"})
    assert response.status_code == 400
    assert response.json()["detail"] == "Molecule with this identifier or SMILES already exists"

def test_list_molecules(client, create_molecule):
    create_molecule("mol1", "C1=CC=CC=C1")
    create_molecule("mol2", "CCO")
    response = client.get("/molecules/")
    assert response.status_code == 200
    assert len(response.json()) == 2

def test_get_molecule(client, create_molecule):
    molecule = create_molecule("mol1", "C1=CC=CC=C1")
    response = client.get(f"/molecules/{molecule.identifier}")
    assert response.status_code == 200
    assert response.json()["identifier"] == "mol1"

def test_update_molecule(client, create_molecule):
    create_molecule("mol1", "C1=CC=CC=C1")
    response = client.put("/molecules/mol1", json={"identifier": "mol1", "smiles": "CCO"})
    assert response.status_code == 200
    assert response.json()["smiles"] == "CCO"

def test_delete_molecule(client, create_molecule):
    molecule = create_molecule("mol1", "C1=CC=CC=C1")
    response = client.delete(f"/molecules/{molecule.identifier}")
    assert response.status_code == 200
    assert response.json()["detail"] == "Molecule deleted"
    response = client.get(f"/molecules/{molecule.identifier}")
    assert response.status_code == 404

def test_substructure_search(client, create_molecule):
    create_molecule("mol1", "C1=CC=CC=C1")
    create_molecule("mol2", "CCO")
    response = client.post("/molecules/substructure/", json={"smiles": "C1=CC=CC=C1"})
    assert response.status_code == 200
    assert len(response.json()) == 1
    assert response.json()[0]["identifier"] == "mol1"

def test_cache_functionality(client, create_molecule, redis_client):
    create_molecule("mol1", "C1=CC=CC=C1")
    
    response = client.post("/molecules/substructure/", json={"smiles": "C1=CC=CC=C1"})
    assert response.status_code == 200
    assert "id" in response.json()[0]

    redis_client.flushdb()  
    start_time = time.time()
    response = client.post("/molecules/substructure/", json={"smiles": "C1=CC=CC=C1"})
    elapsed_time = time.time() - start_time
    assert response.status_code == 200
    assert elapsed_time < 1
