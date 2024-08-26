import pytest
from fastapi.testclient import TestClient
from src.main import app
from src.models import Molecule as MoleculeModel

client = TestClient(app)

@pytest.fixture(autouse=True)
def setup_function():
    # Clear existing molecules before each test
    with client.get("/molecules/"):
        response = client.get("/molecules/")
        for molecule in response.json():
            client.delete(f"/molecules/{molecule['identifier']}")

def test_add_molecule():
    response = client.post("/molecules/", json={"identifier": "mol", "smiles": "Cc1ccccc1"})
    assert response.status_code == 200
    assert response.json() == {"identifier": "mol", "smiles": "Cc1ccccc1"}

def test_get_molecule():
    client.post("/molecules/", json={"identifier": "mol", "smiles": "CCO"})
    response = client.get("/molecules/mol")
    assert response.status_code == 200
    assert response.json() == {"identifier": "mol", "smiles": "CCO"}

def test_update_molecule():
    client.post("/molecules/", json={"identifier": "mol", "smiles": "Cc1ccccc1"})
    response = client.put("/molecules/mol", json={"identifier": "mol", "smiles": "CCO"})
    assert response.status_code == 200
    assert response.json() == {"identifier": "mol", "smiles": "CCO"}
    response = client.get("/molecules/mol")
    assert response.json() == {"identifier": "mol", "smiles": "CCO"}

def test_delete_molecule():
    client.post("/molecules/", json={"identifier": "mol", "smiles": "CC(=O)O"})
    response = client.delete("/molecules/mol")
    assert response.status_code == 200
    assert response.json() == {"detail": "Molecule deleted"}
    response = client.get("/molecules/mol")
    assert response.status_code == 404

def test_list_molecules():
    client.post("/molecules/", json={"identifier": "mol1", "smiles": "c1ccccc1"})
    client.post("/molecules/", json={"identifier": "mol2", "smiles": "CC(=O)O"})
    response = client.get("/molecules/")
    assert response.status_code == 200
    molecules = response.json()
    assert len(molecules) == 2
    assert any(mol['identifier'] == 'mol1' for mol in molecules)
    assert any(mol['identifier'] == 'mol2' for mol in molecules)

def test_substructure_search():
    client.post("/molecules/", json={"identifier": "mol1", "smiles": "CC(=O)O"})
    client.post("/molecules/", json={"identifier": "mol2", "smiles": "c1ccccc1"})
    
    response = client.post('/molecules/substructure/', json="c1ccccc1")
    assert response.status_code == 200
    results = response.json()
    assert len(results) == 1
    assert results[0]['identifier'] == 'mol2'
    assert results[0]['smiles'] == 'c1ccccc1'

def test_get_server():
    response = client.get("/")
    assert response.status_code == 200
    assert "server_id" in response.json()
