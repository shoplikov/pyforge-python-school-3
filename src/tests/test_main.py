import pytest
from fastapi.testclient import TestClient
from main import app, molecule_db

client = TestClient(app)

@pytest.fixture(autouse=True)
def setup_function():
    molecule_db.clear()

def test_add_molecule():
    response = client.post("/molecule/add", json={"identifier": "mol", "smiles": "Cc1ccccc1"})
    assert response.status_code == 200
    assert response.json() == {"message": "molecule was added"}

def test_get_molecule():
    client.post("/molecule/add", json={"identifier": "mol", "smiles": "CCO"})
    response = client.get("/molecule/mol")
    assert response.status_code == 200
    assert response.json() == {"identifier": "mol", "smiles": "CCO"}

def test_update_molecule():
    client.post("/molecule/add", json={"identifier": "mol", "smiles": "Cc1ccccc1"})
    response = client.put("/molecule/mol", json={"identifier": "mol", "smiles": "CCO"})
    assert response.status_code == 200
    assert response.json() == {"message": "molecule was updated"}
    response = client.get("/molecule/mol")
    assert response.json() == {"identifier": "mol", "smiles": "CCO"}

def test_delete_molecule():
    client.post("/molecule/add", json={"identifier": "mol", "smiles": "CC(=O)O"})
    response = client.delete("/molecule/mol")
    assert response.status_code == 200
    assert response.json() == {"message": "molecule was deleted"}
    response = client.get("/molecule/mol1")
    assert response.status_code == 404

def test_list_molecules():
    client.post("/molecule/add", json={"identifier": "mol1", "smiles": "c1ccccc1"})
    client.post("/molecule/add", json={"identifier": "mol2", "smiles": "CC(=O)O"})
    response = client.get("/molecules")
    assert response.status_code == 200
    assert len(response.json()["molecules"]) == 2

def test_substructure_search():
    client.post("/molecule/add", json={"identifier": "mol1", "smiles": "CC(=O)O"})
    client.post("/molecule/add", json={"identifier": "mol2", "smiles": "c1ccccc1"})

    response = client.get("/molecules")
    print("Molecules added:", response.json())
    assert response.status_code == 200
    assert len(response.json()["molecules"]) == 2

    # Perform substructure search
    # response = client.post("/search", json={"mol": "C"})
    response = client.post('/search', params={"mol": "c1ccccc1"}) 
    assert response.status_code == 200
    results = response.json()
    assert len(results) == 1
    assert results[0]['identifier'] == 'mol2'
    assert results[0]['smiles'] == 'c1ccccc1'

def test_get_server():
    response = client.get("/")
    assert response.status_code == 200
    assert "server_id" in response.json()
