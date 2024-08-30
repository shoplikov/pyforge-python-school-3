# test_main.py
import pytest
from fastapi.testclient import TestClient
from src.main import app
from src.models import Molecule as MoleculeModel

client = TestClient(app)

def test_get_server():
    response = client.get("/")
    assert response.status_code == 200
    assert "server_id" in response.json()

def test_add_molecule():
    molecule_data = {"identifier": "mol", "smiles": "Cc1ccccc1"}
    response = client.post("/molecules/", json=molecule_data)
    assert response.status_code == 201
    assert response.json() == molecule_data

def test_get_molecule():
    molecule_data = {"identifier": "mol", "smiles": "CCO"}
    client.post("/molecules/", json=molecule_data)
    response = client.get("/molecules/mol")
    assert response.status_code == 200
    assert response.json() == molecule_data

def test_update_molecule():
    molecule_data = {"identifier": "mol", "smiles": "Cc1ccccc1"}
    client.post("/molecules/", json=molecule_data)
    updated_data = {"identifier": "mol", "smiles": "CCO"}
    response = client.put("/molecules/mol", json=updated_data)
    assert response.status_code == 200
    assert response.json() == updated_data

def test_delete_molecule():
    molecule_data = {"identifier": "mol", "smiles": "CC(=O)O"}
    client.post("/molecules/", json=molecule_data)
    response = client.delete("/molecules/mol")
    assert response.status_code == 204
    response = client.get("/molecules/mol")
    assert response.status_code == 404

def test_get_all_molecules():
    molecule_data1 = {"identifier": "mol1", "smiles": "Cc1ccccc1"}
    molecule_data2 = {"identifier": "mol2", "smiles": "CCO"}
    client.post("/molecules/", json=molecule_data1)
    client.post("/molecules/", json=molecule_data2)
    response = client.get("/molecules/")
    assert response.status_code == 200
    assert len(response.json()) == 2
    assert molecule_data1 in response.json()
    assert molecule_data2 in response.json()