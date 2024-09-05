import pytest
from sqlalchemy.orm import Session
from src.models import Molecule

def test_get_server(client):
    response = client.get("/")
    assert response.status_code == 200
    assert "server_id" in response.json()

def test_add_molecule(client):
    molecule_data = {"identifier": "mol1", "smiles": "CCO"}
    response = client.post("/molecules/", json=molecule_data)
    assert response.status_code == 200
    data = response.json()
    assert data["identifier"] == molecule_data["identifier"]
    assert data["smiles"] == molecule_data["smiles"]

def test_add_duplicate_molecule(client):
    molecule_data = {"identifier": "mol1", "smiles": "CCO"}
    client.post("/molecules/", json=molecule_data)

    response = client.post("/molecules/", json=molecule_data)
    assert response.status_code == 400
    assert response.json() == {"detail": "Molecule with this identifier or SMILES already exists"}

def test_list_molecules(client):
    molecule_data1 = {"identifier": "mol1", "smiles": "CCO"}
    molecule_data2 = {"identifier": "mol2", "smiles": "CCC"}

    client.post("/molecules/", json=molecule_data1)
    client.post("/molecules/", json=molecule_data2)

    response = client.get("/molecules/?page=0&limit=10")
    assert response.status_code == 200
    data = response.json()
    assert len(data) == 2
    assert data[0]["identifier"] == molecule_data1["identifier"]
    assert data[1]["identifier"] == molecule_data2["identifier"]

def test_get_molecule(client):
    molecule_data = {"identifier": "mol1", "smiles": "CCO"}
    client.post("/molecules/", json=molecule_data)

    response = client.get("/molecules/mol1")
    assert response.status_code == 200
    data = response.json()
    assert data["identifier"] == "mol1"
    assert data["smiles"] == "CCO"

def test_update_molecule(client):
    molecule_data = {"identifier": "mol1", "smiles": "CCO"}
    client.post("/molecules/", json=molecule_data)

    updated_data = {"identifier": "mol1", "smiles": "CCC"}
    response = client.put("/molecules/mol1", json=updated_data)
    assert response.status_code == 200
    data = response.json()
    assert data["identifier"] == "mol1"
    assert data["smiles"] == "CCC"

def test_delete_molecule(client):
    molecule_data = {"identifier": "mol1", "smiles": "CCO"}
    client.post("/molecules/", json=molecule_data)

    response = client.delete("/molecules/mol1")
    assert response.status_code == 200
    assert response.json() == {"detail": "Molecule deleted"}

    response = client.get("/molecules/mol1")
    assert response.status_code == 404
