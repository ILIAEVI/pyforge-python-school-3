from fastapi.testclient import TestClient
from src.main import app

client = TestClient(app)


def test_add_molecule():
    response = client.post("/add-molecule", json={"id": 1, "smiles": "CCO"})
    assert response.status_code == 200
    assert response.json() == {"id": 1, "smiles": "CCO"}

    # Test invalid SMILE
    response = client.post("/add-molecule", json={"id": 2, "smiles": "something"})
    assert response.status_code == 400
    assert response.json() == {"detail": "Invalid SMILES"}

    # Test if ID already exists
    response = client.post("/add-molecule", json={"id": 1, "smiles": "CO"})
    assert response.status_code == 400
    assert response.json() == {"detail": "Molecule ID already exists"}


def test_list_molecules():
    response = client.get("/molecule")
    assert response.status_code == 200
    assert len(response.json()) == 1


def test_get_molecule():
    response = client.get("/molecule/1")
    assert response.status_code == 200
    assert response.json() == {"id": 1, "smiles": "CCO"}

    # Test non-existing one
    response = client.get("/molecule/2")
    assert response.status_code == 404
    assert response.json() == {"detail": "Molecule not found"}


def test_update_molecule():
    response = client.put("molecule/update/1", json={"id": 1, "smiles": "CO"})
    assert response.status_code == 200
    assert response.json() == {"id": 1, "smiles": "CO"}

    # Test Invalid SMILE
    response = client.put("molecule/update/1", json={"id": 1, "smiles": "invalid"})
    assert response.status_code == 400
    assert response.json() == {"detail": "Invalid SMILES"}

    # If ID does not exist
    response = client.put("molecule/update/2", json={"id": 2, "smiles": "CO"})
    assert response.status_code == 404
    assert response.json() == {"detail": "Molecule not found"}


def test_delete_molecule():
    response = client.delete("/molecule/delete/1")
    assert response.status_code == 200
    assert response.json() == {"message": "Molecule deleted successfully"}

    # If ID does not exist
    response = client.delete("/molecule/delete/111")
    assert response.status_code == 404
    assert response.json() == {"detail": "Molecule not found"}


def test_search_substructure():
    client.post("/add-molecule", json={"id": 2, "smiles": "CCOCC"})
    client.post("/add-molecule", json={"id": 3, "smiles": "CCNCC"})

    response = client.post("/molecule/search", json={"smiles": "CC"})
    assert response.status_code == 200
    assert len(response.json()["matching_molecules"]) == 2

    # Test case where molecule does not match
    response = client.post("/molecule/search", json={"smiles": "CCCNO"})
    assert response.status_code == 404
    assert response.json() == {"detail": "No substructure matches"}

    # Test Invalid SMILE
    response = client.post("/molecule/search", json={"smiles": "Invalid"})
    assert response.status_code == 400
    assert response.json() == {"detail": "Invalid substructure SMILES string"}


def test_molecule_invalid_id():
    response = client.post("/add-molecule", json={"id": -1, "smiles": "CCO"})
    assert response.status_code == 422

    response = client.post("/add-molecule", json={"id": 'str', "smiles": "CCO"})
    assert response.status_code == 422

    response = client.put("molecule/update/-1", json={"id": -1, "smiles": "CO"})
    assert response.status_code == 422