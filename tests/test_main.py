from unittest.mock import patch
import pytest
from fastapi.testclient import TestClient
from sqlalchemy import create_engine, StaticPool
from sqlalchemy.orm import sessionmaker
from src.main import app, get_db
from src.main import Base
import fakeredis


DATABASE_URL = "sqlite:///:memory:"

engine = create_engine(
    DATABASE_URL,
    connect_args={"check_same_thread": False},
    poolclass=StaticPool,
)
TestingSessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)

@pytest.fixture(scope="module")
def client():
    with TestClient(app) as client:
        yield client


def override_get_db():
    db = TestingSessionLocal()
    try:
        yield db
    finally:
        db.close()

app.dependency_overrides[get_db] = override_get_db

def setup():
    Base.metadata.create_all(bind=engine)

def teardown():
    Base.metadata.drop_all(bind=engine)

@pytest.fixture(scope="module", autouse=True)
def setup_and_teardown():
    setup()
    yield
    teardown()


def test_add_molecule(client):
    response = client.post("/add-molecule", json={"id": 1, "smiles": "COO"})
    assert response.status_code == 201
    assert response.json() == {"id": 1, "smiles": "COO"}

    # Test invalid SMILE
    response = client.post("/add-molecule", json={"id": 2, "smiles": "something"})
    assert response.status_code == 400
    assert response.json() == {"detail": "Invalid SMILES"}

    # Test if ID already exists
    response = client.post("/add-molecule", json={"id": 1, "smiles": "CO"})
    assert response.status_code == 400
    assert response.json() == {"detail": "Molecule ID already exists"}


def test_list_molecules(client):
    response = client.get("/molecule")
    assert response.status_code == 200


def test_get_molecule(client):
    response = client.get("/molecule/1")
    assert response.status_code == 200

    # Test non-existing one
    response = client.get("/molecule/2")
    assert response.status_code == 404


def test_update_molecule(client):
    response = client.put("molecule/update/1", json={"id": 1, "smiles": "CO"})
    assert response.status_code == 200
    assert response.json() == {"id": 1, "smiles": "CO"}

    # Test Invalid SMILE
    response = client.put("molecule/update/1", json={"id": 1, "smiles": "invalid"})
    assert response.status_code == 400
    assert response.json() == {"detail": "Invalid SMILES"}


def test_delete_molecule(client):
    response = client.delete("/molecule/delete/1")
    assert response.status_code == 200
    assert response.json() == {"message": "Molecule deleted successfully"}

    # If ID does not exist
    response = client.delete("/molecule/delete/111")
    assert response.status_code == 404
    assert response.json() == {"detail": "Molecule not found"}


def test_molecule_invalid_id(client):
    response = client.post("/add-molecule", json={"id": -1, "smiles": "CCO"})
    assert response.status_code == 422

    response = client.post("/add-molecule", json={"id": 'str', "smiles": "CCO"})
    assert response.status_code == 422

    response = client.put("molecule/update/-1", json={"id": -1, "smiles": "CO"})
    assert response.status_code == 422


@pytest.fixture
def redis_client():
    # Create a fakeredis instance
    fake_redis = fakeredis.FakeStrictRedis()
    with patch("src.redis.redis_client", fake_redis):
        yield fake_redis

def test_search_substructure_success(redis_client, client):
    with patch("src.tasks.substructure_search_task.delay") as mock_substructure_search_task:
        mock_substructure_search_task.return_value.ready.return_value = True
        mock_substructure_search_task.return_value.state = 'SUCCESS'
        mock_substructure_search_task.return_value.result = [
            {'id': 1, 'smiles': 'CCO'},
            {'id': 2, 'smiles': 'CCN'}
        ]

        response = client.post("/molecule/search", json={"smiles": "C", "limit": 10})

        assert response.status_code == 200
        assert response.json() == {
            "matching_molecules": [{'id': 1, 'smiles': 'CCO'}, {'id': 2, 'smiles': 'CCN'}]
        }
        mock_substructure_search_task.assert_called_once_with("C", 10)

def test_search_substructure_failure(redis_client, client):
    with patch("src.tasks.substructure_search_task.delay") as mock_substructure_search_task:
        mock_substructure_search_task.return_value.ready.return_value = True
        mock_substructure_search_task.return_value.state = 'FAILURE'

        response = client.post("/molecule/search", json={"smiles": "C", "limit": 10})

        assert response.status_code == 400
        assert "Task failed with state: FAILURE" in response.json()["detail"]
