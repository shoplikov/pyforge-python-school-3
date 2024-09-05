import pytest
import os
from fastapi.testclient import TestClient
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from src.main import app
from src.database import Base, get_db

# Override DATABASE_URL for tests
os.environ["DATABASE_URL"] = "sqlite:///:memory:"

# Use an in-memory SQLite database for tests
engine = create_engine(os.environ["DATABASE_URL"], connect_args={"check_same_thread": False})
TestingSessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)

# Override the get_db dependency to use the test database
def override_get_db():
    db = TestingSessionLocal()
    try:
        yield db
    finally:
        db.close()

app.dependency_overrides[get_db] = override_get_db

@pytest.fixture(scope="function", autouse=True)
def setup_database():
    # Create all the tables
    Base.metadata.create_all(bind=engine)
    yield
    # Drop all the tables after the test
    Base.metadata.drop_all(bind=engine)

@pytest.fixture(scope="module")
def client():
    # Return a test client for making requests
    return TestClient(app)
