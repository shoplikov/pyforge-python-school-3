import pytest
from fastapi.testclient import TestClient
from src.main import app
from src.database import engine, Base, SessionLocal
import redis

# Setup Redis for testing
@pytest.fixture(scope='session')
def redis_client():
    client = redis.StrictRedis(host='localhost', port=6379, decode_responses=True)
    yield client
    client.flushdb()

# Setup database for testing
@pytest.fixture(scope='function')
def db_session():
    Base.metadata.create_all(bind=engine)
    session = SessionLocal()
    yield session
    session.close()
    Base.metadata.drop_all(bind=engine)

@pytest.fixture(scope='function')
def client(db_session):
    with TestClient(app) as test_client:
        yield test_client
