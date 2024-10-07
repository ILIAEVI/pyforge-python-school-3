from sqlalchemy.orm import sessionmaker, scoped_session
from sqlalchemy import create_engine


DATABASE_URL = f"postgresql://smiles:smiles@db:5432/smiles"

engine = create_engine(DATABASE_URL)
SessionLocal = scoped_session(sessionmaker(autocommit=False, autoflush=False, bind=engine))


def get_db():
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()
