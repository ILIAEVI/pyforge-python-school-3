from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, scoped_session


DATABASE_URL = "postgresql://test:test@localhost:5432/test"

engine = create_engine(DATABASE_URL)
SessionLocal = scoped_session(sessionmaker(autocommit=False, autoflush=False, bind=engine))

def init_db():
    import models
    models.Base.metadata.create_all(bind=engine)
