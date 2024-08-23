from sqlalchemy import Column, Integer, String
from sqlalchemy.ext.declarative import declarative_base


Base = declarative_base()

class Molecule(Base):
    __tablename__ = 'molecules'
    id = Column(Integer, primary_key=True, index=True)
    smiles = Column(String, index=True)
