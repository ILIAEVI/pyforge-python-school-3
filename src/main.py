from rdkit import Chem
from fastapi import FastAPI, HTTPException, Depends
from pydantic import BaseModel, Field, ValidationError, field_validator
from os import getenv
from .database import SessionLocal
from sqlalchemy.orm import Session
from .models import Molecule
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

app = FastAPI()

molecules = {}


class MoleculeSchema(BaseModel):
    id: int = Field(..., gt=0)
    smiles: str

    @field_validator('id')
    def validate_id(cls, i):
        if i <= 0:
            raise ValidationError("Id must be greater than 0")
        return i


class SubstructureSearch(BaseModel):
    smiles: str
    limit: int = Field(default=100, gt=0)

    @field_validator('limit')
    def validate_limit(cls, value):
        if value <= 0:
            raise ValueError("Limit must be a positive integer")
        return value


def get_db():
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()


@app.get("/")
def get_server():
    return {"server_id": getenv("SERVER_ID", "1")}


@app.post("/add-molecule")
def add_molecule(molecule: MoleculeSchema, db: Session = Depends(get_db)):
    if not Chem.MolFromSmiles(molecule.smiles):
        raise HTTPException(status_code=400, detail="Invalid SMILES")
    if db.query(Molecule).filter(Molecule.id == molecule.id).first():
        raise HTTPException(status_code=400, detail="Molecule ID already exists")
    new_molecule = Molecule(id=molecule.id, smiles=molecule.smiles)
    db.add(new_molecule)
    db.commit()
    db.refresh(new_molecule)
    logger.info(f"Added molecule with ID: {new_molecule.id} and SMILES: {new_molecule.smiles}")
    return {"id": new_molecule.id, "smiles": new_molecule.smiles}


@app.get("/molecule")
def list_molecules(db: Session = Depends(get_db)):
    molecules_lst = db.query(Molecule).all()
    if molecules_lst:
        return molecules_lst
    else:
        raise HTTPException(status_code=404, detail="No molecules found")


@app.get("/molecule/{molecule_id}")
def get_molecule(molecule_id: int, db: Session = Depends(get_db)):
    molecule = db.query(Molecule).filter(Molecule.id == molecule_id).first()
    if molecule is None:
        raise HTTPException(status_code=404, detail="Molecule not found")
    return molecule


@app.put("/molecule/update/{molecule_id}")
def update_molecule(molecule_id: int, mol: MoleculeSchema, db: Session = Depends(get_db)):
    if not Chem.MolFromSmiles(mol.smiles):
        raise HTTPException(status_code=400, detail="Invalid SMILES")
    molecule = db.query(Molecule).filter(Molecule.id == molecule_id).first()
    if molecule is None:
        raise HTTPException(status_code=404, detail="Molecule not found")
    molecule.smiles = mol.smiles
    db.commit()
    db.refresh(molecule)
    logger.info(f"Updated molecule with ID: {molecule.id} to new SMILES: {molecule.smiles}")
    return {"id": molecule.id, "smiles": molecule.smiles}


@app.delete("/molecule/delete/{molecule_id}")
def delete_molecule(molecule_id: int, db: Session = Depends(get_db)):
    molecule = db.query(Molecule).filter(Molecule.id == molecule_id).first()
    if molecule is None:
        raise HTTPException(status_code=404, detail="Molecule not found")
    db.delete(molecule)
    db.commit()
    logger.info(f"Deleted molecule with ID: {molecule_id}")
    return {"message": "Molecule deleted successfully"}


@app.post("/molecule/search")
def search_substructure(query: SubstructureSearch, db: Session = Depends(get_db)):
    substructure_smiles = query.smiles
    limit = query.limit
    molecules_lst = db.query(Molecule).all()
    if not molecules_lst:
        raise HTTPException(status_code=404, detail="No molecules in the database")
    try:
        matching_molecules = list(substructure_search(molecules_lst, substructure_smiles, limit))
    except Exception as e:
        logger.error(f"Error during substructure search: {e}")
        raise HTTPException(status_code=400, detail=str(e))
    if not matching_molecules:
        raise HTTPException(status_code=404, detail="No substructure matches")
    logger.info(f"Substructure search returned {len(matching_molecules)} molecules")
    return {"matching_molecules": matching_molecules}


def substructure_search(mols, mol, lim):
    matching_molecules = []
    substructure_mol = Chem.MolFromSmiles(mol)
    if substructure_mol is None:
        raise ValueError("Invalid substructure SMILES string")

    for molecule in mols:
        if len(matching_molecules) >= lim:
            break
        molecule_smiles = molecule.smiles
        molecule_mol = Chem.MolFromSmiles(molecule_smiles)
        if molecule_mol is None:
            continue
        if molecule_mol.HasSubstructMatch(substructure_mol):
            matching_molecules.append({'id': molecule.id, 'smiles': molecule.smiles})
            yield {'id': molecule.id, 'smiles': molecule.smiles}


if __name__ == "__main__":
    import uvicorn
    uvicorn.run("main:app", reload=True, host="0.0.0.0", port=8000)
