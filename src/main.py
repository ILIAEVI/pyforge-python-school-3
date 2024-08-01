from rdkit import Chem
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
from os import getenv


app = FastAPI()

molecules = {}


class Molecule(BaseModel):
    id: int
    smiles: str


class SubstructureSearch(BaseModel):
    smiles: str


@app.get("/")
def get_server():
    return {"server_id": getenv("SERVER_ID", "1")}


@app.post("/add-molecule")
def add_molecule(molecule: Molecule):
    if not Chem.MolFromSmiles(molecule.smiles):
        raise HTTPException(status_code=400, detail="Invalid SMILES")
    if molecule.id in molecules:
        raise HTTPException(status_code=400, detail="Molecule ID already exists")
    molecules[molecule.id] = molecule.smiles
    return {"id": molecule.id, "smiles": molecule.smiles}


@app.get("/molecule")
def list_molecules():
    if molecules:
        return [Molecule(id=m_id, smiles=smiles) for m_id, smiles in molecules.items()]
    else:
        raise HTTPException(status_code=404, detail="No molecules found")


@app.get("/molecule/{molecule_id}")
def get_molecule(molecule_id: int):
    if molecule_id not in molecules:
        raise HTTPException(status_code=404, detail="Molecule not found")
    return Molecule(id=molecule_id, smiles=molecules[molecule_id])


@app.put("/molecule/update/{molecule_id}")
def update_molecule(molecule_id: int, molecule: Molecule):
    if not Chem.MolFromSmiles(molecule.smiles):
        raise HTTPException(status_code=400, detail="Invalid SMILES")
    if molecule_id not in molecules:
        raise HTTPException(status_code=404, detail="Molecule not found")
    molecules[molecule_id] = molecule.smiles
    return {"id": molecule.id, "smiles": molecule.smiles}


@app.delete("/molecule/delete/{molecule_id}")
def delete_molecule(molecule_id: int):
    if molecule_id not in molecules:
        raise HTTPException(status_code=404, detail="Molecule not found")
    del molecules[molecule_id]
    return {"message": "Molecule deleted successfully"}


@app.post("/molecule/search")
def search_substructure(query: SubstructureSearch):
    substructure_smiles = query.smiles
    try:
        matching_molecules = substructure_search(molecules, substructure_smiles)
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))
    if not matching_molecules:
        raise HTTPException(status_code=404, detail="No substructure matches")
    return {"matching_molecules": matching_molecules}


def substructure_search(mols, mol):
    matching_molecules = []

    substructure_mol = Chem.MolFromSmiles(mol)
    if substructure_mol is None:
        raise ValueError("Invalid substructure SMILES string")

    for mol_id, smiles in mols.items():
        molecule = Chem.MolFromSmiles(smiles)
        if not molecule:
            raise ValueError(f"Invalid SMILES string {smiles}")
        if molecule.HasSubstructMatch(substructure_mol):
            matching_molecules.append({'id': mol_id, 'smiles': smiles})
    return matching_molecules


if __name__ == "__main__":
    import uvicorn
    uvicorn.run("main:app", reload=True, host="0.0.0.0", port=8000)
