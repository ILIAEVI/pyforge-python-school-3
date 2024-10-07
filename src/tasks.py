from src.celery_worker import celery
from rdkit import Chem
from src.models import Molecule
from sqlalchemy.orm import sessionmaker
from src.database import engine


SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)



@celery.task
def substructure_search_task(substructure_smiles, limit):
    db = None
    try:
        db = SessionLocal()
        molecules_lst = db.query(Molecule).all()
        if not molecules_lst:
            raise ValueError("No molecules in the database")
        try:
            matching_molecules = list(substructure_search(molecules_lst, substructure_smiles, limit))
        except Exception as e:
            raise ValueError(f"Error during substructure search: {e}")
        return matching_molecules
    except Exception as e:
        print(f"Error occurred: {e}")
    finally:
        if db:
            db.close()


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
