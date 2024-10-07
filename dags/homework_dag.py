import pendulum
from airflow import DAG
from airflow.operators.python import PythonOperator
from minio import Minio
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from datetime import datetime, timedelta
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
import os



MINIO_URL = 'localhost:9000'
MINIO_ACCESS_KEY = 'minio'
MINIO_SECRET_KEY = 'minio'
BUCKET_NAME = 'smiles'
minio_client = Minio(MINIO_URL, access_key=MINIO_ACCESS_KEY, secret_key=MINIO_SECRET_KEY, secure=False)

DATABASE_URL = f"postgresql://smiles:smiles@db:5432/smiles"
engine = create_engine(DATABASE_URL)
SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)


def extract_data(**kwargs):
    session = SessionLocal()
    today = datetime.today().date()

    molecules = session.execute(f"SELECT id, smiles FROM molecules WHERE DATE(created_at) = '{today}'").fetchall()
    session.close()

    df = pd.DataFrame(molecules, columns=["id", "smiles"])

    kwargs['ti'].xcom_push(key='extracted_data', value=df.to_dict())


def transform_data(**kwargs):
    extracted_data = kwargs['ti'].xcom_pull(key='extracted_data')
    df = pd.DataFrame.from_dict(extracted_data)

    def compute_properties(smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            mw = rdMolDescriptors.CalcExactMolWt(mol)
            logp = Descriptors.MolLogP(mol)
            tpsa = Descriptors.TPSA(mol)
            h_donors = Descriptors.NumHDonors(mol)
            h_acceptors = Descriptors.NumHAcceptors(mol)
            lipinski = (mw <= 500 and logp <= 5 and h_donors <= 5 and h_acceptors <= 10)
            return mw, logp, tpsa, h_donors, h_acceptors, lipinski
        else:
            return None, None, None, None, None, None

    df[['Molecular_Weight', 'LogP', 'TPSA', 'H_Donors', 'H_Acceptors', 'Lipinski_Pass']] = df['smiles'].apply(
        lambda smiles: pd.Series(compute_properties(smiles)))

    kwargs['ti'].xcom_push(key='transformed_data', value=df.to_dict())



def load_data(**kwargs):
    transformed_data = kwargs['ti'].xcom_pull(key='transformed_data')
    df = pd.DataFrame.from_dict(transformed_data)

    file_name = 'molecule_data.xlsx'
    df.to_excel(file_name, index=False)

    with open(file_name, 'rb') as file_data:
        minio_client.put_object(
            bucket_name=BUCKET_NAME,
            object_name=file_name,
            data=file_data,
            length=os.path.getsize(file_name),
            content_type='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'
        )
    os.remove(file_name)



with DAG(
    dag_id='molecule_etl_pipeline',
    start_date=pendulum.today(),
    schedule=None,
    tags=['homework']
) as dag:
    start_op = extract = PythonOperator(
        task_id='extract_data',
        python_callable=extract_data,
        provide_context=True
    )

    transform = PythonOperator(
        task_id='transform_data',
        python_callable=transform_data,
        provide_context=True
    )

    load = PythonOperator(
        task_id='load_data',
        python_callable=load_data,
        provide_context=True
    )

    extract >> transform >> load
