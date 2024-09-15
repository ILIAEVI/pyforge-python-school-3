FROM continuumio/miniconda3
LABEL authors="giorgi"

RUN conda install -c conda-forge rdkit -y

COPY ./requirements.txt ./requirements.txt

RUN pip install -r requirements.txt

COPY . .

EXPOSE 8000

CMD ["uvicorn", "src.main:app", "--host", "0.0.0.0", "--port", "8000"]
