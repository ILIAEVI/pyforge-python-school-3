FROM continuumio/miniconda3
LABEL authors="giorgi"

# Install RDKit
RUN conda install -c conda-forge rdkit -y

COPY src /app

WORKDIR /app

COPY ./requirements.txt ./requirements.txt

RUN pip install -r requirements.txt

EXPOSE 8000

CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]
