FROM continuumio/miniconda3
LABEL authors="giorgi"

# Install RDKit
RUN conda install -c conda-forge rdkit -y

COPY src /app

WORKDIR /app

RUN pip install fastapi uvicorn pydantic

EXPOSE 8000

CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]
