# Python Summer School 2024

Deployment Process

This repository uses GitHub Actions to automate deployment to AWS EC2 via Docker and AWS Elastic Container Registry (ECR). The process is triggered every time changes are pushed to the deployment branch.
Prerequisites

Check: .github/workflows/deploy.yml

Before deploying, ensure the following:

    You have created and configured an AWS account.
    AWS ECR is set up to store the Docker image.
    You have an EC2 instance configured with Docker.
    GitHub Secrets are properly set for AWS_ACCESS_KEY_ID, AWS_SECRET_ACCESS_KEY, SSH_PRIVATE_KEY, and EC2_HOST (public IP of the EC2 instance).

GitHub Actions Workflow

This repository includes a deploy.yml workflow that performs the following steps:

    Trigger: The workflow is triggered when code is pushed to the deployment branch.

    AWS and Docker Environment Variables:
        AWS region, access keys, EC2 public IP, and SSH private key are pulled from GitHub Secrets.

    Checkout the Repository:
        The workflow checks out the latest code from the repository.

    Login to AWS ECR:
        Uses the aws-actions/amazon-ecr-login@v1 action to authenticate Docker with AWS Elastic Container Registry (ECR).

    Build and Push Docker Image:
        The Docker image is built from the current repository using the latest commit hash as the tag.
        The built image is pushed to AWS ECR.

    Install AWS CLI:
        The AWS CLI is installed to allow further interaction with AWS services from the workflow.

    Authenticate Docker with AWS ECR:
        Docker is authenticated using AWS ECR credentials to ensure the image can be pulled.

    Verify ECR Image:
        AWS CLI is used to verify the Docker image pushed to AWS ECR.

    Pull Docker Image:
        The Docker image is pulled from AWS ECR using the image tag based on the commit hash.

    Deploy to EC2:
        The Docker image is deployed to the EC2 instance using SSH.
        SSH connection is established with the EC2 instance, and the following tasks are performed:
            Install Docker (if not already installed).
            Start the Docker service.
            Stop and remove any running instance of the pyforge-python-school-3 container.
            Pull the latest Docker image from ECR.
            Run the new container with the specified configuration.


## Substructure search

Substructure search of chemical compounds is a crucial tool in cheminformatics, enabling researchers to identify and analyze chemical structures containing specific substructures. This method is widely applied in various fields of chemistry, including drug discovery, materials science, and environmental research. Substructure search helps scientists and engineers identify compounds with desired properties, predict reactivity, and understand the mechanisms of chemical reactions.

Modern chemical compound databases contain millions of entries, making traditional search methods inefficient and time-consuming. Substructure search utilizes algorithms that allow for quick and accurate identification of compounds with specified structural fragments. These algorithms are based on graph theory and the use of SMARTS (SMiles ARbitrary Target Specification) codes, ensuring high performance and precision in the search process.

## SMILES

A key element in the representation of chemical structures is the Simplified Molecular Input Line Entry System (SMILES). SMILES is a notation that allows a user to represent a chemical structure in a way that can be easily processed by computers. It encodes molecular structures as a series of text strings, which can then be used for various computational analyses, including substructure searches. The simplicity and efficiency of SMILES make it a widely adopted standard in cheminformatics.

Here are some examples of SMILES notation:

- Water (H₂O): O

- Methane (CH₄): C

- Ethanol (C₂H₅OH): CCO

- Benzene (C₆H₆): c1ccccc1

- Acetic acid (CH₃COOH): CC(=O)O

- Aspirin (C₉H₈O₄): CC(=O)Oc1ccccc1C(=O)O

## Examples of Substructure Search


Below are some examples of substructure searches with visual representations:

1. Searching for the Benzene Ring:
   - Substructure (Benzene): c1ccccc1
   
   <img title="a title" alt="Alt text" src="./images/c1ccccc1.png">
   
   - Example of Found Compound (Toluene): Cc1ccccc1
   
   <img title="a title" alt="Alt text" src="./images/Cc1ccccc1.png">
   
2. Searching for a Carboxylic Acid Group:
   - Substructure (Carboxylic Acid): C(=O)O
   
   <img title="a title" alt="Alt text" src="./images/C(=O)O.png">
   
   - Example of Found Compound (Acetic Acid): CC(=O)O
   
   <img title="a title" alt="Alt text" src="./images/CC(=O)O.png">

These examples illustrate how substructure searches can be used to find compounds containing specific functional groups or structural motifs. By using SMILES notation and cheminformatics tools, researchers can efficiently identify and study compounds of interest.

## Homework

Every week we will add tasks to the folder **hw**. 
Delivery details are available in the **README.md** in the **hw**

## Roadmap

As part of our homeworks, we will try to build a web service for storing and substructural search of chemical compounds

- Use the **RDKit** library to implement substructure search

<img title="a title" alt="Alt text" src="./images/1.png">

- Build RESTful API using **FastApi**

<img title="a title" alt="Alt text" src="./images/2.png">

- Containerizing our solution using **Docker**

<img title="a title" alt="Alt text" src="./images/3.png">

- Adding tests using **pytest**

<img title="a title" alt="Alt text" src="./images/4.png">

- **CI / CD**

<img title="a title" alt="Alt text" src="./images/5.png">

- We will add a **database** for storing molecules

<img title="a title" alt="Alt text" src="./images/6.png">

- **Logging**

<img title="a title" alt="Alt text" src="./images/7.png">

- We will add caching using **Redis** to optimize queries

<img title="a title" alt="Alt text" src="./images/8.png">

- **Celery** to speed up queries
<img title="a title" alt="Alt text" src="./images/9.png">
