import boto3
import json


session = boto3.Session(region_name='eu-north-1')


lambda_client = session.client('lambda')

event = {
    "name": "Student"
}

response = lambda_client.invoke(
    FunctionName='HelloStudentFunction',
    InvocationType='RequestResponse',
    Payload=json.dumps(event)
)

response_payload = json.loads(response['Payload'].read())

print(f"Response from Lambda: {response_payload}")