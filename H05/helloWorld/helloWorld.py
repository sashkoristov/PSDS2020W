import json

def lambda_handler(event, context):
    # TODO implement
    message = 'Hello ' + event['name']
    return {
        'message': message 
    }