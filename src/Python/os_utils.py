import os
from notebook import notebookapp
import urllib
import ipykernel
import json

def get_notebook_path():
    # Get the connection file and kernel ID
    connection_file = os.path.basename(ipykernel.get_connection_file())
    kernel_id = connection_file.split('-', 1)[1].split('.')[0]

    # Find the notebook server info
    for srv in notebookapp.list_running_servers():
        try:
            response = urllib.request.urlopen(srv['url'] + 'api/sessions')
            sessions = json.load(response)
            for sess in sessions:
                if sess['kernel']['id'] == kernel_id:
                    return sess['notebook']['path']
        except Exception:
            pass
    return None

notebook_path = get_notebook_path()
print("Notebook path:", notebook_path)