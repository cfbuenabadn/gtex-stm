import pandas as pd
import os

genes = []
with open('config/genes.txt', 'r') as fh:
    for line in fh:
        genes.append(line.rstrip())
    

import json
gtex_manifest = json.load(open('../data/file-manifest.json'))
gtex_samples = sorted([x['file_name'].split('.')[0] for x in gtex_manifest])

