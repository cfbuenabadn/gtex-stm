import pandas as pd
import os

genes = []
with open('config/genes.txt', 'r') as fh:
    for line in fh:
        genes.append(line.rstrip())
    

import json
gtex_manifest_all = json.load(open('../data/file-manifest.json'))
#gtex_manifest = json.load(open('../data/file-manifest_20.json'))
gtex_samples_all = sorted([x['file_name'].split('.')[0] for x in gtex_manifest_all])

# gtex_manifest_balanced = json.load(open('../data/file-manifest_20-balanced.json'))
# gtex_samples_balanced = sorted([x['file_name'].split('.')[0] for x in gtex_manifest_balanced])

gtex_manifest_30_1 = json.load(open('../data/file-manifest_30-1.json'))
gtex_samples_30_1 = sorted([x['file_name'].split('.')[0] for x in gtex_manifest_30_1])


#genes_50 = ['SRSF3', 'TPM3']
gtex_manifest = json.load(open('../data/file-manifest_50.json'))
gtex_samples = sorted([x['file_name'].split('.')[0] for x in gtex_manifest])

# gtex_samples = gtex_samples_50
# processed_gtex_samples = os.listdir('BamTranscriptome')
# gtex_samples_to_download = [x for x in gtex_samples_50 if (x + '.transcriptome.bed') not in processed_gtex_samples]
