import numpy as np
import pandas as pd
import json
import argparse




def get_tissue_samples(annotation, tissue_id):
        
    all_samples = annotation.loc[(annotation.tissue_id == tissue_id) ].index.sort_values()
    
    return all_samples



def select_json(tissue_samples):
    f = open('manifests/file-manifest.json')
    data = json.load(f)

    selected_json = []
    for s in tissue_samples:
        for i in range(len(data)):
            if s == data[i]['file_name'].split('.')[0]:
                selected_json.append(data[i])

    return selected_json


def write_manifest(selected_json, file_manifest_out):

    with open(file_manifest_out, 'w') as outfile:
        json.dump(selected_json, outfile)
        
        
parser = argparse.ArgumentParser()

parser.add_argument('--tissue_samples', type=str, required=True)
parser.add_argument('--output', type=str, required=True)


if __name__ == '__main__':
    args = parser.parse_args()
    tissue_samples = args.tissue_samples.split('.')
    output = args.output
    
    manifest_json = select_json(tissue_samples)
    write_manifest(manifest_json, output)
