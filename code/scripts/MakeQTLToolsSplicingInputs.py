import gzip

tissue_list = ['Brain_Anterior_cingulate_cortex_BA24',
               'Brain_Cortex',
               'Brain_Frontal_Cortex_BA9',
               'Brain_Putamen_basal_ganglia',
               'Skin_Not_Sun_Exposed_Suprapubic',
               'Liver',
               'Lung', 
               'Heart_Atrial_Appendage', 
               'Muscle_Skeletal',
               'Whole_Blood']

for tissue_id in tissue_list:

    counter = 0
    with gzip.open(f'QTLs/{tissue_id}/leafcutter.qqnorm.bed.gz', 'wb') as out_fh:
        with gzip.open(f'GTEx_Analysis_v8_sQTL_phenotype_matrices/{tissue_id}.v8.leafcutter_phenotypes.bed.gz', 'rb') as fh:
            for line in fh:
                row = line.decode().rstrip().split('\t')
                if counter == 0:
                    out_row = ['#Chr', 'start', 'end', 'pid', 'gid', 'strand'] + row[4:]
                elif counter >= 1:
                    gid = ':'.join(row[3].split(':')[3:])
                    strand = '.'
                    pid = row[3]
                    out_row = row[:3] + [pid, gid, strand] + row[4:]

                out_line = ('\t'.join(out_row) + '\n').encode()
                out_fh.write(out_line)
                    
                counter += 1