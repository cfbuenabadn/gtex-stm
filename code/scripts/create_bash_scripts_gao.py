from itertools import cycle
import pandas as pd

output = [f"ebpmf_models/gao_models/bash_scripts/{str(n)}.sh" for n in range(0,2001)]
gao_genes = list(pd.read_csv('Annotations/gencode.v44.primary_assembly.protein_coding.genes.bed.gz', sep='\t').gene_id)

snmf_bashscript_pairs = zip(gao_genes, cycle(output))

snmf_bashscript_pairs = zip(gao_genes, cycle(output))
for gene, out_f in snmf_bashscript_pairs:
    n = out_f.split('/')[-1].removesuffix('.sh')
    with open(out_f, 'a') as f:
        out_dir_name = f'ebpmf_models/gao_models/RDS/{n}/'
        _ = f.write(f'[ -d "{out_dir_name}" ] || mkdir "{out_dir_name}"\n')
            
        _ = f.write(f'Rscript scripts/run_ebpmf_ad.R {gene} {n}\n')
        _ = f.write(f'touch ebpmf_models/gao_models/RDS/{n}/token.txt\n')
        #lines = ['if [ $? -ne 0 ]; then\n', f'echo "sNMF failed for gene {gene}, passing to the next command."\n', 
        #         'else\n', 'echo "sNMF succeeded for gene {gene}."\n', 'fi\n']
        #_ = f.writelines(lines)

# If there are more output files than accession numbers, the extra
# output files won't get made in the previous loop and snakemake will
# complain of missing output files. as a fail safe, let's append to
# each file in output, in effect making an empty file if a file wasn't
# made in the for loop above
for f in output:
    open(f, 'a').close()