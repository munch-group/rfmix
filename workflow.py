
from gwf import Workflow, AnonymousTarget
import os
from pathlib import Path
import pandas as pd
from collections import defaultdict

# conda env create -f binder/environment.yml
conda_env = 'rfmix-workflow'

def prep_rfmix(analysis, vcf_files, ref_samples, query_samples, out_suffix, output_dir):

    inputs = vcf_files

    output_ref = output_dir + analysis + "/" + out_suffix+"_ref.bcf"
    output_query = output_dir + analysis + "/" + out_suffix+"_query.bcf"

    outputs = {'reference': output_ref, 'query': output_query}
    options = {
        "cores": 2,
        "memory": "30g",
        "walltime": "4:00:00",
        "account": "baboondiversity"
    }
    spec = f"""
    eval "$(conda shell.bash activate {conda_env})"

    bcftools concat {' '.join(vcf_files)} | bcftools view  -s {",".join(ref_samples)}  -q 0.01:minor \
    -O b -o {output_ref} --force-samples

    bcftools index {output_ref}

    bcftools concat {' '.join(vcf_files)} | bcftools view -s {",".join(query_samples)} -q 0.01:minor \
    -O b -o {output_query} --force-samples

    bcftools index {output_query}
    """    
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def rfmix(chrom, query, reference, sample_map, genetic_map, output_path, e=3, G=100, reanalyze_ref=True):
    # Run with "low" memory, then with higher memory to complete all jobs
    output = output_path + chrom
    genetic_map
    inputs = [query, reference, sample_map, genetic_map]
    outputs = [output + ".msp.tsv"]
    options = {
        "cores": 10,
        "memory": "200g",
        "walltime": "12:00:00",
        "account": "baboondiversity"
    }
    spec = f"""
    eval "$(conda shell.bash activate {conda_env})"
    rfmix -f {query} -r {reference} -m {sample_map} -g {genetic_map} -o {output} \
        --chromosome={chrom} -e {e} -G {G} {'--reanalyze-reference' if reanalyze_ref else ''}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# def rfmix_workflow(working_dir=os.getcwd(), 
def rfmix_workflow(gwf=Workflow(working_dir=os.getcwd()), 
                   analyzes=None, output_dir=None, 
                   vcf_files=None,
                   autosome_rec_map=None,
                   x_rec_map=None):

    # gwf = Workflow(working_dir=working_dir)

    targets = defaultdict(list)

    # make a big VCF with all autosomes
    autosomes = [f'chr{x}' for x in range(1, 21)]
    target_list = gwf.map(prep_rfmix, analyzes, name="autosomes",
            extra={"vcf_files": [vcf_files[x] for x in autosomes],
                   "output_dir": output_dir,
                    "out_suffix": "aut"})
    targets['prep_auto'].extend(target_list)

    # make VCF with female X chromosomes
    target_list = gwf.map(prep_rfmix, analyzes, name="chrX_female",
            extra={"vcf_files": [vcf_files['chrX']],
                   "output_dir": output_dir,
                   "out_suffix": "X_female"
                })
    targets['prep_x_female'].extend(target_list)

    # make VCF with all X chromosomes
    target_list = gwf.map(prep_rfmix, analyzes, name="chrX_all",
            extra={"vcf_files": [vcf_files['chrX_all']],
                   "output_dir": output_dir, 
                   "out_suffix": "X_all"
                })
    targets['prep_x_all'].extend(target_list)


    for i in range(len(analyzes)):
        n = analyzes[i]["analysis"]
        file_name = output_dir + n

        gwf.map(rfmix, autosomes, name="rfmix_"+n,
                    extra={"query": file_name+"/aut_query.bcf", 
                           "reference": file_name+"/aut_ref.bcf",
                        "sample_map": file_name+"/ref_names.txt",
                        "genetic_map": autosome_rec_map,
                        "output_path": file_name+"/"})
        targets['rfmix_auto'].extend(target_list)

        target_list = gwf.map(rfmix, ["chrX"], name="rfmix_X_female_"+n,
                    extra={"query": file_name+"/X_female_query.bcf", 
                           "reference": file_name+"/X_female_ref.bcf",
                        "sample_map": file_name+"/ref_names.txt",
                        "genetic_map": x_rec_map,
                        "output_path": file_name+"/female_"})
        targets['rfmix_x_female'].extend(target_list)        

        target_list = gwf.map(rfmix, ["chrX"], name="rfmix_X_all_"+n,
                    extra={"query": file_name+"/X_all_query.bcf", 
                           "reference": file_name+"/X_all_ref.bcf",
                        "sample_map": file_name+"/ref_names.txt",
                        "genetic_map": x_rec_map,
                        "output_path": file_name+"/all_"})
        targets['rfmix_x_all'].extend(target_list)        
        

    return gwf, targets


# we need to assign the workflow to the gwf variable to allow the workflow to be
# run separetely with 'gwf run' in the submoduleB dir:

# gwf, A_targets  = rfmix_workflow(
#                           working_dir=working_dir,
#                         #   working_dir='./relate',
#                           input_files=['./input.txt'],
#                           output_dir='./A_outputs',
#                           analyzes=analyzes, 
#                           output_dir=output_dir, 
#                           ref_name_list=ref_name_list,
#                           path_to_vcfs=path_to_vcfs, x_all_path=x_all_path
#                           )