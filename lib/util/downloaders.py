#python3
import os
import logging
import sys
from Bio import SeqIO

# We download Genome Files: gfu is Genome File Util
def DownloadGenomeToFNA(gfu, genome_ref, scratch_dir):
    """
    Inputs: GFU Object, str (A/B/C), str path
    Outputs: genome_fna_fp (str) Filepath to genome fna
    """

    GenomeToGenbankResult = gfu.genome_to_genbank({'genome_ref': genome_ref})

    logging.info("Genome File Util download result:")
    logging.info(GenomeToGenbankResult)

    genome_fna_fp = get_fa_from_scratch(scratch_dir)

    if genome_fna_fp is None:
        raise RuntimeError("GFU Genome To Genbank did not download Assembly file in expected Manner.")

    return genome_fna_fp


def get_fa_from_scratch(scratch_dir):
    """
    Careful... May not work in the Future
    Inputs:
        scratch_dir: (str) Path to work dir/ tmp etc..
    Outputs:
        FNA fp: (str) Automatic download through GenbankToGenome
    """
    
    fna_fp = None
    scratch_files = os.listdir(scratch_dir)
    all_fna_fps = []
    for f in scratch_files:
        if f[-2:] == "fa":
            all_fna_fps.append(f)

    if len(all_fna_fps) > 1:
        raise RuntimeError("Multiple .fa files in scratch directory. Expecting only one: " + \
                        ", ".join(all_fna_fps))
    elif len(all_fna_fps) == 0:
        raise RuntimeError("No .fa files in scratch directory. Program needs genome fna file to run.")
    else:
        fna_fp = all_fna_fps[0]

    return os.path.join(scratch_dir, fna_fp)


def download_genes_table(ref, dfu, op_fp):
    """
    Args:
        ref (str): ID of object within KBase, looks like 'A/B/C'
                    where A, B, and C are integers
        dfu (Object): DataFileUtil Object
        op_fp (str): Path to where to write the file
    """
    GetObjectsParams = {
            'object_refs': [ref]
            }
    
    ResultantData = dfu.get_objects(GetObjectsParams)['data'][0]['data']
    logging.info(f"Resultant data for genes table ref {ref}:")
    logging.info(ResultantData)


    if "input_genes_table" not in ResultantData:
        raise KeyError("Object must contain data key 'input_genes_table' (the file handle.)"
                        "Existing keys: " + ", ".join(ResultantData.keys()))
    else:
        KB_fh = ResultantData["input_genes_table"] 

    # Set params for shock to file
    ShockToFileParams = {
            "handle_id": KB_fh,
            "file_path": op_fp, 
            "unpack": "uncompress"
            }
    ShockToFileOutput = dfu.shock_to_file(ShockToFileParams)
    logging.info(ShockToFileOutput)
    
    logging.info(f"Downloaded genes table file for ref {ref}")


def download_model(ref, dfu, op_dir):
    """
    Args:
        ref (str): ID of object within KBase, looks like 'A/B/C'
                    where A, B, and C are integers
        dfu (Object): DataFileUtil Object
        op_dir (str): Path to dir in which we write the model file (1 or 2 lines)
    """
    GetObjectsParams = {
            'object_refs': [ref]
            }
    
    ResultantData = dfu.get_objects(GetObjectsParams)['data'][0]['data']
    logging.info(f"Resultant data for model ref {ref}:")
    logging.info(ResultantData)
    
    full_model_str = ""


    for x in ["model_string", "past_end_string", "standard_model_name"]:
        if x not in ResultantData:
            raise KeyError(f"Model should contain key {x}, but does not. Try recreating Model."
                            " Existing keys: " + ", ".join(ResultantData.keys()))
    full_model_str = ResultantData["model_string"] + "\n" + ResultantData["past_end_string"]

    logging.info("Full Model String: " + full_model_str)
    
    smn = ResultantData["standard_model_name"].replace(" ","_").replace("/","_") 
    op_fp = os.path.join(op_dir, smn + ".txt")

    with open(op_fp, 'w') as g:
        g.write(full_model_str)

    logging.info(f"Wrote model at location {op_fp}")

    return op_fp



def download_table_from_ref_to_dir(ref, ret_dp, dfu):
    """
    Args:
        ref (str): ID of object within KBase, looks like 'A/B/C'
                    where A, B, and C are integers
        ret_dp (str): Path to directory where we place
                    the downloaded tables
        dfu (Object): DataFileUtil Object
    """
    GetObjectsParams = {
            'object_refs': [ref]
            }
    
    ResultantData = dfu.get_objects(GetObjectsParams)['data'][0]['data']
    logging.info(f"Resultant data for ref {ref}:")
    logging.info(ResultantData)

    if "input_genes_table" not in ResultantData:
        raise KeyError("Object must contain data key 'input_genes_table' (the file handle.)"
                        "Existing keys: " + ", ".join(ResultantData.keys()))
    else:
        KB_fh = ResultantData["input_genes_table"] 

    # Set params for shock to file
    ShockToFileParams = {
            "handle_id": KB_fh,
            "file_path": op_fp, 
            "unpack": "uncompress"
            }
    ShockToFileOutput = dfu.shock_to_file(ShockToFileParams)
    logging.info(ShockToFileOutput)
    
    logging.info(f"Downloaded file for ref {ref}")





def DownloadFASTQs(dfu, fastq_ref_list, output_dir):
    """
    Args:
        dfu: DataFileUtil Object
        fastq_ref_list: (list<s>) list of refs 'A/B/C' A,B,C are integers
        output_dir: (s) Path to scratch directory or tmp_dir
    Returns:
        list<fastq_fp_l, orig_fq_fns>
            fastq_fp_l (list<str>): List of fastq download locations
            orig_fq_fns (list<str>): List of original fastq names
    """
    fastq_fp_l = []
    orig_fq_fns = []


    for i in range(len(fastq_ref_list)):
        crnt_fastq_ref = fastq_ref_list[i]
        logging.critical("crnt fq ref: " + crnt_fastq_ref)
        #Naming and downloading fastq/a file using DataFileUtil
        fastq_fn = "FQ_" + str(i)
        fastq_fp = os.path.join(output_dir, fastq_fn)
        get_shock_id_params = {"object_refs": [crnt_fastq_ref], "ignore_errors": False}
        get_objects_results = dfu.get_objects(get_shock_id_params)
        logging.info(get_objects_results)
        original_file_name = get_objects_results['data'][0]['info'][1]
        orig_fq_fns.append(original_file_name)
        fq_shock_id = get_objects_results['data'][0]['data']['lib']['file']['id']
        fastq_download_params = {'shock_id': fq_shock_id,'file_path': fastq_fp, 'unpack':'unpack'}
        #Here we download the fastq file itself:
        logging.info("DOWNLOADING FASTQ FILE NUMBER " + str(i+1))
        file_info = dfu.shock_to_file(fastq_download_params)
        logging.info(file_info)
        fastq_fp_l.append(fastq_fp)

    return [fastq_fp_l, orig_fq_fns]


def GetGenomeOrganismName(ws, genome_ref):
    """
    # Getting the organism name using WorkspaceClient
    ws: workspace client object
    genome_ref: (str) A/B/C
    """
    res = ws.get_objects2(
        {
            "objects": [
                {
                    "ref": genome_ref,
                    "included": ["scientific_name"],
                }
            ]
        }
    )
    scientific_name = res["data"][0]["data"]["scientific_name"]
    return scientific_name





# We want scaffold_name and description_name
# Following function is out of use right now
def get_gene_table_config_dict(genbank_fp):
    record = SeqIO.read(genbank_fp, "genbank") 
    print(record.description)
    print(record.id)
    
    gene_table_config_dict = {
            "keep_types": ["1"],
            "scaffold_name": record.id,
            "description": record.description}

    return gene_table_config_dict

def test():
    get_gene_table_config_dict(
            '/Users/omreeg/tmp/Test_Files/E_Coli_BW25113.gbk')


def download_fastq(dfu, fastq_refs_list, scratch_dir, output_fp):
    # We get multiple shock objects at once.
    get_shock_id_params = {"object_refs": fastq_refs_list, 
            "ignore_errors": False}
    get_objects_results = dfu.get_objects(get_shock_id_params)
    logging.debug(get_objects_results['data'][0])
    logging.debug(len(get_objects_results['data']))
    
    # We want to associate a ref with a filename and get a dict that has this
    # association

    raise Exception("STOP - INCOMPLETE")
    fq_shock_id = get_objects_results['data'][0]['data']['lib']['file']['id']
    fastq_download_params = {'shock_id': fq_shock_id,'file_path': fastq_fp, 'unpack':'unpack'}
    #Here we download the fastq file itself:
    logging.info("DOWNLOADING FASTQ FILE " + str(i))
    file_info = dfu.shock_to_file(fastq_download_params)
    logging.info(file_info)

