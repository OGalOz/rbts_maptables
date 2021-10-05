#!python3


import os
import sys
import logging
import json
import shutil
from util.validate import validate_init_params
from util.downloaders import DownloadGenomeToFNA , DownloadFASTQs, GetGenomeOrganismName


"""
This function needs to do a number of things:
    1. Validate Parameters and make sure they're all of the right type, etc.
    2. Download genbank file and fastq file(s).
    3. Create gene table file
    4. Create config files for the function "RunFullProgram" 
"""
def PrepareProgramInputs(params, cfg_d):
    """
    Args:
        params: (d) As imported by spec file
        cfg_d: (d)
            gfu: GenomeFileUtil object
            tmp_dir: (s) path to tmp dir
            dfu: DataFileUtil object
            gt_obj: GeneTable Object 
            model_dir: (s) Path to model's directory: (Should always be scratch_dir)
            blat_cmd: (s) Path to blat loc
            unmapped_fp: (s) Path to unmapped_fp (write)
            tmpFNA_fp: (s) Path to tmp fna (write)
            trunc_fp: (s) Path to Trunc (write)
            mts_tables_dir: Path to directory (not yet created) 
            endFNA_fp: (s) Path to endFNA (write)
            MTS_cfg_fp: (s) Path to write MapTnSeq Config
            ws: Workspace Object

    Adds the following keys:
        genome_fna_fp: (File path to the Genome FNA file)

    Description:
        This function needs to do a number of things:
            1. Validate Parameters and make sure they're all of the right type, etc.
            2. Download genbank file and fastq file(s).
            3. Create gene table file
            4. Create config files for the function "RunFullProgram" 
            5. Get the TnSeq Model (transposon model)
            6. 
    """

    # validated params
    vp = validate_init_params(params, cfg_d)


    # Download genome in GBK format and convert it to fna:
    # gt stands for genome table
    genome_fna_fp = DownloadGenomeToFNA(cfg_d['gfu'], vp['genome_ref'], cfg_d['tmp_dir'])
    cfg_d['genome_fna_fp'] = genome_fna_fp
    genome_scientific_name = GetGenomeOrganismName(cfg_d['ws'], vp['genome_ref'])


    # Note that this gene table will be created at workdir/g2gt_results/genes.GC
    g2gt_results = cfg_d['gt_obj'].genome_to_genetable({'genome_ref': vp['genome_ref']})
    logging.info(g2gt_results)
    cfg_d['gene_table_fp'] = os.path.join(os.path.join(cfg_d['tmp_dir'], 'g2gt_results'), 'genes.GC')


    model_str, past_end_str = get_model_and_pastEnd_strs(vp['tnseq_model_name'])
    model_fp = os.path.join(cfg_d['model_dir'], vp['tnseq_model_name'] + '.txt')
    write_model_to_file(model_fp, model_str, past_end_str)
    cfg_d["model_fp"] = model_fp


    # FASTQs output dir
    fq_dir = os.path.join(cfg_d['tmp_dir'], "FASTQs")
    fastq_fp_l, orig_fq_fns = DownloadFASTQs(cfg_d['dfu'], vp['fastq_ref_list'], fq_dir )

    cfg_d['fastq_fp_l'] = fastq_fp_l
    cfg_d['orig_fq_fns'] = orig_fq_fns 


    MTS_cfg_d = Create_MTS_DRP_config(cfg_d, vp)

    # Creating Map Tn Seq tables directory
    os.mkdir(cfg_d["mts_tables_dir"])

    """
    # Here we write the MapTNSeq config dict 
    # and DRP config dict out to file
    with open(cfg_d["MTS_cfg_fp"], "w") as f:
        f.write(json.dumps(MTS_cfg_d, indent=4))

    logging.info("Wrote config files. MTS at {}.".format(
                cfg_d["MTS_cfg_fp"]))
    """

    return [vp, genome_scientific_name, MTS_cfg_d]


def PrepareUserOutputs(vp, cfg_d):
    """
    cfg_d:
        username: s,
        ws_id:
        ws_obj:
        workspace_name: s,
        pool_fp: s,
        dfu:
        css_style_fp: s
        Main_HTML_report_fp: s
        mts_tables_dir: s 
        model_test: b
        gene_table_fp: s

    vp:
        genome_ref: s,
        fastq_ref_list: list<s>,
       
    """

   

    # We make a directory containing the resultant files
    res_dir = os.path.join(cfg_d['tmp_dir'], "results")
    os.mkdir(res_dir)
    # We move files to this directory (mutant pool & gene table)
    #shutil.copy(cfg_d['pool_fp'], res_dir)
    shutil.copy(cfg_d["gene_table_fp"], res_dir)
    # Move the Map TnSeq Tables to the output directory
    shutil.copytree(cfg_d["mts_tables_dir"], os.path.join(res_dir, "MTS_Tables"))

    # Returning file in zipped format:-------------------------------
    file_zip_shock_id = cfg_d['dfu'].file_to_shock({'file_path': res_dir,
                                          'pack': 'zip'})['shock_id']

    dir_link = {
            'shock_id': file_zip_shock_id, 
           'name':  'results.zip', 
           'label':'rbts_mts_output_dir', 
           'description': 'The directory of outputs from running' \
            + ' Map TnSeq only, contains tables'
           }
    

    report_params = {
            'workspace_name' : cfg_d['workspace_name'],
            "message": "Finished Running MapTnSeq"
            }

    report_params["file_links"] = [dir_link]


    return report_params


# MTS - Map Tn Seq, DRP - Design Random Pool
def Create_MTS_DRP_config(cfg_d, vp):
    """
    Args:
        cfg_d: (as above in PrepareProgramInputs)
            Plus:
            fastq_fp_l: (list<s>) List of file paths
            genome_fna_fp: (File path to the Genome FNA file)

        vp: (d) must contain all used cases below
    Returns:
        MTS_cfg
    """
    # Here we create the config dicts
    map_tnseq_config_dict = {

        "values": {
            "debug": False,
            "keepblat8": True,
            "keepTMPfna": True,
            "flanking": 5,
            "wobbleAllowed": 2,
            "tmp_dir": cfg_d["tmp_dir"],
            "tileSize": 11,
            "stepSize": 11,
            "blatcmd": cfg_d["blat_cmd"],
            "model_fp": cfg_d["model_fp"],
            "maxReads": vp["maxReads"],
            "minQuality": vp["minQuality"],
            "minIdentity": vp["minIdentity"],
            "minScore": vp["minScore"],
            "delta": vp["delta"],
            "fastq_fp_list":  cfg_d['fastq_fp_l'],
            "orig_fq_fns":  cfg_d['orig_fq_fns'],
            "genome_fp": cfg_d['genome_fna_fp']
        }
    }


    return map_tnseq_config_dict 





def get_model_and_pastEnd_strs(standard_model_name):
    """
    Args:
        standard_model_name (str): Name of the model. Should be one of below.
    Description:
        In this function we get the two parts of the model-
        The model string, which is the part of the transposon in which the barcode sits.
        And the past end string, which is after the transposon.
    Returns:
        model_str, past_end_str 
        
    """

    if standard_model_name ==  "Sc_Tn5":
        model_str = "nnnnnGATGTCCACGAGGTCTCTNNNNNNNNNNNNNNNNNNNNCGTACGCTGCAGGTCGACCAGCAGCTATGACATGAAGATGTGTATAAGAGACAG" 
        past_end_str = "GGAAGGGCCCGACGTCGCATGCTCCCGGCCGCCATGGCGGCCGCGGGAATTCGATTGGGCCCAGGTACCAACTACGTCAGGTGGCACTTT"
    elif standard_model_name == "ezTn5_Tet_Bifido":
        model_str = "nnnnnnGATGTCCACGAGGTCTCTNNNNNNNNNNNNNNNNNNNNCGTACGCTGCACCTCGACAGATGTGTATAAGAGACAG" 
        past_end_str = ""
    elif standard_model_name == "ezTn5_kan1":
        model_str = "nnnnnnCTAAGATGTCCACGAGGTCTCTNNNNNNNNNNNNNNNNNNNNCGTACGCTGCAGGTCGACAGATGTGTATAAGAGACAG"
        past_end_str = ""
    elif standard_model_name == "ezTn5_kanU":
        model_str = "nnnnnnGATGTCCACGAGGTCTCTNNNNNNNNNNNNNNNNNNNNCGTACGCTGCAGGTCGACAGATGTGTATAAGAGACAG"
        past_end_str = ""
    elif standard_model_name == "magic_Tn5":
        model_str = 'nnnnnnGATGTCCACGAGGTCTCTNNNNNNNNNNNNNNNNNNNNCGTACGCTGCAGGTCGACCAGCGGCCGGCCGGTTGAGATGTGTATAAGAGACAG'
        past_end_str = 'TCGACGGCTTGGTTTCATAAGCCATCCGCTTGCCCTCATCTGTTACGCCGGCGGTAGCCGGCCAGCCTCGCAGAGCAGGATTCCCGTTGA'
    elif standard_model_name == "magic_mariner":
        model_str = 'nnnnnnGATGTCCACGAGGTCTCTNNNNNNNNNNNNNNNNNNNNCGTACGCTGCAGGTCGACCAGCGGCCGGCCAGACCGGGGACTTATCAGCCAACCTGT' 
        past_end_str = 'TATGTGTTGGGTAACGCCAGGGTTTTCCCAGTCACGACGTTGTAAAACGACGGCCAGTGAATTAATTCTTGCTTATCGGCCAGCCTCGCAGAGCAGGATTCCCGTTGAGCACCGCCAGGTGCGAATAAGGGACAGTGAAGAAG'
    elif standard_model_name == "magic_mariner.2":
        model_str = 'nnnnnnnnnnnnnGATGTCCACGAGGTCTCTNNNNNNNNNNNNNNNNNNNNCGTACGCTGCAGGTCGACCAGCGGCCGGCCAGACCGGGGACTTATCAGCCAACCTGT' 
        past_end_str = 'TATGTGTTGGGTAACGCCAGGGTTTTCCCAGTCACGACGTTGTAAAACGACGGCCAGTGAATTAATTCTTGCTTATCGGCCAGCCTCGCAGAGCAGGATTCCCGTTGAGCACCGCCAGGTGCGAATAAGGGACAGTGAAGAAG'
    elif standard_model_name == "pHIMAR_kan":
        model_str = 'nnnnnnCGCCCTGCAGGGATGTCCACGAGGTCTCTNNNNNNNNNNNNNNNNNNNNCGTACGCTGCAGGTCGACGGCCGGCCAGACCGGGGACTTATCAGCCAACCTGT' 
        past_end_str = 'TATGTGTTGGGTAACGCCAGGGTTTTCCCAGTCACGACGTTGTAAAACGACGGCCAGTGAATTAATTCTTGAAGA' 
    elif standard_model_name == "pKMW3":
        model_str = 'nnnnnnCGCCCTGCAGGGATGTCCACGAGGTCTCTNNNNNNNNNNNNNNNNNNNNCGTACGCTGCAGGTCGACGGCCGGCCAGACCGGGGACTTATCAGCCAACCTGT' 
        past_end_str = 'TATGTGTTGGGTAACGCCAGGGTTTTCCCAGTCACGACGTTGTAAAACGACGGCCAGTGAATTAATTCTTGAAGA'
    elif standard_model_name == "pKMW3_universal":
        model_str = 'nnnnnnGATGTCCACGAGGTCTCTNNNNNNNNNNNNNNNNNNNNCGTACGCTGCAGGTCGACGGCCGGCCAGACCGGGGACTTATCAGCCAACCTGT' 
        past_end_str = 'TATGTGTTGGGTAACGCCAGGGTTTTCCCAGTCACGACGTTGTAAAACGACGGCCAGTGAATTAATTCTTGAAGA'
    elif standard_model_name == "pKMW7":
        model_str = 'nnnnnnCGCCCTGCAGGGATGTCCACGAGGTCTCTNNNNNNNNNNNNNNNNNNNNCGTACGCTGCAGGTCGACGGCCGGCCGGTTGAGATGTGTATAAGAGACAG' 
        past_end_str = 'TCGACGGCTTGGTTTCATCAGCCATCCGCTTGCCCTCATCTGTTACGCCGGCGGTAGCCGGCCAGCCTCGCAGAGC'
    elif standard_model_name == "pKMW7_U":
        model_str = 'nnnnnnGATGTCCACGAGGTCTCTNNNNNNNNNNNNNNNNNNNNCGTACGCTGCAGGTCGACGGCCGGCCGGTTGAGATGTGTATAAGAGACAG'
        past_end_str = 'TCGACGGCTTGGTTTCATCAGCCATCCGCTTGCCCTCATCTGTTACGCCGGCGGTAGCCGGCCAGCCTCGCAGAGC'
    else:
        raise Exception(f"Model name {standard_model_name} not recognized.")


    logging.info(f"Model String: '{model_str}'."
                f" Past End String: '{past_end_str}'.")

    return model_str, past_end_str


def write_model_to_file(op_fp, model_str, past_end_str):
    """
    Args:
        op_fp (str): Path to output model
        model_str (str): Complete model string
        past_end_str (str): Complete past End str (short sequence of vector after transposon)
    Description:
        We take the two components of the model and write out to a file.
    """

    with open(op_fp, "w") as g:
        g.write(model_str + "\n" + past_end_str)
    
    logging.info(f"Wrote TnSeq model to file {op_fp}.")


def clear_dir(dir_path):
    for filename in os.listdir(dir_path):
        file_path = os.path.join(dir_path, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            print('Failed to delete %s. Reason: %s' % (file_path, e))


def test():

    return None

def main():

    test()
    return None

if __name__ == "__main__":
    main()
