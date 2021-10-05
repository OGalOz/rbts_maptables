#python3
import os
import logging
import re

def validate_init_params(params, cfg_d):
    """
    Args:
        params: (d)
            Must contain all the parameters passed in
            as shown in function.
        cfg_d: (d)
            #models_dir: (s) Path to all models
    Returns:
        vp: (d) "Validated Params"
            genome_ref: (s)
            fastq_ref_list: (list<s>)

            (MapTnSeq vars)
            maxReads: (i)
            minQuality: (i)
            minIdentity: (i)
            minScore: (i)
            delta: (i)

    """

    #Validated params dict
    vp = {}

    if 'genome_ref' in params:
        vp['genome_ref'] = params['genome_ref']
    else:
        raise Exception("Genome Ref not passed in params.")
    if 'tnseq_model_name' in params:
        vp['tnseq_model_name'] = params['tnseq_model_name']
    else:
        raise Exception("Model Name not passed in params.")
    # We create this dict below (ref -> 1) to check that there are no repeats of fastqs
    existing_fq = {}
    if 'fastq_ref_list' in params:
        #fastq_ref will be a list since there can be multiple.
        fq_ref_list = params['fastq_ref_list']
        if len(fq_ref_list) == 0:
            raise Exception("There must be at least 1 FASTQ file. None found.")
        for ref in fq_ref_list:
            if ref in existing_fq:
                raise Exception("Do not repeat the same FASTQ reads files: Found repeat ref: " + \
                                ref)
            else:
                existing_fq[ref] = 1
        vp['fastq_ref_list'] = params['fastq_ref_list']
    else:
        raise Exception("Fastq Ref not passed in params.")


    #MapTnSeq variables
    if "maxReads" in params:
        mR = params["maxReads"]
        if mR is None or mR == "":
            # There is no limit to number of reads to go through
            # We set maxReads to 10 trillion
            vp["maxReads"] = 10**13
        elif not isinstance(mR, int):
            vp["maxReads"] = int(mR)
        else:
            vp["maxReads"] = mR
    else:
        raise Exception("maxReads not passed in params")

    if "minQuality" in params:
        mQ = params["minQuality"]
        if mQ is None or mQ == "":
            # there is no minQuality
            vp["minQuality"] = 0
        elif not isinstance(mQ, int):
            vp["minQuality"] = int(mQ)
        else:
            vp["minQuality"] = mQ
    else:
        raise Exception("minQuality not passed in params")

    if "minIdentity" in params:

        mI = params["minIdentity"]
        if mI is None or mI == "":
            # there is no minIdentity - set to default 90
            vp["minIdentity"] = 90
        elif not isinstance(mI, int):
            # It's string not int?
            vp["minIdentity"] = int(mI)
        else:
            vp["minIdentity"] = mI
    else:
        raise Exception("minIdentity not passed in params")

    if "minScore" in params:
        mS = params["minScore"]
        if mS is None or mS == "":
            vp["minScore"] = 15 
        elif not isinstance(mS, int):
            # It's string not int?
            vp["minScore"] = int(mS)
        else:
            vp["minScore"] = mS
    else:
        raise Exception("minScore not passed in params")
    if "delta" in params:
        delta = params["delta"]
        if delta is None or delta == "":
            vp["delta"] = 5
        elif not isinstance(delta, int):
            # It's string not int?
            vp["delta"] = int(delta)
        else:
            vp["delta"] = delta
    else:
        raise Exception("delta not passed in params")
    """
    if "minN" in params:
        # For Design Random Pool: minimum number of good reads for a barcode with
        # a specific mapping
        mN = params["minN"]
        if mN is None or mN == "":
            vp["minN"] = 5 
        elif not isinstance(mN, int):
            # It's string not int?
            mN = int(mN)
            vp["minN"] = mN
        else:
            vp["minN"] = mN
    else:
        raise Exception("minN not passed in params")

    if "minFrac" in params:
        # For Design Random Pool: Minimum fraction of reads w/ barcode that 
        # agree with preferred mapping
        mF = params["minFrac"]
        if mF is None or mF == "":
            vp["minFrac"] = 0.75 
        elif not isinstance(mF, float):
            # It's string not float?
            mF = float(mF)
            vp["minFrac"] = float(mF)
        else:
            vp["minFrac"] = mF
    else:
        raise Exception("minFrac not passed in params")

    if "minRatio" in params:
        # For Design Random Pool: Minimum ratio of reads for preferred mapping over
        # 2nd-most-frequent mapping (2nd most-preferred mapping)

        mR = params["minRatio"]
        if mR is None or mR == "":
            vp["minRatio"] = 8.0 
        elif not isinstance(mR, float):
            # It's string not float?
            mR = float(mR)
            vp["minRatio"] = float(mR)
        else:
            vp["minRatio"] = mR
    else:
        raise Exception("minRatio not passed in params")

    if "maxQBeg" in params:
        #  Maximum location index in query where hit to genome can occur
        mN = params["maxQBeg"]
        if mN is None or mN == "":
            vp["maxQBeg"] = 3 
        elif not isinstance(mN, int):
            # It's string not int?
            mN = int(mN)
            vp["maxQBeg"] = mN
        else:
            vp["maxQBeg"] = mN
    else:
        raise Exception("maxQBeg not passed in params")

    if "pool_description" in params:
        if params["pool_description"] == '' or params["pool_description"] is None:
            vp['pool_description'] = 'No description given.'
        else:
            vp['pool_description'] = params['pool_description']
    else:
        vp['pool_description'] = 'No description given.'

    if params['KB_Pool_Bool'] == 'yes':
        vp['KB_Pool_Bool'] = True
    else:
        vp['KB_Pool_Bool'] = False
    
    """

    if 'output_name' in params:
        if (params['output_name'] != '' and params['output_name'] is not None):
            vp['output_name'] = check_output_name(params['output_name'])
        else:
            vp['output_name'] = "Untitled"
    else:
        vp['output_name'] = "Untitled"

    return vp



# op_name is string, (output_name for app)
def check_output_name(op_name):
    op_name = op_name.replace(' ', '_')
    rgx = re.search(r'[^\w]', op_name)
    if rgx:
        logging.warning("Non-alphanumeric character in output name: " + rgx[0])
        op_name = "Default_Name_Check_Chars"
    return op_name



def validate_custom_model(custom_model_string):
    """
    Args: 
        custom_model_string (str) String of custom model. 
        Should look like the other models (2 lines, etc)
    Returns: 
        tested_model_string (str) String of custom model.
    """

    if len(custom_model_string) < 2:
        raise Exception("Custom Model form incorrect, it contains fewer than 2 "
        "characters while it should be at least 20 bp long.")
        
    if len(custom_model_string.split("\n")) > 3:
        raise Exception("Custom Model form incorrect- max amount of lines is 2.")
    
    #For now there is minimal testing. Eventually add specific tests.
    tested_model_string = custom_model_string

    return tested_model_string


