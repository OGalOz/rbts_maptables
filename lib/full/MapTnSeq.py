#!python3

import os
import sys
import json
import logging
import argparse
import re
import subprocess
import copy
from collections import namedtuple
from typing import Optional, Tuple

"""
Notes:
    Program called through "RunMapTnSeq"
    Overall scratch dir will be scratch/MapTStmp
    temporary fna file will be called "TMP.fna"
    trunc file will be called "TRUNC.fastq"
    PastEndFNA will be called "PastEnd.fna"
    FASTQ file will be unzipped.
"""


def RunMapTnSeq(input_args, DEBUGPRINT=False):
    """
    Args:
        input_args: (dict) All necessary 
            keys will be listed under the function "ValidateInputs".
        DEBUGPRINT: (bool) Used to debug and decipher the program
            by printing variables to json files in the tmp dir

    This function runs the entire process of MapTnSeq
        The variable "output_fp" in input_args is where 
        the output of the program is written to.
        The output file is tab-delimited and contains a list
        of the usable reads, where, by column, the information
        is presented: read name, barcode, which scaffold, pos (base-pair),
        strand, unique or not flag, beginning and end of hit to genome in 
        read after trimming the transposon seq, the bit score, and
        the % identity (10 overall fields)
    """
    
    # We check that inputs are correct and add "model" and "pastEnd" (both str):
    parsed_vars = ValidateInputs(input_args)
    parsed_vars["DEBUGPRINT"] = DEBUGPRINT

    # parses input fastq file and finds barcodes, generates the dict "mapnames", 
    # writes TMPFNA and TRUNC
    ret_d = find_barcodes_and_end_of_transposon(parsed_vars)


    # It is possible the model wasn't found at all, and then ret_d = -1
    if isinstance(ret_d, int) and ret_d == -1:
        logging.warning("Stopping current run on:\n"
        "Model {}\n".format(parsed_vars['model_fp']) + \
        "FASTQ {}\n".format(parsed_vars['fastq_fp']) + \
        "Genome {}".format(parsed_vars['genome_fp']))
        return [-1,0]
    else:
        parsed_vars.update(ret_d)
        if parsed_vars['modeltest']:
            # If we are only testing for correctness of the model
            # we exit the program at this phase
            logging.info("Model: {} passed.".format(parsed_vars['model_fp']))
            return [0, ret_d['nTryToMap']]


    if parsed_vars['pastEnd'] is not None and parsed_vars['pastEnd'] != '':
        # generates the dict "hitsPastEnd", runs BLAT8 on PastEnd
        hitsPastEnd = pastEndBLAT8(parsed_vars)
    else:
        hitsPastEnd = {}

    parsed_vars['hitsPastEnd'] = hitsPastEnd
  

    # Creates Output File Handle, to be closed later.
    # Also writes much to MapTnSeq table output
    # Runs BLAT8 on the genome
    OUTPUT, HG_d = Map_to_genomeBLAT8(parsed_vars)
    parsed_vars['OUTPUT'] = OUTPUT
    parsed_vars.update(HG_d)


    # Writes the file "unmapped_fp" in parsed_vars
    # also, unlinks the file tmpFNA_fp if debug isn't True
    nUnmapped = WriteUnMappedFile(parsed_vars)
    parsed_vars['nUnmapped'] = nUnmapped

    # We add to "OUTPUT" file handle all the pastEnd values
    PrintOutHitsPastEnd(parsed_vars)

    parsed_vars['OUTPUT'].close()

    # Generate text report (returns string)
    text_report, text_list, report_dict = CreateTextReport(
                    parsed_vars['nLong'], 
                    parsed_vars['nReads'],
                    parsed_vars 
                    )

    logging.info(text_report)


    return_dict = {
            "fastq_fp": parsed_vars["fastq_fp"],
            "text_report_dict": report_dict
            }

    return return_dict



def ValidateInputs(input_args):
    """
    input_args should be a dict with the following:
        debug: (bool)
        keepblat8: (bool) Whether or not to keep the blat8 files
        keepTMPfna: (bool) Whether or not to keep the TMP fna file
        modeltest: (bool)
        maxReads: (int or None)
        minQuality: (int)
        flanking: (int)
        wobbleAllowed: (int)
        minIdentity: (float)
        minScore: (float)
        delta: (int)
        tileSize: (int)
        stepSize: (int)
        tmp_dir: (str) Path to working directory
        blatcmd: (str) BLAT command
        unmapped_fp: (str) Path to unmapped file to write to
        tmpFNA_fp: (str) Path to Temp FNA file to write to
        endFNA_fp: (str) Path to a pastEnd fna file to write to
        trunc_fp: (str) Path to write TRUNC file to write to
        genome_fp: (str) Path to genome FASTA file to read from
        model_fp: (str) filepath to model to read from
        fastq_fp: (str) Path to fastq file to read from
                Note, fastq file should be unzipped. If it contains 
                paired-end reads, the second read will be ignored.
        output_fp: (str) Path to output file to write to
    """

    logging.debug(input_args.keys())

    # Checking debug and modeltest
    for v in ["debug", "modeltest", "keepblat8", "keepTMPfna"]:
        if v not in input_args or not isinstance(input_args[v], bool):
            raise Exception(v + " must be included in inputs and must be bool.")

    # Checking ints
    for v in ["maxReads", "minQuality", "flanking",
            "wobbleAllowed", "delta",
            "tileSize", "stepSize",
            "minIdentity", "minScore"]:
        if v not in input_args or not isinstance(input_args[v], int):
            if v == "maxReads" and input_args[v] is None:
                continue
            else:
                raise Exception(v + " must be in inputs and must be int.")

    # Checking strs
    for v in ["tmp_dir", "unmapped_fp", "trunc_fp", "genome_fp", "endFNA_fp",
            "model_fp", "fastq_fp", "output_fp", "blatcmd", "tmpFNA_fp"]:
        if v not in input_args or not isinstance(input_args[v], str):
            raise Exception(v + " must be in inputs and must be str.")


    # Parsing Model File into "model" and "pastEnd":
    # Model
    model, pastEnd = ParseModel(input_args['model_fp'])

    input_args['model'] = model
    input_args['pastEnd'] = pastEnd

    input_args['nMapped'] = 0 
    input_args['nMapUnique'] = 0 
    input_args['nPastEndIgnored'] = 0     
    input_args['nPastEndTrumps'] = 0  
    input_args['nameToBarcode'] = {}

    return input_args


def ParseModel(model_fp):
    """
    Args:
        model_fp: (str) is path to model file, should be 1 or 2 lines
    Returns:
        model: (str) the string of the first line of the model file.
        pastEnd: (str) string of the second line of the model file. Could be ''
    """
    with open(model_fp, "r") as g:
        model = g.readline()
        model = model.rstrip()
        if not re.match(r'^n*[ACGT]+N+[ACGT]+$', model):
            raise Exception("Invalid model: " + model + "\n for file " + \
                    input_args['model_fp'])
        pastEnd = g.readline()
        pastEnd = pastEnd.rstrip()
        if pastEnd != '':
            if not re.match(r'^[ACGT]+$', pastEnd):
                raise Exception("Invalid past-end sequence: " + pastEnd + \
                        "\n for file " + input_args['model_fp'])

    return [model, pastEnd]




def CreateTextReport(nLong, nReads, inp_dict):
    """
    We create a string with info about program:

    nReads (int)
    nLong (int)
    inp_dict: (dict) Contains:
        fastq_fp: (str)
        nameToBarcode (dict)
        nTryToMap (int) 
        nMapped (int)
        nMapUnique (int)
        hitsPastEnd (dict)
        nPastEndIgnored (int)
        nPastEndTrumps (int)
        minReadLength (int)
        minScore (int)
    """

    # text report dict - for HTML table later
    trd = {
            "fastq_fn": inp_dict["fastq_fp"].split("/")[-1],
            "R_proc": nReads,
            "long_enough": nLong,
            "minRlen": inp_dict['minReadLength'],
            "tot_BC": len(inp_dict['nameToBarcode'].keys()),
            "nTryToMap": inp_dict['nTryToMap'],
            "nMapped": inp_dict['nMapped'],
            "Uniquely_Mapped": inp_dict['nMapUnique'],
            "total_hits_pastEnd": len(inp_dict['hitsPastEnd'].keys()),
            "Hits_pE": len(inp_dict['hitsPastEnd'].keys()) - inp_dict['nPastEndIgnored'],
            "nPastEndIgnored": inp_dict["nPastEndIgnored"],
            "nPastEndTrumps": inp_dict['nPastEndTrumps']
            }

    text_str = "A good run has a Mapped % that's greater than about 10%\n"
    text_str += "    Mapped %: {:3f}\n".format(inp_dict['nMapped']/nReads)
    text_str += "Counts:\n "
    text_str += "    Reads processed (nReads): " + str(nReads) + "\n" \
              + "    Long-enough (nLong): " + str(nLong) + "\n" \
              + "        (Long-enough means Each read must have at least length " \
              + str(inp_dict['minReadLength']) + ")\n" \
              + "    Barcodes found : " + str(len(inp_dict['nameToBarcode'].keys())) + "\n" \
              + "    Mapping attempted for (nTryToMap): " + str(inp_dict['nTryToMap']) + "\n" \
              + "        (nTryToMap counts the reads where the transposon ends before the end of the read\n" \
              + "        and there are at least {} nucleotides between).".format(inp_dict['minScore']) + "\n" \
              + "    Mapped (nMapped): " + str(inp_dict['nMapped']) + "\n" \
              + "        nMapped is the total number of good enough reads written to the output." + "\n" \
              + "    Uniquely (Only one good loc) (nMapUnique): " + str(inp_dict['nMapUnique']) + "\n" \
              + "    Hits past end of transposon (hitsPastEnd - nPastEndIgnored): " \
              +  str(len(inp_dict['hitsPastEnd'].keys()) \
                        - inp_dict['nPastEndIgnored']) + "\n" \
              + "    Weak/ambiguous (nPastEndIgnored): " + str(inp_dict['nPastEndIgnored']) + "\n" \
              + "    Trumped hit to genome (nPastEndTrumps):" + str(inp_dict['nPastEndTrumps']) \
              + " times\n"
    if nReads > 0:
        text_str += "Proportions: \n" \
                  + "    Long-enough % (nLong/nReads): {:.3f}\n".format(nLong/nReads) \
                  + "    Barcode % (barcodesFound/nReads): {:3f}\n".format(
                          len(inp_dict['nameToBarcode'].keys())/nReads) \
                  + "    Attempted % (nTryToMap/nReads): {:3f}\n".format(inp_dict['nTryToMap']/nReads) \
                  + "    Mapped % (nMapped/nReads): {:3f}\n".format(inp_dict['nMapped']/nReads) \
                  + "    Past-end % (hits past end/ nReads): {:3f}\n".format(
                    (len(inp_dict['hitsPastEnd'].keys())-inp_dict['nPastEndIgnored'])/nReads)

    text_list = text_str.split("\n") 



    return [text_str, text_list, trd]









def PrintOutHitsPastEnd(inp_dict):
    """
    This function writes "hits past end" to OUTPUT file handle

    Args:
        hitsPastEnd (dict)
        nameToBarcode (dict)
        OUTPUT: file handle to output

    """
    for k in inp_dict['hitsPastEnd']:
        read, score = k, inp_dict['hitsPastEnd'][k]
        # read name (short name), barcode, scaffold
        if score > 0.0:
            inp_dict['OUTPUT'].write("\t".join([
                read, 
                str(inp_dict['nameToBarcode'][read]),
                "pastEnd",
                "",
                "",
                "",
                "",
                "",
                "",
                "\n"
                ]))






def WriteUnMappedFile(inp_dict):
    """
    Args:

    inp_dict: (dict) must contain:
        unmapped_fp: (str)
        tmpFNA_fp: (str)
        mapnames: (dict)
        nameToBarcode: (dict)
        debug: (bool)
        keepTMPfna: (bool) Whether or not to unlink tmpFNA file

    vars created:
        nUnmapped (int) number of lines in unMapped file

    Writing to unmapped file at "unmapped_fp", we check each barcode in 
        tmpFNA_fp and see if it's in the dict "mapnames"
    """

    nUnmapped = 0
    tmpFNA_fp = inp_dict['tmpFNA_fp']
    UNMAPPED = open(inp_dict['unmapped_fp'], "w")
    FNA = open(tmpFNA_fp, "r")
    line_num = 1
    header = FNA.readline()
    while header != "":
        header = header.rstrip()
        if not header[0] == ">":
            raise Exception("Cannot parse {} at line {} in {}".format(header, 
                line_num, tmpFNA_fp))
        name = header[1:]
        seq = FNA.readline()
        line_num += 1
        seq = seq.rstrip()
        if not re.match(r'[ACGTN]+$', seq):
            raise Exception("Cannot parse {} in {} at line {}".format(
                            seq,tmpFNA_fp, line_num))
        if name in inp_dict['mapnames'] and inp_dict['mapnames'][name] == 0:
            UNMAPPED.write('>{} {}\n{}\n'.format(
                inp_dict['nameToBarcode'][name],
                name,
                seq
                ))
            nUnmapped += 1
        header = FNA.readline()
        line_num += 1
    

    FNA.close()
    UNMAPPED.close()
    logging.info("Wrote {} unmapped reads to {} in fasta format\n".format(
        nUnmapped,
        inp_dict['unmapped_fp']))

    if not inp_dict['keepTMPfna']:
        os.unlink(tmpFNA_fp)

    return nUnmapped







def Map_to_genomeBLAT8(inp_dict):
    """
    Description:
        This function first runs BLAT8 on the Good Reads file (TMPFNA) as a query
        and the genome FNA as a database.
        If no pastEnd in the model, hitsPastEnd is an empty dict.
        ot takes the results of BLAT8; for each line in the BLAT8 file, it either
        adds it to a list to be analyzed through HandleGenomeBlat later, or analyzes
        the current list with HandleGenomeBlat.

    inp_dict:
        debug: (bool)
        unmapped_fp: (str)
        hitsPastEnd: (dict) (shortname to blat8 score of query)
        mapnames: (dict)
        tmpFNA_fp: (str)
        genome_fp: (str) One of the inputs to the program (as genome) should
            be a FA file
        blatcmd: (str) (Location of blat8 in OS)
        tmp_dir: (str) path to working dir
        blatcmd: (str) Path to BLAT command, i.e. Bin/blat
        minScore: (int) minimum score for mapping to genome or past-end
        minIdentity: (float) minimum %identity for mapping to genome or past-end
        tileSize: (int) size of an alignment tile
        stepSize: (int) distance between the starting bases of alignment tiles
            (will overlap if stepSize<tileSize)
        output_fp: (str) path to output_file

        The below are used for Handle Genome BLAT

        nPastEndTrumps: (int)  hit to past-end (close to) as good as hit to genome
        nPastEndIgnored: (int) weak hit to past-end ignored 
        nMapUnique: (int)
        nMapped: (int)
        nameToBarcode: (dict)
        delta: (int) minimum difference in score for considering a hit unique.
           
    Returns:
        OUTPUT file handle for output file

    """
    b8_d = inp_dict
    b8_d['queriesFile'] = inp_dict['tmpFNA_fp']
    b8_d['dbFile'] = inp_dict['genome_fp']
    blat8_fp = os.path.join(inp_dict['tmp_dir'], 'blat8_genome')
    b8_d['blat8_fp'] = blat8_fp
    RunBLAT8(b8_d)

    logging.info("Parsing " + blat8_fp)

    b = open(blat8_fp, "r")
    # OUTPUT is file handle for output
    OUTPUT = open(inp_dict['output_fp'], "w")

    # Handle Genome BLAT Dict
    HG_d = {
        "nPastEndTrumps": inp_dict[ "nPastEndTrumps"],
        "nPastEndIgnored": inp_dict[ "nPastEndIgnored"],
        "nMapUnique": inp_dict[ "nMapUnique"],
        "nMapped": inp_dict[ "nMapped"],
        "nameToBarcode": inp_dict[ "nameToBarcode"],
        "delta": inp_dict[ "delta"],
        "OUTPUT": OUTPUT
        }


    lines = []
    c_line = b.readline()

    while c_line != "":
        # Removing newline \n
        c_line = c_line.rstrip()
        F = c_line.split('\t')

        # query name, subject name, identity, length, mismatches, number gaps, 
        # query Beginning, query End, subject Beginning loc, subject End loc, e-value, bit score
        query, subject, identity, l, mm, gs, qB, qE, sB, sE, evl, score = F

        inp_dict['mapnames'][query] = 1

        # Sometimes there is the same query twice in a row 
        if len(lines) == 0 or query == lines[0][0]:
            lines.append(F)
        else:
            # We have a new query, so existing query will be "Handled", and new 
            # query will be added to "lines"
            HandleGenomeBLAT(lines, inp_dict['hitsPastEnd'], HG_d, 
                    inp_dict['debug'])
            # We make lines contain the new query
            lines = [F]

        c_line = b.readline()

    # We run it a last time on the remaining lines
    HandleGenomeBLAT(lines, inp_dict['hitsPastEnd'], HG_d, inp_dict['debug'])

    b.close()

    if not inp_dict['keepblat8']:
        os.unlink(blat8_fp)

    # We return file handle, and HG_d (Handle Genome Dict)
    return [OUTPUT, HG_d]
           



# This only runs if there is a "pastEnd" part to the model (a second line)
def pastEndBLAT8(inp_dict):
    """
    We run BLAT8 to determine the 'bitscore' of each sequence in TMPFNA 
        (good sequences) against the pastEnd and keep the highest score 
        of each readname to compare later against hits into the genome.
    We keep the results in a dict called "hitsPastEnd"


    Args:

      inp_dict: (dict) must contain
          endFNA_fp: (str) filepath to End FNA
          tmpFNA_fp: (str) filepath to temp FNA
          pastEnd: (str) an optional part of the end of the model
              that contains the sequence "past the end" of the 
              transposon that might arise from residual intact plasmid.
          blatcmd: (str) BLAT command
          tmp_dir: (str) path to working directory
          minScore: (int) minimum score for mapping to genome or past-end
          minIdentity: (int) minimum %identity for mappinf to genome or past-end
          tileSize: (int) size of an alignment tile
          stepSize: (int) distance between the starting bases of alignment tiles
              (will overlap if stepSize<tileSize)
          debug: (bool)
      Optional:
          unmapped
          mapnames: (dict)

    Returns:
        hitsPastEnd: (dict) read (str) to score (float) of match to past-end sequence
            of a specific query

    """

    # read to score of match to past-end sequence
    hitsPastEnd = {}

    with open(inp_dict['endFNA_fp'], "w") as g:
        g.write(">pastend\n{}\n".format(inp_dict['pastEnd']))

    # Making the input to RunBLAT8:
    blat8_inp_d = inp_dict
    # QueriesFile is all the succesful sequences with pieces in genome
    blat8_inp_d['queriesFile'] = inp_dict['tmpFNA_fp']
    # dbFile will be just a single sequence: the pastEnd
    blat8_inp_d['dbFile'] = inp_dict['endFNA_fp']
    pastEnd_blat8_fp = os.path.join(inp_dict['tmp_dir'], 'blat8_pastEnd')
    blat8_inp_d['blat8_fp'] = pastEnd_blat8_fp 

    logging.info("Running PastEnd BLAT8 test")
    RunBLAT8(blat8_inp_d)

    logging.info("Parsing past-end hits to {}\n".format(pastEnd_blat8_fp))

    # We open the blat8 file we just created
    B = open(pastEnd_blat8_fp, "r")

    crnt_line =  B.readline()
    while crnt_line != '':
        F = crnt_line.rstrip().split('\t')

        # query, subject, identity, length, mismatches, number gaps, query Beginning,
        #   query End, subject Beginning, subject End, e-value, bit score
        query, subject, identity, l, mm, gs, qB, qE, sB, sE, e_v, score = F
        # We keep the highest pastEnd hit score for each query - which is a shortname
        if query not in hitsPastEnd: 
            hitsPastEnd[query] = float(score)
        elif hitsPastEnd[query] < float(score):
            hitsPastEnd[query] = float(score)

        # This part can be reformatted
        if 'unmapped_fp' in inp_dict: 
            if 'mapnames' in inp_dict:
                inp_dict['mapnames'][query] = 1

        crnt_line = B.readline()

    B.close()

    #In case of debugging or some other reason, we can keep the blat8 file for checking later.
    if not inp_dict['keepblat8']:
        os.unlink(pastEnd_blat8_fp)

    return hitsPastEnd 



# This function finds barcodes and end of transposon and writes remaining 
# sequence to TMP FNA
# FASTQ must already be unzipped
def find_barcodes_and_end_of_transposon(inp_dict):
    """
    This function does the following:
        Reads either the entire FASTQ input file or up to "MaxReads"*4 lines (*4 accounts for the
        4 lines per read).
        First it makes sure the FASTQ sequence is properly formatted - every sequence
            has at least a single nucleotide. Additionally
            the sequence and quality must be the same length, and 
            the sequence must be longer than the model and the 'minScore'
            variable.
        Then it takes the sequence and the model and looks for the "flanking"
        parts of the model in the sequence, so for example: 5 nucleotides
        to the left and right of the barcode in the model in the sequence.
        It then checks the quality of the barcode reading and compares each nt
        to the variable "minQuality". If the barcode length and quality
        are good, then it keeps the barcode and its start location.
        Then it maps the "short" name (name of seq before space) to the barcode
        in the dict "nameToBarcode".
        Then it looks for the end of the model within the sequence in the function
        FindModelEnd. If we find the end of the model in the sequence, then we 
        write the shortname and the part of the sequence after the model to the
        file "TMPFNA". So TMPFNA contains shortname of sequence and part of sequence
        after the model in a FASTA file format:
            >Name\nSequenceAfterModel\n>Name\nSeq...
        unmapped keeps track of unmapped sequences
        TRUNC is a truncated version of the fastq file with only good sequences
            from the end of the Model and the barcode is maintained in the name

    Args:
      inp_dict must contain:
        fastq_fp (str) Read from
        trunc_fp (str) Write to
        tmpFNA_fp (str) Write to
        unmapped_fp (str or None) Write to (fp)
        maxReads (int or None) (maxReads) 
        model (str) Model string
        nameToBarcode: (dict)
        minScore (int) The minimum amount of sequence beyond model needed.
        wobbleAllowed (int) uncertainty in location of barcode or end of transposon,
            on either side of expectation
        minQuality (int) every nucleotide in a barcode must be at least this quality
        flanking (int) number of nucleotides on each side that must match
        debug (bool)

    Returns:
        ret_d: (dict) containing
            nLong
            nTryToMap
            nReads
            mapnames

    """
    
    # We un-dict this variable to allow for easier reading
    nameToBarcode = inp_dict['nameToBarcode']
    
    # FILE HANDLES (Note: read, write, write)
    FASTQ = open(inp_dict['fastq_fp'], "r")
    TMPFNA = open(inp_dict['tmpFNA_fp'], "w")
    TRUNC = open(inp_dict['trunc_fp'], "w")
    

    # Current number of reads read
    nReads = 0
    # How many reads are candidates for mapping
    nTryToMap = 0
    # Current number of Long hits -
    # this is defined by: the length of the sequence vs length of model and 
    # minScore
    nLong = 0

    # Start and End of barcode within model:
    barcodeStart = inp_dict['model'].find("N")
    barcodeEnd = inp_dict['model'].rfind("N")
    barcodeLen = barcodeEnd - barcodeStart + 1


    if barcodeStart == 0:
        raise Exception("Barcode starts at beginning of model, should be in the middle.")
    for i in range(barcodeLen):
        if inp_dict['model'][barcodeStart + i] != "N":
            raise Exception("Barcode should be consecutive Ns, found a non-N in between edge 'N's."
                            " Supposed barcode range: " + inp_dict['model'][barcodeStart:barcodeEnd + 1])
   
    # Dict of read (shortname) to 1 if mapped or pastEnd, 0 otherwise.
    # mapnames only fills out when unmapped option is set to True
    mapnames = {}

   
    # The following dict is used for the "Find" functions:
    # FindBarcode & FindModelEnd. We use fastq_fp to print
    # out debug statements
    #find_info_tup = namedtuple("FindCFG", "flanking wobbleAllowed minQuality fastq_fp debug")
    #f_cfg = find_info_tup(inp_dict['flanking'], inp_dict['wobbleAllowed'], 
    #                      inp_dict['minQuality'], inp_dict['fastq_fp'],
    #                      inp_dict['debug'])
    Find_cfg_d = {
            "flanking": inp_dict['flanking'],
            "wobbleAllowed": inp_dict['wobbleAllowed'],
            "minQuality": inp_dict['minQuality'],
            "fastq_fp": inp_dict['fastq_fp'],
            "debug": inp_dict['debug']
            }

    # This is to keep track of location in FASTQ files for errors
    line_num = 0
    
    minReadLength = len(inp_dict['model']) + inp_dict['minScore']

    # We break out of while loop if file ends (name == '')
    while (inp_dict['maxReads'] is None or (nReads < inp_dict['maxReads'])):
        
        name = FASTQ.readline()
        line_num += 1


        if name == '':
            # If the file ends we break out of the loop
            break
        name = name.rstrip()
        if not name[0] == '@':
            raise Exception("Sequence name line does not start with @. File: "
            "{}, Line no. {}".format(inp_dict['fastq_fp'], line_num))

        seq = FASTQ.readline()
        line_num += 1
        seq = seq.rstrip()
        if not len(seq) > 0:
            raise Exception("Sequence line is empty. File {}, Line {}".format(
                inp_dict['fastq_fp'], line_num))
        if not re.match(r'^[A-Z]+$', seq):
            raise Exception("Sequence line contains invalid chars: " + seq \
                    + " \n File {}, Line {}".format(
                inp_dict['fastq_fp'], line_num) )


        break_line = FASTQ.readline()
        line_num += 1
        if not break_line[0] == '+':
            raise Exception("Third line does not start with +")

        quality = FASTQ.readline()
        quality = quality.rstrip()
        line_num += 1
        
        if len(quality) != len(seq):
            raise Exception("Quality line is wrong length. "
                    " File {}, Line {}".format(
                inp_dict['fastq_fp'], line_num) )

        # Ignore second side of paired-end reads
        if re.match(r'^\S+ 2:', name):
            continue

        nReads += 1

        # Short sequences are unmappable
        if len(seq) < minReadLength:
            continue
        
        nLong += 1

        
        # We keep track of location within FASTQ file for Debugging purposes
        #Find_cfg_d['line_num'] = line_num

        # obsStart is start of barcode within sequence
        # str, int. This function returns [None, None] if the 
        # quality or length fails.
        barcode, obsStart = FindBarcode(seq, quality, inp_dict['model'],
                                        barcodeStart, barcodeEnd, 
                                        Find_cfg_d, line_num)
        if barcode is None:
            continue


        # We create a shortname which removes " " to end of "name" of sequence
        # e.g. "@M00361:58:000000000-C9BPW:1:1102:11648:1000 1:N:0:GCCAAT"
        # becomes "@M00361:58:000000000-C9BPW:1:1102:11648:1000"
        shortname = name.split(' ')[0]

        if shortname in nameToBarcode:
            raise KeyError("Duplicate read name: {}\nFile {} line no. {}".format(
                            shortname, inp_dict['fastq_fp'], line_num -3))

        nameToBarcode[shortname] = barcode

        # We take start point of barcode in sequence and start point of 
        # barcode in the Model and take their difference to get where
        # the model would start in the sequence
        expModelStart = obsStart - barcodeStart
        # transposonEnd is an int (or None) which tells you where the model Ends
        # within the sequence
        transposonEnd = FindModelEnd(seq, inp_dict['model'], expModelStart,
                                    Find_cfg_d)

        if (transposonEnd is not None) and (len(seq) >= transposonEnd + inp_dict[
                'minScore']):
            # We write out to info files

            # The part of sequence in the Genome
            inGenome = seq[transposonEnd + 1:]
            TMPFNA.write(">{}\n{}\n".format(shortname, inGenome))
            nTryToMap += 1

            if inp_dict['unmapped_fp'] is not None:
                # i.e. if we are keeping track of unmapped names
                mapnames[shortname] = 0

            words = name.split(' ')
            words[0] += ":" + barcode
            TRUNC.write(" ".join(words) + "\n" + seq[transposonEnd:] \
                    + "\n+\n" + quality[transposonEnd:] + "\n")

        if line_num % 1000000 == 0:
            logging.info("So far, read {} lines".format(line_num))


    
    FASTQ.close()
    TMPFNA.close()
    TRUNC.close()
    logging.info("Read " + str(nReads) + " reads\n")
    
    if nTryToMap == 0:
        logging.warning("None of the reads are candidates for mapping."
                " None match the model and are long enough.\n")
        os.unlink(inp_dict['tmpFNA_fp'])
        return -1 

    ret_d = {
            "nReads": nReads,
            "nTryToMap": nTryToMap,
            "nLong": nLong,
            "mapnames": mapnames,
            "nameToBarcode": nameToBarcode,
            "minReadLength": minReadLength
            }

    return ret_d

    

def FindBarcode(seq, quality, model, expStart, expEnd, cfg_d, line_num) -> Optional[Tuple[str,int]]:
    """
    Args:
        seq (str) DNA sequence
        quality (str) quality of sequence
        model (str) model string
        expStart (int) Expected Start of barcode within sequence (It is the start of barcode within the model)
        expEnd (int) Expected end of barcode within sequence (It is the end of the barcode within the model)
        cfg_d: (dict) (config_dict)
            wobbleAllowed (int) uncertainty in location of barcode or end of transposon,
                on either side of expectation
            minQuality (int) every nucleotide in a barcode must be at least this quality
            flanking (int) number of nucleotides on each side that must match
            debug (bool)
            fastq_fp: (str) The path to the FASTQ file
        line_num: (int) line number within FASTQ file of quality 

    Returns:
        list<barcode, obsStart>
            barcode: (str) sequence of length ~20 which is the random barcode
            obsStart: (int) position of first nucleotide of barcode within seq
    """

    # undict for readability
    flanking = cfg_d['flanking']

    pre = model[expStart - flanking:expStart]
    # We search for "pre" within the sequence at the location 'barcode Start' within model - flanking
    preLoc = FindSubstr(pre, seq, expStart - flanking, cfg_d['wobbleAllowed'])
    if preLoc is None:
        return [None, None]

    post = model[expEnd + 1: expEnd + 1 + flanking]

    # We search for "post" within the sequence at the location 'barcode End' within model - flanking
    postLoc = FindSubstr(post, seq, expEnd + 1, cfg_d['wobbleAllowed'])
    if postLoc is None:
        return [None, None]

    # positions of 1st and last nucleotides of barcode
    start = preLoc + flanking
    end = postLoc - 1

    # expected End and Start are related to loc. of barcode in Model
    if not end - start == expEnd - expStart:
        if cfg_d['debug']:
            logging.info("Wrong barcode length: {} {} not {} {} in {}\n".format(
            start, end, expStart, expEnd, seq))
        return [None, None]

    barcode = seq[start:end + 1]

    if cfg_d['minQuality'] > 0:
        barqual = quality[start: end + 1]
        #quality score is encoded as 33
        scores = [ord(c) for c in barqual]
        for score in scores:
            if score < 0 or score > 100:
                raise Exception("Invalid score {} from barcode {}".format(
                score, barcode) + " quality {}".format(
                barqual) + " in file {} line {}".format(
                    cfg_d['fastq_fp'], line_num))
            if score < cfg_d['minQuality']:
                if cfg_d['debug']:
                    logging.info("Low quality {} for barcode {} in {}\n".format(
                        score, barcode, seq) + "File {} line {}".format(
                            cfg_d['fastq_fp'], line_num))
                return [None, None]
    
    return [barcode, start] # str, int


def FindSubstr(subseq, seq, expAt, wobble):
    """
    Description:
        Here we search for model parts pre and post barcode within a sequence
        from reads.

    Args:
        subseq (str) Subsequence within string we're looking for
        seq (str) Full sequence to search subseq within
        expAt (int) expected at location of subseq
        wobble (int) possible interval around expAt for which subseq can be found.
    """
    L = len(seq)
    K = len(subseq)
    for i in range(expAt - wobble, expAt + wobble + 1):
        if i >= 0 and i < L and seq[i:i+K] == subseq:
            return i

    return None



def FindModelEnd(seq, model, expOff, cfg_d):
    """
    Args: 

    seq: (str) String of Sequence from FASTQ
    model: (str) String of entire model (no PastEnd)
    expOff: (int) expected location model starts
    cfg_d:
        wobbleAllowed (int) uncertainty in location of barcode or end of transposon,
            on either side of expectation
        minQuality (int) every nucleotide in a barcode must be at least this quality
        flanking (int) number of nucleotides on each side that must match
        debug: (bool)
        fastq_fp: (str) The name of the FASTQ file
        line_num: (int) line number within FASTQ file of quality 

    Returns:
    None or transposonEnd (int) index in seq where transposon ends.

    """
    # expected End index of model is where it starts + length of model
    # expOff could be 0
    expEnd = expOff + len(model)

    # We look for the end of the model within the sequence
    # at is an int if found, None if not found
    at = FindSubstr(model[-cfg_d['flanking']:],seq, expEnd - cfg_d['flanking'],
                    cfg_d['wobbleAllowed'])

    if at is None:
        if cfg_d['debug']:
            logging.info("No end sequence {} near position {} of \n{}\n".format(
                model[-cfg_d['flanking']:], expEnd - cfg_d['flanking'],
                seq) + " in file {} line no. {}".format(
                    cfg_d['fastq_fp'], cfg_d['line_num']))
        return None
    
    # last position of transposon
    return at + cfg_d['flanking'] - 1



def RunBLAT8(inp_dict):
    """
    Run this UCSC tool on a few cases:
    Case A:
        The queries file is full of end of sequences that had the  barcodes in
            them (tmpFNA sequences).
        The database is the pastEnd sequence (only two lines)
    Case B:
        The queries is the tmpFNA sequences - the barcodes reads
        The database is the genome FASTA file

    Args:
    
    inp_dict: (dict) must contain
        queriesFile (str) In case A & B: (tmpFNA_fp) 
        dbFile (str) case A: (endFNA_fp) caseB: (genome FASTA File)
        blat8_fp: (str) Path to output file
        blatcmd: (str) Path to BLAT command, i.e. Bin/blat
        minScore: (int) minimum score for mapping to genome or past-end
        minIdentity: (int) minimum %identity for mappinf to genome or past-end
        tileSize: (int) size of an alignment tile
        stepSize: (int) distance between the starting bases of alignment tiles
            (will overlap if stepSize<tileSize)
    

    """
    
    blat_args = [
            inp_dict['blatcmd'], 
            "-out=blast8", 
            "-t=dna",
            "-minScore=" + str(inp_dict['minScore']),
            "-minIdentity=" + str(inp_dict['minIdentity']),
            "-maxIntron=0",
            "-noTrimA",
            "-tileSize=" + str(inp_dict['tileSize']),
            "-stepSize=" + str(inp_dict['stepSize']),
            inp_dict['dbFile'],
            inp_dict['queriesFile'],
            inp_dict['blat8_fp']
            ]

    logging.info("Running blat8 with following command:\n{}".format(
        " ".join(blat_args)))
    blat8_result = subprocess.run(blat_args)
    logging.info("blat8 result: " + str(blat8_result))

    return None 


def HandleGenomeBLAT(rows, hitsPastEnd, HG_d, debug):
    """
    Description:
        This function writes out reads to pre-pool table file
        Each time this is run is on a single read name- which may correspond to a number
            of rows in the blat output. If there are multiple hits, the length of rows 
            input will be > 1, otherwise it is 1.
            The rows are sorted in decreasing bit score order (already by blat8).
        IF pastEnd is in the model, then we will have already run blat8 on all the tmpFNA sequences
            and we compare the scores to the pastEnd scores, which are in "hitsPastEnd" dict.
        We compare the scores for the queries with the PastEnd to the genome. "Trumping"
            means that the score for one is greater than the other (*no political affiliation
            implied in the term). We update hitsPastEnd dict for the read to '0' if the 
            best genomic hit score trumps the hit-past-end score. 

    Args:
        rows: list of lists, each sub list a blat8 split line (split by tab)
        hitsPastEnd: (dict) (shortname to score against PastEnd sequence)
            hitsPastEnd keeps track of "inGenome" sequences that are close 
            or the same as the pastEnd sequence from the plasmid.
        debug: (bool)
        HG_d (dict): Handle Genome Blat dict
            nPastEndTrumps: (int)  hit to past-end (close to) as good as hit to genome
            nPastEndIgnored: (int) weak hit to past-end ignored 
            nMapUnique: (int)
            nMapped: (int)
            nameToBarcode: (dict)
                short_name (str) -> barcode (str)
            delta: (int) minimum difference in score for considering a hit unique.
            OUTPUT: file handle for output

    """
    if len(rows) == 0:
        return None

    # read is a string of a readname, i.e. '@K00364:102:HWV25BBXX:1:1101:5335:1314'
    read = rows[0][0]

    if read == '' or read not in HG_d['nameToBarcode']:
        logging.info(rows)
        logging.info(HG_d['nMapped'])
        error_msg = "BLAT Compare Genome to TmpFNA Read Name error: " + read
        error_msg += "\nThis is an internal error running BLAT8, Contact helpdesk."
        raise RuntimeError(error_msg)


    # indexes within besthits entries lists:
    SCAFFOLD, POSITION, STRAND, SCORE, IDENTITY, QBEG, QEND = range(7)

    # besthits is a list of lists, where each list corresponds to a row from reads
    # if the score passes a threshold, each sublist is 
    # [subject name (scaffold), subject Begin, strand, score, identity, 
    # query Begin, query End]
    besthits = []

    for row in rows:
        # read2, subject, identity, length, mismatches, num gaps,
        #   query Begin index, query End index, subject Begin index,
        #   subject End index, e-value, bit score (12 in total).
        read2, sbj, idn, l, mm, gaps, qB, qE, sB, sE, evl, score = row
        if read2 != read:
            raise Exception("new read not equal: \n{}!=\n{}".format(read2, read))
        if len(besthits) == 0 or (
                float(score) >= float(besthits[0][SCORE]) - HG_d['delta']):

            # Often delta is as small as 5, so if the difference is greater than 5 between
            # the top scorer and the one checking, there will only be 1 best hit

            # convert from 0-based to 1-based pos, and note that sB always < sE
            # so flip if stranded
            strand = "+" if int(sB) < int(sE) else "-"

            besthits.append([sbj, sB, strand, score, idn, qB, qE])
            if debug:
                logging.info("Hit for {}:{}:{} to {}:{}:{} score {}".format(
                    read, qB, qE, sbj, sB, sE, score) + " identity {}".format(
                        idn))

    # Scores always ordered from greatest to least in BLAT8 output file 
    # This is always the value with the best bit score
    # pos is position in genome scaffold
    scaffold, pos, strand, score, idn, qB, qE = besthits[0]

    # and output a mapping row (or none)
    if read in hitsPastEnd:
        if hitsPastEnd[read] >= float(score) - HG_d['delta']:
            HG_d['nPastEndTrumps'] += 1
            if debug:
                logging.info("Past end trumps for " + read)
            return None
        else:
            HG_d['nPastEndIgnored'] += 1
            if debug:
                logging.info("Ignoring hit past end of transposon for " + read)
            # This is so we don't print out a line for it later on...
            hitsPastEnd[read] = 0.0

    # else if no trumping, we write to output (MTS Table)
    # read name, barcode, scaffold (from besthit (genome.fna)), position in scaffold,
    # strand (+ or -), uniqueness: whether or not there are multiple loc for good reads,
    # query Beginning, query End, bit score, % identity 
    HG_d['OUTPUT'].write("\t".join([
        read, 
        HG_d['nameToBarcode'][read],
        scaffold,
        pos, 
        strand,
        "1" if len(besthits) == 1 else "0",
        qB,
        qE,
        score,
        idn + "\n"
        ]))

    if len(besthits) == 1:
        HG_d['nMapUnique'] += 1

    HG_d['nMapped'] += 1

    return None






def test():

    return None


def main():
    test()

    return None

if __name__ == "__main__":
    main()

