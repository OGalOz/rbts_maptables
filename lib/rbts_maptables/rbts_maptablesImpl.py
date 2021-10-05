# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.GenomeFileUtilClient import GenomeFileUtil
from installed_clients.AssemblyUtilClient import AssemblyUtil
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace
from installed_clients.rbts_genome_to_genetableClient import rbts_genome_to_genetable
from full.FullProgram import CompleteRun
from util.PrepareIO import PrepareProgramInputs, PrepareUserOutputs, clear_dir
from installed_clients.KBaseReportClient import KBaseReport
#END_HEADER


class rbts_maptables:
    '''
    Module Name:
    rbts_maptables

    Module Description:
    A KBase module: rbts_maptables
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = ""
    GIT_COMMIT_HASH = ""

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.shared_folder = config['scratch']
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        self.ws_url = config['workspace-url']
        #END_CONSTRUCTOR
        pass


    def run_rbts_maptables(self, ctx, params):
        """
        Args:
            params (d) contains keys:
                'workspace_name' (str): ,
                'genome_ref' (str): ,
                'tnseq_model_name' (str): ,
                'fastq_ref_list' list<str>: ,
                'maxReads' (int): ,
                'minQuality' (int): ,
                'minIdentity' (int): ,
                'minScore' (int): ,
                'delta' (int): ,
                'output_name' (str)
        This example function accepts any number of parameters and returns results in a KBaseReport
        :param params: instance of mapping from String to unspecified object
        :returns: instance of type "ReportResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_rbts_maptables

        logging.basicConfig(level=logging.DEBUG)
        if len(os.listdir(self.shared_folder)) > 0:
            logging.info("Clearing scratch directory")
            clear_dir(self.shared_folder)

        # Preparing main classes - dfu, gfu, genetableobj, workspace
        dfu = DataFileUtil(self.callback_url)
        gfu = GenomeFileUtil(self.callback_url)
        genetable_obj = rbts_genome_to_genetable(self.callback_url)
        # We need the workspace object to get info on the workspace the app is running in.
        token = os.environ.get('KB_AUTH_TOKEN', None)
        ws = Workspace(self.ws_url, token=token)
        ws_info = ws.get_workspace_info({'workspace': params['workspace_name']})
        workspace_id = ws_info[0]
        td = self.shared_folder

        # Initializing config info 
        cfg_d = {
                "gfu": gfu,
                "dfu": dfu,
                "ws": ws,
                "ws_id": workspace_id,
                "gt_obj": genetable_obj,
                "username" : ctx['user_id'],
                "tmp_dir": td,
                "model_dir": td, 
                "gene_table_fp": os.path.join(td, "genes.GC"),
                "mts_tables_dir": os.path.join(td, "MTS_Tables"),
                "blat_cmd": "/kb/module/lib/rbts_maptables/blat"
        }


        # We divide the program into 3 parts:
        # Part 1: Prepare to run program: Download necessary files, create configs
        vp, genome_scientific_name, MTS_cfg_d = PrepareProgramInputs(params, cfg_d)


        # Part 1.5: Prepare vars 
        cfg_d["workspace_name"] = params["workspace_name"]
        #cfg_d["Main_HTML_report_fp"] = html_fp

        # Part 2: Run the central part of the program using recently created config files
        pre_HTML_d = CompleteRun(MTS_cfg_d, 
                                 cfg_d["tmp_dir"], genome_scientific_name,
                                 cfg_d, vp)

        logging.debug(os.listdir(td))

        PrepareUserOutputs(vp, cfg_d)


        report = KBaseReport(self.callback_url)
        report_info = report.create({'report': {'objects_created':[],
                                                'text_message': ""},
                                                'workspace_name': params['workspace_name']})
        output = {
            'report_name': report_info['name'],
            'report_ref': report_info['ref'],
        }
        #END run_rbts_maptables

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_rbts_maptables return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
