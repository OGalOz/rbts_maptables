# -*- coding: utf-8 -*-
import os
import time
import unittest
import logging
import copy
from configparser import ConfigParser

from rbts_maptables.rbts_maptablesImpl import rbts_maptables
from rbts_maptables.rbts_maptablesServer import MethodContext
from rbts_maptables.authclient import KBaseAuth as _KBaseAuth

from installed_clients.WorkspaceClient import Workspace


class rbts_maptablesTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = os.environ.get('KB_AUTH_TOKEN', None)
        config_file = os.environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('rbts_maptables'):
            cls.cfg[nameval[0]] = nameval[1]
        # Getting username from Auth profile for token
        authServiceUrl = cls.cfg['auth-service-url']
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'rbts_maptables',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = Workspace(cls.wsURL)
        cls.serviceImpl = rbts_maptables(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']
        suffix = int(time.time() * 1000)
        cls.wsName = "mapTnSeq_Test_" + str(suffix)
        cls.MTS_Test_Defaults = {
           'workspace_name': cls.wsName,
           'maxReads': None,
           'minQuality': 5,
           'minIdentity': 90,
           'minScore': 15,
           'delta': 5
        }
        ret = cls.wsClient.create_workspace({'workspace': cls.wsName})  # noqa

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    # NOTE: According to Python unittest naming rules test method names should start from 'test'. # noqa
    def test1(self):
        ##
        ## Default Test (Quick)
        ##
        #"" 
        # Dyella Japonica
        genome_ref = "63063/9/1" 
        fastq_ref_list = ["63063/7/1"]
        tnseq_model_name = "pKMW3_universal"
        pool_description = "Testing"
        output_name = "Test1"
        test_d = copy.deepcopy(self.MTS_Test_Defaults)
        test_d["genome_ref"] = genome_ref
        test_d["fastq_ref_list"] = fastq_ref_list
        test_d["tnseq_model_name"] = tnseq_model_name
        test_d["output_name"] = output_name 
        #test_d["yy"] = yy 
        ret = self.serviceImpl.run_rbts_maptables(self.ctx, test_d)
        # Check ret:
        logging.info("Finished running test 1. Results:")
        logging.info(ret)

    def test2(self):
        ## Test Purpose:
        ##     Multiple FASTQ files
        ##     & minQuality = 0 
        #""

        genome_ref = "63063/9/1" 
        fastq_ref_list = ["63063/7/1", "63063/11/1", "63063/13/1"]
        tnseq_model_name = "pKMW3_universal"
        pool_description = "Testing"
        output_name = "Test1"
        test_d = copy.deepcopy(self.MTS_Test_Defaults)
        test_d["genome_ref"] = genome_ref
        test_d["fastq_ref_list"] = fastq_ref_list
        test_d["tnseq_model_name"] = tnseq_model_name
        test_d["output_name"] = output_name 
        test_d["minQuality"] = 0 
        #test_d["yy"] = yy 
        ret = self.serviceImpl.run_rbts_maptables(self.ctx, test_d)
        # Check ret:
        logging.info("Finished running test 2. Results:")
        logging.info(ret)
        #""
        #pass
    def test3(self):
        ## Test Purpose:
        ##    New Genome (E Coli Keio) and FASTQs
        ##    & minIdentity = 99

        # E Coli
        genome_ref = "63063/3/1" 
        fastq_ref_list = ["63063/2/1"]
        tnseq_model_name = "Sc_Tn5"
        pool_description = "Testing"
        output_name = "Test1"
        test_d = copy.deepcopy(self.MTS_Test_Defaults)
        test_d["genome_ref"] = genome_ref
        test_d["fastq_ref_list"] = fastq_ref_list
        test_d["tnseq_model_name"] = tnseq_model_name
        test_d["output_name"] = output_name 
        test_d["minIdentity"] = 99
        #test_d["yy"] = yy 
        ret = self.serviceImpl.run_rbts_maptables(self.ctx, test_d)
        # Check ret:
        logging.info("Finished running test 3. Results:")
        logging.info(ret)
        #pass

    def test4(self):
        ## Test Purpose:
        ##    New Genome (E Coli Keio) and FASTQs
        ##    & minIdentity = 1

        # E Coli
        genome_ref = "63063/3/1" 
        fastq_ref_list = ["63063/2/1"]
        tnseq_model_name = "Sc_Tn5"
        pool_description = "Testing"
        output_name = "Test1"
        test_d = copy.deepcopy(self.MTS_Test_Defaults)
        test_d["genome_ref"] = genome_ref
        test_d["fastq_ref_list"] = fastq_ref_list
        test_d["tnseq_model_name"] = tnseq_model_name
        test_d["output_name"] = output_name 
        test_d["minIdentity"] = 1
        #test_d["yy"] = yy 
        ret = self.serviceImpl.run_rbts_maptables(self.ctx, test_d)
        # Check ret:
        logging.info("Finished running test 4. Results:")
        logging.info(ret)
        #pass
    def test5(self):
        ## Test Purpose:
        ## Break early no genome_ref

        # E Coli
        genome_ref = "63063/3/1" 
        fastq_ref_list = ["63063/2/1"]
        tnseq_model_name = "Sc_Tn5"
        pool_description = "Testing"
        output_name = "Test1"
        test_d = copy.deepcopy(self.MTS_Test_Defaults)
        test_d["fastq_ref_list"] = fastq_ref_list
        test_d["tnseq_model_name"] = tnseq_model_name
        test_d["output_name"] = output_name 
        test_d["minIdentity"] = 1
        #test_d["yy"] = yy 
        ret = self.serviceImpl.run_rbts_maptables(self.ctx, test_d)
        # Check ret:
        logging.info("Finished running test 5. Should have failed:")
        logging.info(ret)
    """
    def test_your_method(self):
        # Prepare test objects in workspace if needed using
        # self.getWsClient().save_objects({'workspace': self.getWsName(),
        #                                  'objects': []})
        #
        # Run your method by
        # ret = self.getImpl().your_method(self.getContext(), parameters...)
        #
        # Check returned data with
        # self.assertEqual(ret[...], ...) or other unittest methods
        ret = self.serviceImpl.run_rbts_maptables(self.ctx, {'workspace_name': self.wsName,
                                                             'parameter_1': 'Hello World!'})
    """
