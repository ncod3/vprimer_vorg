#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import errno
import re

# global variants
import vprimer.glv as glv
import vprimer.utils as utl

# https://qiita.com/mimitaro/items/3506a444f325c6f980b2
import configparser
# https://docs.python.org/ja/3/library/configparser.html

# class
from vprimer.param import Param
from vprimer.logging_config import LogConf


class Conf(object):

    def __init__(self):

        self.init_version = 0.0

        self.ini_file = ''
        self.ini_file_path = ''
        self.ini = configparser.ConfigParser()
        # don't convert to lower case
        self.ini.optionxform = str

        # log
        self.log = LogConf()

        self.thread = 0

        self.parallel_blast_cnt = 0
        self.parallele_full_thread = 0
        self.blast_num_threads = 0

        # default 1000
        self.blast_word_size = 0
        self.use_joblib_threading = 'no'
        self.parallel = False

        # indel_len
        self.min_indel_len = 0
        self.max_indel_len = 0
        # product_size
        self.min_product_size = 0
        self.max_product_size = 0

        # group
        self.regions_dict = dict()
        self.distin_g_list = list()
        self.g_members_dict = dict()

        self.cwd = os.getcwd()

        self.ref_dir = ''
        self.out_dir = ''
        self.out_bak_dir = ''
        self.log_dir = ''

        self.pick_mode = ''
        self.progress = ''
        self.stop = 'no'

        self.ref = ''
        self.vcf = ''

        self.min_indel_len = 0
        self.max_indel_len = 0
        self.min_product_size = 0
        self.max_product_size = 0
        self.fragment_pad_len = 0

        self.PRIMER_THERMODYNAMIC_PARAMETERS_PATH = ''
        self.PRIMER_PRODUCT_SIZE_RANGE = ''
        self.PRIMER_NUM_RETURN = ''
        self.PRIMER_MIN_SIZE = ''
        self.PRIMER_OPT_SIZE = ''
        self.PRIMER_MAX_SIZE = ''
        self.PRIMER_MIN_GC = ''
        self.PRIMER_OPT_GC = ''
        self.PRIMER_MAX_GC = ''
        self.PRIMER_MIN_TM = ''
        self.PRIMER_OPT_TM = ''
        self.PRIMER_MAX_TM = ''
        self.PRIMER_MAX_POLY_X = ''
        self.PRIMER_PAIR_MAX_DIFF_TM = ''

        self.alternate_distance = 0

        # ------------------------
        self.ref_fasta = ''
        self.ref_fasta_chrom_list = ''
        self.ref_fasta_fai = ''
        self.ref_fasta_pickle = ''
        self.ref_fasta_slink_system = ''
        self.ref_fasta_user = ''

        self.vcf_file_user = ''
        self.vcf_file_slink_system = ''
        self.vcf_file = ''

        #-- sample_nickname
        self.vcf_basename_to_fullname = dict()
        self.vcf_nickname_to_fullname = dict()
        self.nickname_to_basename = dict()
        self.basename_to_nickname = dict()

        #-- enzyme
        self.enzyme_files_list = list()
        self.enzyme_files_str = ''

        self.blastdb_title = ''
        self.blastdb = ''
        # ------------------------


    def read_ini(self):

        # ---------------------------------------------------------------
        # ini_file from param.p['config'] absolute or relative path
        self.ini_file = glv.param.p.config
        print("ini_file = {}".format(self.ini_file))
        # ---------------------------------------------------------------

        # ini_file_path
        self.ini_file_path = "{}/{}".format(
            self.cwd, self.ini_file)
        print("self.ini_file_path = {}".format(self.ini_file_path))

        # read ini_file
        if os.path.exists(self.ini_file_path):
            print("found {}".format(self.ini_file_path))

            # https://docs.python.org/ja/3/library/configparser.html
            with open(self.ini_file_path, encoding='utf-8') as fp:
                self.ini.read_file(fp)
 
                # adjustment of variable format
                self._rectify_variable()

                #
                self._ini_into_variable()

                # several path
                self._set_paths()

                # thread
                self._thread_adjusting()

                # nickname
                self._get_nickname()
                # regions
                self._get_regions()
                # distin_g
                self._get_distinguish_groups()
                # g_members
                self._get_members()

        else:
            print("Not found {}, exit.".format(self.ini_file_path))
            #raise FileNotFoundError(errno.ENOENT,
            #    os.strerror(errno.ENOENT),
            #    self.ini_file_path)
            sys.exit(1)

        return self


    def _ini_into_variable(self):

        sect = 'global'

        self.init_version = \
            float(self._set_default(sect, 'init_version', 1.0))

        self.use_joblib_threading = \
             str(self._set_default(sect, 'use_joblib_threading', 'yes'))

        self.parallel = \
            self._conv_joblib(self.use_joblib_threading)

        #log.info("{} {} {}".format(
        #    self.ini[sect]['use_joblib_threading'],
        #    self.use_joblib_threading,
        #    self.parallel))

        self.thread = \
            int(self._set_default(sect, 'thread', 2))

        self.ref = \
            str(self._set_default(sect, 'ref', ''))
        self.vcf = \
            str(self._set_default(sect, 'vcf', ''))

        self.min_indel_len = \
            int(self._set_default(sect, 'min_indel_len', 50))
        self.max_indel_len = \
            int(self._set_default(sect, 'max_indel_len', 200))
        self.min_product_size = \
            int(self._set_default(sect, 'min_product_size', 200))
        self.max_product_size = \
            int(self._set_default(sect, 'max_product_size', 500))

        self.ref_dir = \
            str(self._set_default(sect, 'ref_dir', 'refs'))
        self.log_dir = \
            str(self._set_default(sect, 'log_dir', 'logs'))
        self.out_dir = \
            str(self._set_default(sect, 'out_dir', 'out_dir'))
        self.out_bak_dir = ''

        self.pick_mode = \
            str(self._set_default(sect, 'pick_mode', 'all'))
        self.progress = \
            str(self._set_default(sect, 'progress', 'all'))

        #sect = 'regions'
        #sect = 'sample_nickname'
        #sect = 'groups'

        sect = 'vprimer'
        self.fragment_pad_len = \
            int(self._set_default(sect, 'fragment_pad_len', 500))

        sect = 'caps'

        self.enzyme_files_str = \
            (self._set_default(sect, 'enzyme_files', ''))

        sect = 'primer3'

        self.PRIMER_THERMODYNAMIC_PARAMETERS_PATH = \
            str(self._set_default(
                sect, 'PRIMER_THERMODYNAMIC_PARAMETERS_PATH', ''))
        self.PRIMER_PRODUCT_SIZE_RANGE = "{}-{}".format(
            self.min_product_size, self.max_product_size)

        self.PRIMER_NUM_RETURN = str(1)

        self.PRIMER_MIN_SIZE = \
            str(self._set_default(sect, 'PRIMER_MIN_SIZE', 23))

        self.PRIMER_OPT_SIZE = \
            str(self._set_default(sect, 'PRIMER_OPT_SIZE', 25))
        self.PRIMER_MAX_SIZE = \
            str(self._set_default(sect, 'PRIMER_MAX_SIZE', 27))
        self.PRIMER_MIN_GC = \
            str(self._set_default(sect, 'PRIMER_MIN_GC', 40))
        self.PRIMER_OPT_GC = \
            str(self._set_default(sect, 'PRIMER_OPT_GC', 50))
        self.PRIMER_MAX_GC = \
            str(self._set_default(sect, 'PRIMER_MAX_GC', 60))
        self.PRIMER_MIN_TM = \
            str(self._set_default(sect, 'PRIMER_MIN_TM', 57.0))
        self.PRIMER_OPT_TM = \
            str(self._set_default(sect, 'PRIMER_OPT_TM', 60.0))
        self.PRIMER_MAX_TM = \
            str(self._set_default(sect, 'PRIMER_MAX_TM', 63.0))
        self.PRIMER_MAX_POLY_X = \
            str(self._set_default(sect, 'PRIMER_MAX_POLY_X', 4))
        self.PRIMER_PAIR_MAX_DIFF_TM = \
            str(self._set_default(
                sect, 'PRIMER_PAIR_MAX_DIFF_TM', 4))

        sect = 'blast'

        self.alternate_distance = \
            int(self._set_default(sect, 'alternate_distance', 10000))
        self.blast_word_size = \
            int(self._set_default(
                sect, 'blast_word_size', self.PRIMER_MIN_SIZE))


    def _set_default(self, sect, key, value):

        if key in self.ini[sect]:
            return self.ini[sect][key]
        else:
            return value


    def _conv_joblib(self, string):

        bool_joblib = False
        if string == 'yes':
            bool_joblib = True

        return bool_joblib


    def _thread_adjusting2(self):

        #log.info("thread={} parallel={}".format(
        #    self.thread, self.parallel))

        self.parallel = \
            self._conv_joblib(self.use_joblib_threading)

        if self.thread < 6:
            self.parallel = False

        if self.parallel == True:

            # unit is 4+1=5
            parallel_base = self.thread
            self.parallele_full_thread = parallel_base

            # blast = 4
            self.parallel_blast_cnt = int(parallel_base / 5)
            self.blast_num_threads = 4

        else:
            # 6 = 5
            full_thread = self.thread - 1
            self.parallele_full_thread = full_thread
            self.blast_num_threads = full_thread

        #log.info("thread={}, parallel={}, parallel_blast_cnt={}".format(
        #    self.parallel, self.thread, self.parallel_blast_cnt))
        #log.info("parallele_full_thread={}, blast_num_threads={}".format(
        #    self.parallele_full_thread, self.blast_num_threads))


    def _thread_adjusting(self):
        ''' in Palallel, if there are 10 threads blast cmd will use at least
            2 cores so par 1 parallel.
            So main python always use 1,
            parallel use 1 thread, blast use 2 threads
        '''

            #self.parallel_blast_cnt = 1
            #self.parallele_full_thread = 1 
            #self.blast_num_threads = 1

        #self._thread_adjusting2()
        #return

        #log.info("thread={} use_joblib_threading={}, parallel={}".format(
        #    self.thread, self.use_joblib_threading, self.parallel))

        self.parallel = \
            self._conv_joblib(self.use_joblib_threading)

        if self.thread < 6:
            self.parallel = False

        if self.parallel == True:
            #(7) = 6
            parallel_base = self.thread
            self.parallele_full_thread = parallel_base

            # 6/3=2
            self.parallel_blast_cnt = int(parallel_base / 3)
            self.blast_num_threads = 2

        else:
            # 6 = 5
            full_thread = self.thread - 1
            self.parallele_full_thread = full_thread
            self.blast_num_threads = full_thread

        #log.info("thread={}, use_joblib_threading={}".format(
        #    self.thread, self.use_joblib_threading))
        #log.info("parallel={}, parallel_blast_cnt={}".format(
        #    self.parallel, self.parallel_blast_cnt))
        #log.info("parallele_full_thread={}, blast_num_threads={}".format(
        #    self.parallele_full_thread, self.blast_num_threads))


    def _get_nickname(self):

        # [sample_nickname]
        # nickname    basename
        # sample1   = sample1_sorted.bam
        # sample2   = sample1_sorted.bam

        for nickname in self.ini['sample_nickname']:
            basename = self.ini['sample_nickname'][nickname]

            # fullname will filled in vcf_file _get_sample_name_list
            self.nickname_to_basename[nickname] = basename
            self.basename_to_nickname[basename] = nickname

            #log.debug("nickname={} => basename={}".format(
            #    nickname, basename)) 


    def merge_conf(self):

        # debug by param
        if glv.param.p.stop != None:
            glv.conf.stop = str(glv.param.p.stop)

        # update by param
        if glv.param.p.thread != None:
            glv.conf.thread = int(glv.param.p.thread)
            self.use_joblib_threading = 'yes'
            self._thread_adjusting()

        if glv.param.p.joblib_threading != None:
            self.use_joblib_threading = glv.param.p.joblib_threading
            self._thread_adjusting()

        if glv.param.p.out_dir != None:
            self.out_dir = str(glv.param.p.out_dir)
            self._set_paths()

        if glv.param.p.ref != None:
            self.ref = str(glv.param.p.ref)

        if glv.param.p.vcf != None:
            self.vcf = str(glv.param.p.vcf)

        if glv.param.p.progress != None:
            self.progress = str(glv.param.p.progress)

        # indel_len
        if glv.param.p.min_indel_len != None:
            self.min_indel_len = int(glv.param.p.min_indel_len)

        if glv.param.p.max_indel_len != None:
            self.max_indel_len = int(glv.param.p.max_indel_len)

        # product_size
        if glv.param.p.min_product_size != None:
            self.min_product_size = int(glv.param.p.min_product_size)

        if glv.param.p.max_product_size != None:
            self.conf.max_product_size = int(glv.param.p.max_product_size)

        self.PRIMER_PRODUCT_SIZE_RANGE = "{}-{}".format(
            self.min_product_size, self.max_product_size)

        #log.info("{}".format(glv.param.p))

        self._print_setting()

        return self


    def _print_setting(self):

        pass
        #log.info("{}={}".format('ini_file', self.ini_file))
        #log.info("{}={}".format('ini_file_path', self.ini_file_path))
        #log.info("{}={}".format('thread', self.thread))
        #log.info("{}={}".format('use_joblib_threading', \
        #    self.use_joblib_threading))
        #log.info("{}={}".format('parallel', \
        #    self.parallel))
        #log.info("{}={}".format('parallel_blast_cnt', \
        #    self.parallel_blast_cnt))
        #log.info("{}={}".format('parallele_full_thread', \
        #    self.parallele_full_thread))
        #log.info("{}={}".format('blast_num_threads', \
        #    self.blast_num_threads))
        #log.info("{}={}".format('blast_word_size', \
        #    self.blast_word_size))
        #log.info("{}={}".format('PRIMER_PRODUCT_SIZE_RANGE', \
        #    self.PRIMER_PRODUCT_SIZE_RANGE))


    def _set_paths(self):

        #---------------------------------------------------
        # out_dir
        if glv.param.p.out_dir != None:
            self.ini['global']['out_dir'] = glv.param.p.out_dir

        # result out dir
        self.out_dir =  "{}/{}".format(
            self.cwd, self.ini['global']['out_dir'])

        # logs dir under out dir
        self.log_dir = "{}/{}".format(
            self.out_dir, self.ini['global']['log_dir'])

        # system reference dir
        self.ref_dir = "{}/{}".format(
            self.cwd, self.ini['global']['ref_dir'])

        # out_bak_dir
        self.out_bak_dir =  "{}/{}".format(
            self.out_dir, 'bak')

        # make dir
        self._make_dir_tree()

        # set to LogConf
        global log
        log = self.log.start(__name__, self.out_dir, self.log_dir)

        # log for utils
        utl.start_log()


    def _make_dir_tree(self):

        # already made at conf
        dirs = [
            self.out_dir,
            self.log_dir,
            self.out_bak_dir,
            self.ref_dir,
        ]

        for dir in dirs:
            os.makedirs(dir, exist_ok=True) # refs

    
    def _get_regions(self):
    
        # chrom:1-1000
        for rg in self.ini['regions']:
            rg_list = self.ini['regions'][rg].split(':')
    
            if rg_list[0] == '':
                #log.error("error, region={} is empty, exit.".format(rg))
                sys.exit(1)
    
            if len(rg_list) == 1:
                self.regions_dict[rg] = {
                    'chr': rg_list[0],
                    'start': None,
                    'end': None,
                    'reg': self.ini['regions'][rg]}
    
            else:
                pos = rg_list[1].split('-')
                self.regions_dict[rg] = {
                    'chr': rg_list[0],
                    'start': int(pos[0]),
                    'end': int(pos[1]),
                    'reg': self.ini['regions'][rg]}
    
            #log.info("{}: {}".format(rg, self.regions_dict[rg]))
    
    
    def _get_distinguish_groups(self):
    
        distins = self.ini['groups']['distinguish_groups'].split(';')
        #log.debug("distins={}".format(distins))
    
        for distin in distins:
    
            #log.debug("distin={}".format(distin))
            # [global][pick_mode]
            pick_mode = self.ini['global']['pick_mode']
    
            # group / reg1, reg2 : pick_mode
            g_pick_mode = distin.split(':')
            if len(g_pick_mode) == 1:
                # not exist pick_mode, only group_regions
                g_region = g_pick_mode
                #log.debug("(1) pick_mode don't exist: {}".format(distin))
    
            else:
                # exist pick_mode
                if len(g_pick_mode[1]) != 0:
                    pick_mode = g_pick_mode[1]
                    #log.debug("(3) pick_mode exist: {}".format(distin))
                #else:
                    #log.debug("(2) pick_mode empty: {}".format(distin))
    
            # split into group and regions
            g_region = g_pick_mode[0].split('/')
    
            #log.debug("g_region: {}".format(g_region))
    
            if len(g_region) == 1:
                # there is no region, only groups separated by <>
                groups = g_region[0].split('<>')
                if groups[0] == '' or groups[1] == '':
                    #log.error(
                    #    "error, empty distinguish_groups {} exit.".format(
                    #    distin))
                    sys.exit(1)
    
                g_dict = {
                    0: groups[0],
                    1: groups[1],
                    'region': [''],
                    'pick_mode': pick_mode}
    
            else:
                # all staff gathered
                groups = g_region[0].split('<>')
                if groups[0] == '' or groups[1] == '':
                    #log.error(
                    #    "error, empty distinguish_groups {} exit.".format(
                    #    distin))
                    sys.exit(1)
    
                regions = g_region[1].split(',')
    
                # gGroup1 <> gGroup2 /
                if len(regions[0]) == 0:
                    g_dict = {
                        0: groups[0],
                        1: groups[1],
                        'region': [''],
                        'pick_mode': pick_mode}
                else:
                    g_dict = {
                        0: groups[0],
                        1: groups[1],
                        'region': regions,
                        'pick_mode': pick_mode}
    
            self.distin_g_list.append(g_dict)
            #log.info("{}".format(g_dict))
    
    
    def _get_members(self):
    # group_members =
    #     gHitomebore     : Hitomebore, Sasanishiki
    #     gArroz_da_terra : Arroz_da_terra, Kasalath
    #     gTakanari       : Mochibijin, Tawawakko
    #     gNortai         : Nortai, NERICA1
    
        group_members_str_org = self.ini['groups']['group_members']
        #log.info("{}".format(group_members_str_org))
        group_members_str = re.sub(r";", ",", group_members_str_org)
        #log.info("{}".format(group_members_str))
        g_members_tmp = group_members_str.split(',')

        group_line = ''
        for gm_tmp in g_members_tmp:
            if ':' in gm_tmp:
                group_line = "{};{},".format(group_line, gm_tmp)
            else:
                group_line = "{},{}".format(group_line, gm_tmp)

        group_line = re.sub(r"^;", "", group_line)
        group_line = re.sub(r",+", ",", group_line)
        #log.info("{}".format(group_line))

        g_members = group_line.split(';')
        #log.info("g_members {}".format(g_members))
    
        for g_member in g_members:
            gmem_list = g_member.split(':')
            #log.debug("{}".format(gmem_list))
    
            if len(gmem_list) == 1:
                #log.error("error, empty group_members {} exit.".format(
                #    g_member))
                sys.exit(1)
    
            elif gmem_list[0] == '' or gmem_list[1] == '':
                #log.error("error, empty group_members {} exit.".format(
                #    g_member))
                sys.exit(1)
    
            self.g_members_dict[gmem_list[0]] = gmem_list[1].split(',')
            #log.info("{}: {}".format(
            #    gmem_list[0],
            #    self.g_members_dict[gmem_list[0]]))

    
    def _rectify_variable(self):
    
        for section in self.ini.sections():
            for key in self.ini[section]:

                val = self.ini[section][key]
                # hash comment remove
                val = utl.strip_hash_comment(val)

                # remove \n at the beginning of value
                val = val.lstrip()

                # replace internal \n to semicolons
                val = val.replace('\n', ';')

                # replace white space to one space
                if key == 'group_members':
                    #log.info("{}".format(val))

                    val = re.sub(r"\s+", " ", val)
                    val = re.sub(r"\s*:\s*", ":", val)
                    val = re.sub(r" ", ",", val)
                    val = re.sub(r",+", ",", val)
                    val = re.sub(r",;", ";", val)
                    # g_members ['group0:TDr2946A_Mal,TDr1489A_Fem',
                    # 'TDr2284A_Mon,TDr1499A_Mon,TDr1509A_Fem',
                    # 'group1:TDr1510A_Fem,TDr3782A_Mal',
                    # 'TDr1533A_Mal,TDr1543A_Fem,TDr1858C_Fem']

                else:
                    val = re.sub(r"\s+", "", val)

                # reset
                self.ini[section][key] = val
    
                #log.info("{}:{}={}".format(section, key, val))
    
    
    
