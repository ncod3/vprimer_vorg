#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import errno
import time
import datetime

# global configuration
import vprimer.glv as glv

import subprocess as sbp
import pandas as pd
import pickle

import vprimer.utils as utl

from vprimer.blast import Blast
from vprimer.logging_config import LogConf


class RefFasta(object):

    def __init__(self):

        # as glv.ref.refseq
        self.refseq = dict()

        # as glv.ref.ref_fasta_pickle
        self.ref_fasta_pickle = ''

    def open_log(self):

        global log
        log = LogConf.open_log(__name__)


    def prepare_ref(self):

        # user's fasta: convert relative path to absolute path based on cwd
        if glv.conf.ref.startswith('/'):
            # originally absolute path
            glv.conf.ref_fasta_user = glv.conf.ref
        else:
            # cwd + relative path
            glv.conf.ref_fasta_user = "{}/{}".format(
                glv.conf.cwd, glv.conf.ref)

        log.info("glv.conf.ref_fasta_user {}".format(glv.conf.ref_fasta_user))

        # ref_fasta_user: existence confirmation
        if os.path.isfile(glv.conf.ref_fasta_user):
            log.info("{} found.".format(glv.conf.ref_fasta_user))
        else:
            log.info("{} not found. exit.".format(glv.conf.ref_fasta_user))
            sys.exit(1)
            
        # ext, basename, without_ext
        # https://note.nkmk.me/python-os-basename-dirname-split-splitext/
        basename_user = os.path.basename(glv.conf.ref_fasta_user)
        root_ext_pair = os.path.splitext(glv.conf.ref_fasta_user)
        without_ext = root_ext_pair[0]
        basename_without_ext = os.path.basename(without_ext)
        ext = root_ext_pair[1]

        # ref_fasta_slink_system
        # make symlink user's fasta to sys_ref_dir as .org(.gz)
        if ext == '.gz':
            glv.conf.ref_fasta_slink_system = "{}/{}{}".format(
                glv.conf.ref_dir, basename_user, '.org_slink.gz')
            # for blast
            glv.conf.blastdb_title = basename_without_ext
            
        else:
            glv.conf.ref_fasta_slink_system = "{}/{}{}".format(
                glv.conf.ref_dir, basename_user, '.org_slink')
            # for blast
            glv.conf.blastdb_title = basename_user

        glv.conf.blastdb = "{}/{}{}".format(
            glv.conf.ref_dir, glv.conf.blastdb_title, '.blastdb')


        log.info("glv.conf.blastdb={}".format(glv.conf.blastdb))

        if os.path.isfile(glv.conf.ref_fasta_slink_system):
            log.info("{} exist.".format(glv.conf.ref_fasta_slink_system))
        else:
            log.info("os.symlink {} {}.".format(
                glv.conf.ref_fasta_user, glv.conf.ref_fasta_slink_system))
            os.symlink(
                glv.conf.ref_fasta_user, glv.conf.ref_fasta_slink_system)

        log.info("ext ({}).".format(ext))

        # convert to bgz if ext is .gz and set to ref_fasta
        if ext == '.gz':

            # it should be convert to bgz in ref_dir
            glv.conf.ref_fasta = "{}/{}".format(
                glv.conf.ref_dir, basename_user)

            log.info("ext {}, glv.conf.ref_fasta={}.".format(
                ext, glv.conf.ref_fasta))

            # half of thread?
            cmd1 = 'bgzip -cd -@ {} {} \
                    | bgzip -@ {} > {}'.format(
                        glv.conf.parallele_full_thread,
                        glv.conf.ref_fasta_slink_system,
                        glv.conf.parallele_full_thread,
                        glv.conf.ref_fasta)

        else:
            # it should be convert to bgz in ref_dir
            glv.conf.ref_fasta = "{}/{}{}".format(
                glv.conf.ref_dir,
                basename_user,
                '.gz')

            cmd1 = 'bgzip -c -@ {} {} > {}'.format(
                        glv.conf.parallele_full_thread,
                        glv.conf.ref_fasta_slink_system,
                        glv.conf.ref_fasta)

        # execute
        if os.path.isfile(glv.conf.ref_fasta):
            log.info("{} exist.".format(glv.conf.ref_fasta))

        else:
            log.info("{} not exist. do cmd={}".format(
                glv.conf.ref_fasta,
                cmd1))

            utl.try_exec(cmd1)
 
        # make fai file
        cmd2 = 'samtools faidx {}'.format(
            glv.conf.ref_fasta, glv.conf.log_dir)

        glv.conf.ref_fasta_fai = "{}{}".format(glv.conf.ref_fasta, '.fai')

        if os.path.isfile(glv.conf.ref_fasta_fai):
            log.info("{} exist.".format(glv.conf.ref_fasta_fai))

        else:
            utl.try_exec(cmd2)

        # read fasta to dict vprimer.cnf.refseq
        glv.conf.ref_fasta_pickle = "{}{}".format(
            glv.conf.ref_fasta, '.pickle')
        self._read_fasta()

        # ref to makeblastdb
        Blast.makeblastdb()

        return self

    def pick_refseq(self, chrom, start_coordinate, end_coordinate):
        ''' for refseq substr, etc...
        '''

        slice_stt = start_coordinate - 1
        slice_end = end_coordinate

        #   1   2   3   4   5   6   coordinate   3-5 tho
        # +---+---+---+---+---+---+
        # | 0 | 1 | 2 | 3 | 4 | 5 | idx
        # | P | y | t | h | o | n |
        # +---+---+---+---+---+---+
        # 0   1   2   3   4   5   6 slice        2-5
        #

        return self.refseq[chrom][slice_stt:slice_end]


    def _read_fasta(self):

        if os.path.isfile(glv.conf.ref_fasta_pickle):
            log.info("{} exist.".format(glv.conf.ref_fasta_pickle))
            self._read_fasta_pickle()
        else:
            log.info("{} not exist.".format(glv.conf.ref_fasta_pickle))
            self._read_fasta_first()

    def _read_fasta_pickle(self):

        with open(glv.conf.ref_fasta_pickle, 'rb') as f:
            glv.ref.refseq = pickle.load(f)
        glv.conf.ref_fasta_chrom_list = glv.ref.refseq.keys()
        #log.info("glv.conf.ref_fasta_chrom_list={}.".format(
        #    glv.conf.ref_fasta_chrom_list))

    def _read_fasta_first(self):

# glv.conf.ref_fasta_fai
#chr01   43270923    7   60  61
#chr02   35937250    43992120    60  61
#chr03   36413819    80528332    60  61
#chr04   35502694    117549055   60  61
#chr05   29958434    153643468   60  61
#chr06   31248787    184101217   60  61
#chr07   29697621    215870825   60  61
#chr08   28443022    246063414   60  61
#chr09   23012720    274980494   60  61
#chr10   23207287    298376767   60  61
#chr11   29021106    321970850   60  61
#chr12   27531856    351475649   60  61

        # get chrom list from fai text
        df_fai = pd.read_csv(
            glv.conf.ref_fasta_fai, sep = '\t',
            header = None, index_col = None)

        # ref_fasta_chrom_list from fai column 0 (chrom name)
        glv.conf.ref_fasta_chrom_list = df_fai.loc[:, 0].to_list()

        log.info("fai {}, chrom cnt={}".format(
            glv.conf.ref_fasta_fai, len(glv.conf.ref_fasta_chrom_list)))

        # for each chrom name
        log.info("read refseq by samtools faidx from fasta {}".format(
            glv.conf.ref_fasta))
        start = time.time()

        chrom_seq_list = []
        last_chrom = ''
        for chrom in glv.conf.ref_fasta_chrom_list:

            # get sequence from samtools command
            cmd1 = "samtools faidx {} {}".format(glv.conf.ref_fasta, chrom)
            cmd_list = cmd1.split(' ')

            # log.info("{}".format(cmd_list))

            # using command output by pipe, get sequence into python
            proc = sbp.Popen(
                cmd_list, stdout = sbp.PIPE, stderr = sbp.PIPE)

            # got bytes (b'')
            for byte_line in proc.stdout:
                # bytes to str, strip \n
                b_line = byte_line.decode().strip()

                #print("{}={}(top)".format(chrom, len(chrom_seq_list)))
                # fasta header
                if b_line.startswith('>'):
                    # not the first time
                    if len(chrom_seq_list) != 0:
                        # dictionary
                        glv.ref.refseq[last_chrom] = ''.join(chrom_seq_list)
                        chrom_seq_list = []
                        continue

                else:
                    # append to list
                    chrom_seq_list.append(b_line)
                    #print("{}={}(append)".format(chrom, len(chrom_seq_list)))

            last_chrom = chrom

        if len(chrom_seq_list) != 0:
            glv.ref.refseq[last_chrom] = ''.join(chrom_seq_list)
            print("{},last_len={}".format(
                chrom, len(glv.ref.refseq[last_chrom])))

        elapsed_time = datetime.timedelta(seconds=(time.time() - start))
        log.info("read refseq done. time:{0}".format(elapsed_time))

        # pickle
        with open(glv.conf.ref_fasta_pickle, 'wb') as f:
            pickle.dump(glv.ref.refseq, f)

        log.info('dumped glv.ref.refseq->{}'.format(
            glv.conf.ref_fasta_pickle))

