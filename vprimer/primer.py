# -*- coding: utf-8 -*-

import sys
import os
import errno
import time
import re

# global configuration
import vprimer.glv as glv
import vprimer.utils as utl

from vprimer.logging_config import LogConf
log = LogConf.open_log(__name__)

import pandas as pd
import vcfpy
import subprocess as sbp

from joblib import Parallel, delayed
#import dill as pickle

from vprimer.product import Product
from vprimer.allele_select import AlleleSelect
from vprimer.eval_variant import EvalVariant
from vprimer.variant import Variant
from vprimer.inout_primer3 import InoutPrimer3
from vprimer.blast import Blast

class Primer(object):

    def __init__(self):

        pass

    def construct_primer(self):

        # progress check
        if utl.progress_check('primer') == False:
            log.info("progress={} so skip primer.".format(
                glv.conf.progress))
            return
        log.info("Start processing {}".format('primer'))


        # for each distinguish_groups
        for distin_dict in glv.outlist.distin_files:

            marker_file = distin_dict['marker']['out_path']
            df_distin = pd.read_csv(
                marker_file, sep='\t', header=0, index_col=None)

            out_txt_file = distin_dict['primer']['out_path']
            utl.save_to_tmpfile(out_txt_file)

            with open(out_txt_file, mode='a') as f:
                # write header
                #f.write("{}\n".format(distin_dict['primer']['hdr_text']))

                start = time.time()

                if glv.conf.parallel == True:
                    log.info(
                        "do Parallel cpu {}, parallele {} blast {}".format(
                            glv.conf.thread,
                            glv.conf.parallel_blast_cnt,
                            glv.conf.blast_num_threads))

                    Parallel(
                        n_jobs=glv.conf.parallel_blast_cnt,
                        backend="threading")(
                        [
                            delayed(self._loop_primer3_check_blast) \
                                (distin_dict, marker_df_row, f) \
                                for marker_df_row in df_distin.itertuples()
                        ]
                    )

                else:
                    log.info("do Serial cpu {} / serial {} blast {}".format(
                        glv.conf.thread,
                        1,
                        glv.conf.blast_num_threads))

                    for marker_df_row in df_distin.itertuples():

                        self._loop_primer3_check_blast(
                            distin_dict, marker_df_row, f)

            utl.sort_file(
                'primer', distin_dict, out_txt_file,
                'chrom', 'pos', 'try_cnt', 'number')

            log.info("primer {} {}".format(
                utl.elapsed_time(time.time(), start),
                distin_dict['primer']['base_nam']))


    def _loop_primer3_check_blast(self, distin_dict, marker_df_row, f):

        prinfo = PrimerInfo()
        prinfo.prepare_from_marker_file(distin_dict, marker_df_row)
        prinfo.get_excluded_region()
        prinfo.make_p3_info()
        blast_check_result_list = list()

        line = ''
        try_cnt = 0
        blast_check = ''

        while True:

            #log.debug("<try_cnt {}> SEQUENCE_EXCLUDED_REGION={}".format(
            #    try_cnt,
            #    iopr3.get_sequence_excluded_region()))
            try_cnt += 1

            # primer3
            self._do_primer3_pipe(prinfo.iopr3)

            #log.debug("")
            #log.debug("{}".format(prinfo.iopr3.get_primer3_out()))

            #log.debug("\nP3_COMMENT={}".format(iopr3.get_p3_comment()))

            # break if cann't generate primer
            if prinfo.iopr3.PRIMER_ERROR != '':
                #log.info("try_cnt {} {} PRIMER_ERROR={}".format(
                #    try_cnt,
                #    prinfo.iopr3.get_sequence_id(),
                #    prinfo.iopr3.PRIMER_ERROR))
                break

            elif prinfo.iopr3.PRIMER_PAIR_NUM_RETURNED == 0:
                #log.info("try_cnt {} {} PRIMER_PAIR_NUM_RETURNED={}".format(
                #    try_cnt,
                #    prinfo.iopr3.get_sequence_id(),
                #    prinfo.iopr3.PRIMER_PAIR_NUM_RETURNED))
                break

            # calc abs pos
            abs_left_stt, abs_left_end, \
            abs_right_stt, abs_right_end = \
                self._get_primer_pos_info(prinfo)

            # make fasta id
            left_fasta_id = self._make_primer_name(
                prinfo.chrom, abs_left_stt, abs_left_end, "plus")
            right_fasta_id = self._make_primer_name(
                prinfo.chrom, abs_right_stt, abs_right_end, "minus")

            prinfo.iopr3.set_primer_name(left_fasta_id, right_fasta_id)

            # search my blastn-short
            blast_check_result_list = \
                Blast.primer_blast_check(
                    left_fasta_id, right_fasta_id,
                    prinfo.iopr3.get_primer_left_seq(),
                    prinfo.iopr3.get_primer_right_seq())

            # break if complete
            if len(blast_check_result_list) == 0:

                #log.info("complete try_cnt {} {}".format(
                #    try_cnt, prinfo.iopr3.get_sequence_id()))

                #------------------------------
                complete = 1
                line, blast_check = self._primer_complete_to_line(
                    complete, blast_check_result_list, prinfo, try_cnt)
                f.write('{}\n'.format(line))
                break

            else:
                # go to next chance to add the primer pos to ex
                prinfo.iopr3.add_ex_region(prinfo.iopr3.get_primer_region())

                #log.info("next try_cnt {} {}".format(
                #    try_cnt,
                #    prinfo.iopr3.get_sequence_excluded_region()))

                #------------------------------
                complete = 0
                line, blast_check = self._primer_complete_to_line(
                    complete, blast_check_result_list, prinfo, try_cnt)
                f.write('{}\n'.format(line))


        if len(blast_check_result_list) != 0:
            log.info("skipped try_cnt {} {} {}".format(
                try_cnt,
                prinfo.iopr3.get_sequence_id(),
                blast_check))
           # log.info(iopr3.p3_out)

#        if line != '':
#            f.write('{}\n'.format(line))
#            # if you need
#            f.flush()


    def _primer_complete_to_line(
        self, complete, blast_check_result_list, prinfo, try_cnt):

        # blast_check
        blast_check = "-"
        add_cnt = len(blast_check_result_list) - 1
        if add_cnt == -1:
            pass
        elif add_cnt == 0:
            blast_check = "{}".format(blast_check_result_list[0])
        else:
            blast_check = "{}(+{})".format(
                blast_check_result_list[0], add_cnt)


        l_list = list()

        # to primer out file
        l_list += [prinfo.marker_id]

        l_list += [prinfo.chrom]
        l_list += [prinfo.pos]
        l_list += [prinfo.targ_grp]
        l_list += [prinfo.gts_segr_lens]
        l_list += [prinfo.targ_ano]
        l_list += [prinfo.set_enz_cnt]
        l_list += [prinfo.var_type]
        l_list += [prinfo.marker_info]
        l_list += [prinfo.vseq_lens_ano_str]

        l_list += [try_cnt]
        l_list += [complete]
        l_list += [blast_check]

        l_list += [prinfo.target_gno]
        l_list += [prinfo.target_len]

        l_list += [prinfo.g0_seq_target_len]
        l_list += [prinfo.g0_seq_target]
        l_list += [prinfo.g1_seq_target_len]
        l_list += [prinfo.g1_seq_target]

        l_list += [prinfo.seq_template_ref_len]
        l_list += [prinfo.seq_template_ref_abs_pos]
        l_list += [prinfo.seq_template_ref_rel_pos]

        l_list += [prinfo.iopr3.get_primer_product_size()]

        l_list += [prinfo.iopr3.get_primer_left()]
        l_list += [prinfo.iopr3.get_primer_left_id()]
        l_list += [prinfo.iopr3.get_primer_left_seq()]

        l_list += [prinfo.iopr3.get_primer_right()]
        l_list += [prinfo.iopr3.get_primer_right_id()]
        l_list += [prinfo.iopr3.get_primer_right_seq()]

        l_list += [prinfo.SEQUENCE_TARGET]
        # 
        l_list += [prinfo.iopr3.get_sequence_excluded_region()]
        l_list += [prinfo.seq_template_ref]

        return '\t'.join(map(str, l_list)), blast_check


    def _make_primer_name(
            self, chrom, abs_primer_stt_pos, abs_primer_end_pos, strand):

        # {NC_028450.1}44676.44700.plus
        #primer_name = "{{{}}}{}.{}.{}".format(
        #    chrom,
        #    abs_primer_stt_pos,
        #    abs_primer_end_pos,
        #    strand)

        # NC_028450.1:44676-44700:plus
        primer_name = "{}:{}-{}:{}".format(
            chrom,
            abs_primer_stt_pos,
            abs_primer_end_pos,
            strand)

        return primer_name


    def _do_primer3_pipe(self, iopr3):

        # exec primer3 through pipe
        primer3_in = iopr3.get_p3_input()

        primer3_out_p = sbp.Popen(
            ['primer3_core'],
            stdin=sbp.PIPE,
            stdout=sbp.PIPE)

        primer3_out = primer3_out_p.communicate(
            primer3_in.encode())[0].decode()

        iopr3.set_primer3_out(primer3_out)


    def _get_primer_pos_info(self, prinfo):

        rel_left_stt, left_len = prinfo.iopr3.get_primer_left_info()

        #log.debug("rel_left_stt={}".format(rel_left_stt))
        #log.debug("self.abs_frag_pad_pre_stt={}".format(
        #    self.abs_frag_pad_pre_stt))

        rel_right_end, right_len = prinfo.iopr3.get_primer_right_info()

        rel_right_stt = rel_right_end - right_len + 1

        abs_left_stt = self._get_abspos(
            rel_left_stt, prinfo.abs_frag_pad_pre_stt)


        abs_left_end = abs_left_stt + left_len - 1

        abs_right_stt = self._get_abspos(
            rel_right_stt, prinfo.abs_frag_pad_pre_stt)

        abs_right_end = abs_right_stt + right_len - 1

        return \
            abs_left_stt, abs_left_end, \
            abs_right_stt, abs_right_end


    def _get_abspos(self, ref_pos, abs_template_stt):

        # 11
        #  1234567890
        #        7 11+7=18 -1

        return abs_template_stt + ref_pos - 1


class PrimerInfo(object):

    def __init__(self):

        self.chrom = ''
        self.pos = 0

        self.targ_grp = ''
        self.g0_name = ''
        self.g1_name = ''

        self.gts_segr_lens = ''

        self.targ_ano = ''
        self.g0_ano = -1
        self.g1_ano = -1

        self.set_enz_cnt = ''
        self.var_type = ''
        self.marker_info = ''

        self.vseq_lens_ano_str = ''

        self.marker_id = ''

        self.target_gno = -1
        self.target_len = 0
        self.enzyme_name = ''
        self.digest_pattern = ''


        self.g0_seq_target_len = 0
        self.g0_seq_target = ''
        self.g1_seq_target_len = 0
        self.g1_seq_target = ''
        self.seq_template_ref_len = 0
        self.seq_template_ref_abs_pos = ''
        self.seq_template_ref_rel_pos = ''
        self.SEQUENCE_TARGET = ''
        self.seq_template_ref = ''

        # abs template info
        self.abs_frag_pad_pre_stt = 0
        self.abs_frag_pad_pre_end = 0
        self.abs_around_seq_pre_stt = 0
        self.abs_around_seq_pre_end = 0
        self.abs_pos = 0
        self.abs_around_seq_aft_stt = 0
        self.abs_around_seq_aft_end = 0
        self.abs_frag_pad_aft_stt = 0
        self.abs_frag_pad_aft_end = 0

        # abs template info
        self.rel_frag_pad_pre_stt = 0
        self.rel_frag_pad_pre_end = 0
        self.rel_around_seq_pre_stt = 0
        self.rel_around_seq_pre_end = 0
        self.rel_pos = 0
        self.rel_around_seq_aft_stt = 0
        self.rel_around_seq_aft_end = 0
        self.rel_frag_pad_aft_stt = 0
        self.rel_frag_pad_aft_end = 0

        self.SEQUENCE_EXCLUDED_REGION = ''
   
        self.p3_comment = ''


    def prepare_from_marker_file(self, distin_dict, marker_df_row):

        hdr_dict = distin_dict['marker']['hdr_dict']

        # basic
        self.chrom, self.pos, \
        self.targ_grp, self.g0_name, self.g1_name, \
        self.gts_segr_lens, \
        self.targ_ano, self.g0_ano, self.g1_ano, \
        self.set_enz_cnt, \
        self.var_type, self.marker_info, \
        self.vseq_lens_ano_str, \
        self.enzyme_name, self.target_gno, self.target_len, \
        self.digest_pattern, \
        self.target_gno, self.target_len = \
            PrimerInfo.get_basic_primer_info(marker_df_row, hdr_dict)

        self.marker_id = self._get_marker_id()

        self.g0_seq_target_len = \
            int(marker_df_row[hdr_dict['g0_seq_target_len']])
        self.g0_seq_target = \
            str(marker_df_row[hdr_dict['g0_seq_target']])
        self.g1_seq_target_len = \
            int(marker_df_row[hdr_dict['g1_seq_target_len']])
        self.g1_seq_target = \
            str(marker_df_row[hdr_dict['g1_seq_target']])

        self.seq_template_ref_len = \
            int(marker_df_row[hdr_dict['seq_template_ref_len']])
        self.seq_template_ref_abs_pos = \
            str(marker_df_row[hdr_dict['seq_template_ref_abs_pos']])

        #log.debug("{}".format(self.seq_template_ref_abs_pos))

        self.seq_template_ref_rel_pos = \
            str(marker_df_row[hdr_dict['seq_template_ref_rel_pos']])

        #log.debug("{}".format(self.seq_template_ref_rel_pos))

        self.SEQUENCE_TARGET = \
            str(marker_df_row[hdr_dict['SEQUENCE_TARGET']])
        self.seq_template_ref = \
            str(marker_df_row[hdr_dict['seq_template_ref']])

        # abs template info
        self.abs_frag_pad_pre_stt, \
        self.abs_frag_pad_pre_end, \
        self.abs_around_seq_pre_stt, \
        self.abs_around_seq_pre_end, \
        self.abs_pos, \
        self.abs_around_seq_aft_stt, \
        self.abs_around_seq_aft_end, \
        self.abs_frag_pad_aft_stt, \
        self.abs_frag_pad_aft_end = \
            Product.separate_seq_template_pos(
                self.seq_template_ref_abs_pos)

        # rel template info
        self.rel_frag_pad_pre_stt, \
        self.rel_frag_pad_pre_end, \
        self.rel_around_seq_pre_stt, \
        self.rel_around_seq_pre_end, \
        self.rel_pos, \
        self.rel_around_seq_aft_stt, \
        self.rel_around_seq_aft_end, \
        self.rel_frag_pad_aft_stt, \
        self.rel_frag_pad_aft_end = \
            Product.separate_seq_template_pos(
                self.seq_template_ref_rel_pos)


    def _get_relpos(self, abs_pos):

        #   | self.abs_frag_pad_pre_stt
        #                 self.abs_frag_pad_aft_end|
        #   <-------><========>P<=========><------->
        #            |self.abs_around_seq_pre_stt
        #      self.abs_around_seq_aft_end|
        #   |101 109|
        #   123456789
        #       |105
        #   105-101+1 = 5

        rel_pos = abs_pos - self.abs_frag_pad_pre_stt + 1

        return rel_pos


    def get_excluded_region(self):

        SEQUENCE_EXCLUDED_REGION = list()

        #logf_l = ["{} self.pos={} rel_pos={} rel_end_pos={} "]
        #logf_l += ["region_len={} template_len={}"]
        #logf = "".join(logf_l)

        region = "{}:{}-{}".format(
            self.chrom,
            self.abs_frag_pad_pre_stt,
            self.abs_frag_pad_aft_end)

        reader = vcfpy.Reader.from_path(glv.conf.vcf_file)
        vcf_ittr = reader.fetch(region)

        # access to vcf using iterater
        for record in vcf_ittr:
            sample0 = glv.conf.g_members_dict[self.g0_name][0]
            sample1 = glv.conf.g_members_dict[self.g1_name][0]

            # もし、サンプル間でvariantが見つかった場合は、
            s0_0, s0_1, s1_0, s1_1 = \
                AlleleSelect.record_call_for_sample(record, sample0, sample1)

            if Variant.is_same_gt(s0_0, s0_1, s1_0, s1_1) == False:

                # 20200713 here
                if self.pos != record.POS:
                    rel_pos = self._get_relpos(record.POS) #
                    # そのポジションのrefのvseq分を登録する
                    # 長さがfragment長を超える場合は、
                    # 調整する。
                    # 見つかったのはPOS
                    # REFのlength

                    # self.pos|
                    #     1036|
                    # ATGCATGCA ref_len=1
                    #         T
                    #         C
                    #         1036 + 1 - 1
                    region_len = len(record.REF)
                    rel_end_pos = rel_pos + region_len - 1

                    #  pos   len     end
                    #  1036 (10)     1045
                    #            1041 temp_len

#                    log.debug(logf.format(
#                        1, self.pos, rel_pos, rel_end_pos, region_len,
#                        self.seq_template_ref_len))


                    #log.debug("{}, {}".format(
                    #    rel_end_pos, self.seq_template_ref_len))
                    #log.debug("{}, {}".format(
                    #    type(rel_end_pos), type(self.seq_template_ref_len)))


                    if rel_end_pos > self.seq_template_ref_len:
                        diff_len = rel_end_pos - self.seq_template_ref_len
                        region_len = region_len - diff_len

#                        log.debug(logf.format(
#                            2, self.pos, rel_pos, rel_end_pos, region_len,
#                            self.seq_template_ref_len))

                    SEQUENCE_EXCLUDED_REGION += [
                        "{},{}".format(rel_pos, region_len)]

        self.SEQUENCE_EXCLUDED_REGION = " ".join(SEQUENCE_EXCLUDED_REGION)
        #log.debug("SEQUENCE_EXCLUDED_REGION={}".format(
        #    self.SEQUENCE_EXCLUDED_REGION))


    def make_p3_info(self):

        # save information
        self.p3_comment = self._make_p3_comment()
        #log.debug("{}".format(self.p3_comment))

        self.iopr3 = InoutPrimer3()
        self.iopr3.set_p3_comment(self.p3_comment)
        self.iopr3.set_sequence_target(self.SEQUENCE_TARGET)
        self.iopr3.add_ex_region(self.SEQUENCE_EXCLUDED_REGION)
        self.iopr3.set_sequence_id(self.marker_id)
        self.iopr3.set_sequence_template(self.seq_template_ref)

        #log.debug("")
        #log.debug("{}".format(self.iopr3.get_p3_input()))


    def _make_p3_comment(self):

        p3_comment = list()
        p3_comment += ["{}:{}".format('marker_id', self.marker_id)]
        p3_comment += ["{}:{}".format('var_type', self.var_type)]
        p3_comment += ["{}:{}".format('g0_name', self.g0_name)]
        p3_comment += ["{}:{}".format('g1_name', self.g1_name)]
        p3_comment += ["{}:{}".format('marker_info', self.marker_info)]


        p3_comment += ["{}:{}".format(
            'g0_seq_target_len', self.g0_seq_target_len)]
        p3_comment += ["{}:{}".format(
            'g0_seq_target', self.g0_seq_target)]
        p3_comment += ["{}:{}".format(
            'g1_seq_target_len', self.g1_seq_target_len)]
        p3_comment += ["{}:{}".format(
            'g1_seq_target', self.g1_seq_target)]

        p3_comment += ["{}:{}".format(
            'seq_template_ref_len', self.seq_template_ref_len)]
        p3_comment += ["{}:{}".format(
            'seq_template_ref_abs_pos', self.seq_template_ref_abs_pos)]
        p3_comment += ["{}:{}".format(
            'seq_template_ref_rel_pos', self.seq_template_ref_rel_pos)]

        return ';'.join(map(str, p3_comment))


    def _get_marker_id(self):

        if self.var_type == glv.INDEL:
            enzyme_len = "{}.{}".format(
                self.target_gno,
                self.target_len)
        else:
            enzyme_len = "{}.{}.{}".format(
                self.enzyme_name,
                self.target_gno,
                self.target_len)

        marker_id = "{}.{}.{}.{}.{}".format(
            self.chrom,
            self.pos,
            self.targ_ano,
            self.var_type,
            enzyme_len)

        return marker_id

    # 自分自身のクラスメソッドは、PrimerInfoで呼ぶ
    @classmethod
    def get_basic_primer_info(cls, df_row, hdr_dict):

        chrom = str(df_row[hdr_dict['chrom']])
        pos = int(df_row[hdr_dict['pos']])

        targ_grp = str(df_row[hdr_dict['targ_grp']])
        g0_name, g1_name = targ_grp.split(',')

        gts_segr_lens = df_row[hdr_dict['gts_segr_lens']]

        targ_ano = str(df_row[hdr_dict['targ_ano']])
        g0_ano, g1_ano = map(int, targ_ano.split(','))

        set_enz_cnt = str(df_row[hdr_dict['set_enz_cnt']])
        var_type = str(df_row[hdr_dict['var_type']])
        marker_info = str(df_row[hdr_dict['marker_info']])

        vseq_lens_ano_str = str(df_row[hdr_dict['vseq_lens_ano_str']])

        enzyme_name = '-'
        digest_pattern = '-'

        if var_type == glv.INDEL:
            target_gno, target_len = map(int, marker_info.split(','))
        else:
            enzyme_name, \
            target_gno, target_len, digest_pattern = marker_info.split(',')
            target_gno = int(target_gno)
            target_len = int(target_len)

        return \
            chrom, pos, \
            targ_grp, g0_name, g1_name, \
            gts_segr_lens, \
            targ_ano, g0_ano, g1_ano, \
            set_enz_cnt, \
            var_type, marker_info, \
            vseq_lens_ano_str, \
            enzyme_name, target_gno, target_len, \
            digest_pattern, \
            target_gno, target_len

