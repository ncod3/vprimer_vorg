# -*- coding: utf-8 -*-

import sys
import os
import errno
import time

# global configuration
import vprimer.glv as glv
import vprimer.utils as utl

from vprimer.logging_config import LogConf
log = LogConf.open_log(__name__)

import pandas as pd

from vprimer.eval_variant import EvalVariant
from vprimer.product import Product
from vprimer.primer import PrimerInfo

class FormTxt(object):

    def __init__(self):

        self.line = ''

        self.marker_id = ''

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

        self.target_gno = -1
        self.target_len = 0
        self.enzyme_name = ''
        self.digest_pattern = ''

        self.g0_seq_target_len = ''
        self.g0_seq_target = ''
        self.g1_seq_target_len = ''
        self.g1_seq_target = ''
        self.seq_template_ref_len = ''
        self.seq_template_ref_abs_pos = ''
        self.seq_template_ref_rel_pos = ''
        self.PRIMER_PAIR_0_PRODUCT_SIZE = ''
        self.PRIMER_LEFT_0 = ''
        self.left_primer_id = ''
        self.PRIMER_LEFT_0_SEQUENCE = ''
        self.PRIMER_RIGHT_0 = ''
        self.right_primer_id = ''
        self.PRIMER_RIGHT_0_SEQUENCE = ''
        self.SEQUENCE_TEMPLATE = ''


    def format_text(self):

        # progress check
        if utl.progress_check('formsafe') == False and \
            utl.progress_check('formfail') == False:

            log.info("progress={} so skip form.".format(
                glv.conf.progress))
            return
        log.info("Start processing {}".format('formsafe'))

        # for each distinguish_groups
        for distin_dict in glv.outlist.distin_files:

            # read variant file
            primer_file = distin_dict['primer']['out_path']
            df_distin = pd.read_csv(
                primer_file, sep='\t', header=0, index_col=None)

            # complete == 1 or == 0
            safe = 1
            fail = 0
            for complete, proc in zip([fail, safe], ['formfail', 'formsafe']):
                log.info("{} {}".format(complete, proc))

                df_distin_complete = \
                    df_distin[df_distin['complete'] == complete]

                #------------------------
                # check chrom-pos duplicate marker 
                df_chrom_pos = df_distin_complete.loc[:, ['chrom', 'pos']]
                df_chrom_pos_duplicated = \
                    df_chrom_pos[df_chrom_pos.duplicated()]

                duplicate_pos_dict = dict()
                for c_p_row in df_chrom_pos_duplicated.itertuples():

                    chrom = c_p_row[1]
                    pos = c_p_row[2]

                    if not chrom in duplicate_pos_dict:
                        duplicate_pos_dict[chrom] = dict();

                    if not pos in duplicate_pos_dict[chrom]:
                        duplicate_pos_dict[chrom][pos] = pos

                #------------------------
                # file name to write out result to text
                out_txt_file = distin_dict[proc]['out_path']
                log.info("out_txt_file={}.".format(out_txt_file))

                utl.save_to_tmpfile(out_txt_file)

                with open(out_txt_file, mode='a') as f:

                    # write header
                    f.write("{}\n".format(
                        distin_dict['formsafe']['hdr_text']))

                    # each variant
                    for primer_df_row in df_distin_complete.itertuples():

                        self._prepare_from_primer_file(
                            primer_df_row, distin_dict)

                        self._format_product(duplicate_pos_dict)

                        # 書き出す
                        f.write("{}\n".format(self.line))


    def _get_group_product_size(self):

        # rel template info
        rel_frag_pad_pre_stt, \
        rel_frag_pad_pre_end, \
        rel_around_seq_pre_stt, \
        rel_around_seq_pre_end, \
        rel_pos, \
        rel_around_seq_aft_stt, \
        rel_around_seq_aft_end, \
        rel_frag_pad_aft_stt, \
        rel_frag_pad_aft_end = \
            Product.separate_seq_template_pos(
                self.seq_template_ref_rel_pos)

        # variant sequenceの、REFとの差分
        vseq_lens_ano = [int(x) for x in self.vseq_lens_ano_str.split(',')]
        vseq_ref_len = vseq_lens_ano[0]

        # gnoに対応するanoの対応表
        ano_gno = [self.g0_ano, self.g1_ano]

        # primer3からのprimer positionと長さ情報
        # 左は、productのstart
        product_stt_pos, left_len = map(int, self.PRIMER_LEFT_0.split(','))

        # rithr_stt_posは、revcomのprimerの開始なので、
        # product sizeの末端と考えて良い
        product_end_pos, right_len = map(int, self.PRIMER_RIGHT_0.split(','))

        # group0, group1のproduct sizeの計算
        l_product_size = list()
        l_product_end_pos = list()

        for gno in range(2):
            # グループのvseqと、REFの差分を計算する
            vseq_gno_len = vseq_lens_ano[ano_gno[gno]]
            diff_len = vseq_ref_len - vseq_gno_len

            # productの末端の調整
            my_product_end_pos = product_end_pos - diff_len
            my_product_size = my_product_end_pos - product_stt_pos + 1

            l_product_end_pos.append(my_product_end_pos)
            l_product_size.append(my_product_size)

        # capsのdigest posの位置とsizeの計算
        enzyme = self.enzyme_name
        digested_size = ['-', '-']
        # 切られる側
        d_gno = self.target_gno
        # if glv.INDEL, this is indel length diff. !=INDEL, it is digest pos
        digest_diff = self.target_len

        if self.var_type == glv.INDEL:

            # 長い側にサイズをいれる
            digested_size[d_gno] = l_product_size[d_gno]

        else:

            enzyme = self.marker_info.split(',')[0]
            # vseqとaround_seqの開始基準位置(rel_frag_pad_pre_end)
            target_stt_pos = rel_frag_pad_pre_end

            # productの最終ポジション
            my_product_end_pos = l_product_end_pos[d_gno]

            # digestされたときの、右側の開始ポジション
            digest_right_stt_pos = target_stt_pos + self.target_len - 1

            # 右側のサイズ
            d_size_R = my_product_end_pos - digest_right_stt_pos + 1
            # 左側のサイズ
            d_size_L = l_product_size[d_gno] - d_size_R

            digested_size[d_gno] = "{}/{}".format(d_size_L, d_size_R)

            digest_diff = abs(d_size_R - d_size_L)

            # <-----><=====R=====><----->
            # <-----><=====Raaaa=====><----->
            #       012345678901234567890
            #            /5   /0     /7
            #            /digest_right_stt_pos
            # <.........><..................>
            #            5     10-5+1      10
            #            a                  b
            #        

        return \
            l_product_size[0], \
            l_product_size[1], \
            enzyme, \
            digested_size[0], \
            digested_size[1], \
            digest_diff


    def _add_comment(self, duplicate_pos_dict):

        comment_list = list()

        #log.debug("{}".format(duplicate_pos_dict))

        if len(duplicate_pos_dict) == 0:
            pass
        elif self.pos in duplicate_pos_dict[self.chrom]:
            comment_list = [glv.COMMENT_dup]

        if len(comment_list) != 0:
            comment_list += [self.set_enz_cnt]

#        if self.set_n != 1:
#            if self.set_n == 2:
#                comment = glv.COMMENT_AABC + ','
#            if self.set_n == 4:
#                comment = glv.COMMENT_ABCD + ','

#        if comment == '':
#            comment = glv.COMMENT_nop

        if len(comment_list) == 0:
            comment_list = '-'

        return ','.join(comment_list)

    def _format_product(self, duplicate_pos_dict):

        product_size_0, product_size_1, enzyme, \
        digested_size_0, digested_size_1, digest_diff = \
            self._get_group_product_size()

        comment = self._add_comment(duplicate_pos_dict)

        #--------------------
        line_list = list()

        line_list += [self.chrom]
        line_list += [self.pos]
        line_list += [self.try_cnt]
        line_list += [self.complete]

        line_list += [comment]

        line_list += [self.var_type]
        line_list += [self.gts_segr_lens]
        line_list += [enzyme]

        line_list += [self.g0_name]
        line_list += [product_size_0]

        line_list += [digested_size_0]

        line_list += [digest_diff]

        line_list += [digested_size_1]

        line_list += [product_size_1]
        line_list += [self.g1_name]

        line_list += [self.left_primer_id]
        line_list += [self.PRIMER_LEFT_0_SEQUENCE]

        line_list += [self.right_primer_id]
        line_list += [self.PRIMER_RIGHT_0_SEQUENCE]

        self.line = '\t'.join(map(str, line_list))


    def _prepare_from_primer_file(self, primer_df_row, distin_dict):

        # ほとんど同じなのはなんとかしないと

        hdr_dict = distin_dict['primer']['hdr_dict']

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
            PrimerInfo.get_basic_primer_info(primer_df_row, hdr_dict)

        #log.debug("self.chrom={} pos={}".format(self.chrom, self.pos))


        self.marker_id = str(primer_df_row[hdr_dict['marker_id']])

        self.try_cnt = str(primer_df_row[hdr_dict['try_cnt']])
        self.complete = str(primer_df_row[hdr_dict['complete']])
        self.blast_check = str(primer_df_row[hdr_dict['blast_check']])

        self.g0_seq_target_len = \
            int(primer_df_row[hdr_dict['g0_seq_target_len']])
        self.g0_seq_target = \
            str(primer_df_row[hdr_dict['g0_seq_target']])
        self.g1_seq_target_len = \
            int(primer_df_row[hdr_dict['g1_seq_target_len']])
        self.g1_seq_target = \
            str(primer_df_row[hdr_dict['g1_seq_target']])

        self.seq_template_ref_len = \
            int(primer_df_row[hdr_dict['seq_template_ref_len']])
        self.seq_template_ref_abs_pos = \
            str(primer_df_row[hdr_dict['seq_template_ref_abs_pos']])
        self.seq_template_ref_rel_pos = \
            str(primer_df_row[hdr_dict['seq_template_ref_rel_pos']])

        self.PRIMER_PAIR_0_PRODUCT_SIZE = \
            int(primer_df_row[hdr_dict['PRIMER_PAIR_0_PRODUCT_SIZE']])
        self.PRIMER_LEFT_0 = \
            str(primer_df_row[hdr_dict['PRIMER_LEFT_0']])
        self.left_primer_id = \
            str(primer_df_row[hdr_dict['left_primer_id']])
        self.PRIMER_LEFT_0_SEQUENCE = \
            str(primer_df_row[hdr_dict['PRIMER_LEFT_0_SEQUENCE']])
        self.PRIMER_RIGHT_0 = \
            str(primer_df_row[hdr_dict['PRIMER_RIGHT_0']])
        self.right_primer_id = \
            str(primer_df_row[hdr_dict['right_primer_id']])
        self.PRIMER_RIGHT_0_SEQUENCE = \
            str(primer_df_row[hdr_dict['PRIMER_RIGHT_0_SEQUENCE']])
        self.SEQUENCE_TEMPLATE = \
            str(primer_df_row[hdr_dict['SEQUENCE_TEMPLATE']])


