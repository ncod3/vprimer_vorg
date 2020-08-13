# -*- coding: utf-8 -*-

import sys
import os
import errno

import logging
log = logging.getLogger(__name__)

# global configuration
import vprimer.glv as glv
import vprimer.utils as utl

from Bio import Restriction
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA

from vprimer.product import Product

class EvalVariant(object):

    def __init__(self):

        self.line = ''

        self.chrom = ''
        self.pos = 0

        self.targ_grp = ''
        self.g0_name = ''
        self.g1_name = ''
        self.gname_gno = list()

        self.gts_segr_lens = ''
        self.set_n_str = ''
        self.set_my_no = 0
        self.set_n = 0

        self.targ_ano = ''
        self.g0_ano = -1
        self.g1_ano = -1
        self.ano_gno = list()

        self.len_g0g1_dif_long = ''
        self.g0_len = 0
        self.g1_len = 0
        self.diff_len = 0
        self.long_gno = -1

        self.short_gno = -1

        self.var_type = ''
        self.vseq_ano_str = ''
        self.vseq_ano = list()
        self.vseq_lens_ano = list()

        # --------------------------------
        self.marker_available = False

        self.vseq_ref = ''
        self.vseq_ref_len = ''

        # around_seq
        self.abs_around_seq_pre_stt = 0
        self.abs_around_seq_pre_end = 0
        self.abs_around_seq_aft_stt = 0
        self.abs_around_seq_aft_end = 0
    
        self.around_seq_pre = ''
        self.around_seq_aft = ''

        self.around_seq_pre_len = 0
        self.around_seq_aft_len = 0

        # target ref
        self.seq_target_ref = ''
        self.seq_target_ref_len = 0
        self.seq_target_ref_rel_pos = ''

        self.SEQUENCE_TARGET = ''

        # product
        self.g0_prod = Product()
        self.g1_prod = Product()
        self.gr_prod = [self.g0_prod, self.g1_prod]

        # caps
        self.caps_dict = dict()
        self.caps_found = 0

        # digest_gno:g0_prod_size:g1_prod_size:128/29/20
        #self.digest_info = ''

        # seq_template
        self.fragment_pad_len = 0
        
        self.abs_frag_pad_pre_stt = 0
        self.abs_frag_pad_pre_end = 0
        self.abs_frag_pad_aft_stt = 0
        self.abs_frag_pad_aft_end = 0

        self.seq_template_ref_abs_pos = ''
        self.seq_template_ref_rel_pos = ''

        self.frag_pad_pre = ''
        self.frag_pad_aft = ''
        self.seq_template_ref = ''

        self.frag_pad_pre_len = 0
        self.frag_pad_aft_len = 0
        self.seq_template_ref_len = 0


    def evaluate_for_marker(
        self, variant_df_row, distin_dict, enzyme_list):
        ''' 1 variant 1 exec
        '''


        # variantからデータ収集
        self._prepare_from_variant(variant_df_row, distin_dict)

        # glv.MODE_ALL    = "all"
        # glv.MODE_INDEL  = "indel"    # indel
        # glv.MODE_SNPMNV = "snpmnv"   # snp mnv mind
        # glv.MODE_SNP    = "snp"      # snp

        # skip if pick_mode is different
        if utl.is_my_pick_mode(
            self.var_type, distin_dict['pick_mode']) != True:
            return

        # around_seqを作成
        self._pick_ref_around_seq()

        # g0, g1のproduct作成
        for gno in range(2):
            ano = self.ano_gno[gno]
            self.gr_prod[gno].set_info(
                gno,
                ano,
                self.gname_gno[gno],
                self.chrom,
                self.pos,
                self.var_type,
                self.vseq_ano[ano],
                self.vseq_lens_ano[ano],
                self.around_seq_pre, 
                self.around_seq_aft,
            )

        # typeがindel以外ならcaps調査
        if self.var_type != glv.INDEL:
            self.check_caps(enzyme_list)
        else:
            self.marker_available = True

        #log.debug("{} {} {} {} {}".format(
        #    self.variant_id,
        #    self.marker_available,
        #    self.long_gno,
        #    self.diff_len,
        #    self.caps_dict
        #    ))


    def check_caps(self, enzyme_list):

        caps_result = list()

        for gno in range(2):
            caps_result_dict, caps_result_dict_str = \
                EvalVariant.get_caps_result(
                    self.gr_prod[gno].seq_target, enzyme_list)

            caps_result.append(caps_result_dict)

        # type(BsmFI) RestrictionType
        for enz_res_type in caps_result[0].keys():

            # compare  {EcoRI: []},  {EcoRI: [17]}
            len0 = len(caps_result[0][enz_res_type])
            len1 = len(caps_result[1][enz_res_type])

            if len0 == 1 and len1 == 0 or len0 == 0 and len1 == 1:
                enzyme_name = str(enz_res_type)

                if len0 == 1:
                    digest_side = 0
                else:
                    digest_side = 1

                # from result dict
                digest_pos = int(caps_result[digest_side][enz_res_type][0])

                # {'EcoRI': {
                #   'pattern': 'G^AATT_C', 'diggno': 1, 'digpos': 17}
                # }
                self.caps_dict[enzyme_name] = {
                    'pattern': enz_res_type.elucidate(),
                    'diggno': digest_side,
                    'digpos': digest_pos}

        #log.debug("{}".format(self.caps_dict))

        # 見つかった数
        self.caps_found = len(self.caps_dict)
        if self.caps_found > 0:
            self.marker_available = True
        else:
            # 初期化かな
            pass


    @classmethod
    def get_caps_result(cls, seq_target, enzyme_list):

        # http://biopython.org/DIST/docs/cookbook/Restriction.html
        # 2.6 Analysing sequences with a RestrictionBatch
        ar_seq = Seq(seq_target, IUPACAmbiguousDNA())
        rb = Restriction.RestrictionBatch(enzyme_list)
        # If linear is False, the restriction sites that span over
        # the boundaries will be included.
        caps_result_dict = rb.search(ar_seq, linear=True)

        caps_result_dict_str = dict()

#    log.debug("caps_result_dict {}".format(caps_result_dict))

        # convert enzyme class from RestrictionType to string
        for enzyme_RestrictionType in caps_result_dict.keys():
            enzyme_string = str(enzyme_RestrictionType)
            caps_result_dict_str[enzyme_string] = \
                caps_result_dict[enzyme_RestrictionType]

#        log.debug("{}".format(type(enzyme_RestrictionType)))
#        log.debug("{}".format(str(enzyme_RestrictionType)))
#        log.debug("{}".format(type(str(enzyme_RestrictionType))))

#    sys.exit(1)

#    log.debug("{} {}".format(caps_result_dict, caps_result_dict_str))

        return caps_result_dict, caps_result_dict_str


    def _pick_ref_around_seq(self):

        self.vseq_ref = self.vseq_ano[0]
        self.vseq_ref_len = len(self.vseq_ref)

        # abs pos決め
        # pos=60, 60-1=59
        self.abs_around_seq_pre_end = self.pos - 1

        # 59-10+1 = 50
        self.abs_around_seq_pre_stt = \
            self.abs_around_seq_pre_end - glv.AROUND_SEQ_LEN + 1

        # 60+16 = 76
        # 60+1 = 61
        self.abs_around_seq_aft_stt = self.pos + self.vseq_ref_len

        # 76+10-1=85
        self.abs_around_seq_aft_end = \
            self.abs_around_seq_aft_stt + glv.AROUND_SEQ_LEN - 1

        # preの切り出し
        self.around_seq_pre = glv.ref.pick_refseq(
            self.chrom,
            self.abs_around_seq_pre_stt,
            self.abs_around_seq_pre_end).upper()

        # aftの切り出し
        self.around_seq_aft = glv.ref.pick_refseq(
            self.chrom,
            self.abs_around_seq_aft_stt,
            self.abs_around_seq_aft_end).upper()

        self.around_seq_pre_len = len(self.around_seq_pre)
        self.around_seq_aft_len = len(self.around_seq_aft)

        self.seq_target_ref = "{}{}{}".format(
            self.around_seq_pre,
            self.vseq_ref,
            self.around_seq_aft)

        self.seq_target_ref_len = len(self.seq_target_ref)


    def _prepare_from_variant(self, variant_df_row, distin_dict):

        hdr_dict = distin_dict['variant']['hdr_dict']

        self.chrom = str(variant_df_row[hdr_dict['chrom']])
        self.pos = int(variant_df_row[hdr_dict['pos']])

        # cl1/cl2
        self.targ_grp = str(variant_df_row[hdr_dict['targ_grp']])
        self.g0_name, self.g1_name = self.targ_grp.split(',')
        self.gname_gno = [self.g0_name, self.g1_name]

        # 00/01,hohe_s1,1.1/1.1
        self.gts_segr_lens = str(variant_df_row[hdr_dict['gts_segr_lens']])

        # 1/1
        self.set_n_str = str(variant_df_row[hdr_dict['set_n']])
        self.set_my_no, self.set_n = map(int, self.set_n_str.split('/'))

        # targ_ano
        self.targ_ano = str(variant_df_row[hdr_dict['targ_ano']])
        self.g0_ano, self.g1_ano = map(int, self.targ_ano.split(','))
        # convert table
        self.ano_gno = [self.g0_ano, self.g1_ano]

        # 1,1,0,-1
        self.len_g0g1_dif_long = str(
            variant_df_row[hdr_dict['len_g0g1_dif_long']])
        self.g0_len, self.g1_len, self.diff_len, self.long_gno = \
            map(int, self.len_g0g1_dif_long.split(','))

        if self.long_gno == glv.SAME_LENGTH:
            self.long_gno = 0
            self.short_gno = 1
        else:
            self.short_gno = 0 if self.long_gno == 1 else 1

        self.var_type = str(variant_df_row[hdr_dict['var_type']])

        # T,G access to all vseq
        self.vseq_ano_str = str(variant_df_row[hdr_dict['vseq_ano_str']])
        self.vseq_ano = self.vseq_ano_str.split(',')
        self.vseq_lens_ano = [len(vseq_ano) for vseq_ano in self.vseq_ano]


    def make_seq_template_ref(self):

        #self.seq_template_ref = ''
        #self.seq_template_ref_len = 0

        # frag_padを切り出す
        self.fragment_pad_len = glv.conf.fragment_pad_len

        # abs_posを決める
        self.abs_frag_pad_pre_end = self.abs_around_seq_pre_stt - 1

        self.abs_frag_pad_pre_stt = \
            self.abs_frag_pad_pre_end - self.fragment_pad_len + 1

        self.abs_frag_pad_aft_stt = self.abs_around_seq_aft_end + 1
        self.abs_frag_pad_aft_end = \
            self.abs_frag_pad_aft_stt + self.fragment_pad_len - 1

        # templateの絶対pos stringを最初に作る。
        # 最初はposはすべて完成しているが、今後端を切るために。
        self.seq_template_ref_abs_pos = \
            self._make_abs_pos(
                self.abs_frag_pad_pre_stt,
                self.abs_frag_pad_aft_end)

        # これは最初の
        # templateの相対pos string
        self.seq_template_ref_rel_pos, \
        self.SEQUENCE_TARGET = \
            self._convert_to_rel_pos(self.seq_template_ref_abs_pos)

        # refのfragpadを取り出す。
        self.frag_pad_pre, self.frag_pad_aft, self.seq_template_ref = \
            self._get_seq_template_ref(
                self.chrom, self.seq_template_ref_abs_pos)

        self.frag_pad_pre_len = len(self.frag_pad_pre)
        self.frag_pad_aft_len = len(self.frag_pad_aft)
        self.seq_template_ref_len = len(self.seq_template_ref)

        # groupのproductにもセットする
        for gno in range(2):
            self.gr_prod[gno].set_frag_pad(
                self.frag_pad_pre, self.frag_pad_aft)


    def _make_abs_pos(
        self,
        abs_frag_pad_pre_stt,
        abs_frag_pad_aft_end):

        return "{}/{}/{}/{}/{}/{}/{}/{}/{}".format(
            abs_frag_pad_pre_stt,
            self.abs_frag_pad_pre_end,
            self.abs_around_seq_pre_stt,
            self.abs_around_seq_pre_end,
            self.pos,
            self.abs_around_seq_aft_stt,
            self.abs_around_seq_aft_end,
            self.abs_frag_pad_aft_stt,
            abs_frag_pad_aft_end)


    def _convert_to_rel_pos(self, seq_template_ref_abs_pos):

        abs_frag_pad_pre_stt, abs_frag_pad_pre_end, \
        abs_around_seq_pre_stt, abs_around_seq_pre_end, \
        abs_pos, \
        abs_around_seq_aft_stt, abs_around_seq_aft_end, \
        abs_frag_pad_aft_stt, abs_frag_pad_aft_end = \
            self._separate_pos_str(seq_template_ref_abs_pos)

        # 11......21
        # 1.......11
        dif = abs_frag_pad_pre_stt

        rel_frag_pad_pre_stt = 1
        rel_frag_pad_pre_end = abs_frag_pad_pre_end - dif + 1
        rel_around_seq_pre_stt = abs_around_seq_pre_stt - dif + 1
        rel_around_seq_pre_end = abs_around_seq_pre_end - dif + 1
        rel_pos = abs_pos - dif + 1
        rel_around_seq_aft_stt = abs_around_seq_aft_stt - dif + 1
        rel_around_seq_aft_end = abs_around_seq_aft_end - dif + 1
        rel_frag_pad_aft_stt = abs_frag_pad_aft_stt - dif + 1
        rel_frag_pad_aft_end = abs_frag_pad_aft_end - dif + 1

        #log.debug("{} {} {} {}".format(
        #    self.SEQUENCE_TARGET,
        #    rel_around_seq_pre_stt,
        #    rel_around_seq_aft_end,
        #    rel_around_seq_pre_stt))

        # 100-1=99+1=100
        sequence_target = "{},{}".format(
            rel_around_seq_pre_stt,
            rel_around_seq_aft_end - rel_around_seq_pre_stt + 1)

        return "{}/{}/{}/{}/{}/{}/{}/{}/{}".format(
            rel_frag_pad_pre_stt,
            rel_frag_pad_pre_end,
            rel_around_seq_pre_stt,
            rel_around_seq_pre_end,
            rel_pos,
            rel_around_seq_aft_stt,
            rel_around_seq_aft_end,
            rel_frag_pad_aft_stt,
            rel_frag_pad_aft_end), \
            sequence_target


    def _separate_pos_str(self, pos_str):

        return map(int, pos_str.split('/'))


    def _get_seq_template_ref(self, chrom, seq_template_ref_abs_pos):

        abs_frag_pad_pre_stt, abs_frag_pad_pre_end, \
        abs_around_seq_pre_stt, abs_around_seq_pre_end, \
        abs_pos, \
        abs_around_seq_aft_stt, abs_around_seq_aft_end, \
        abs_frag_pad_aft_stt, abs_frag_pad_aft_end = \
            self._separate_pos_str(seq_template_ref_abs_pos)

        # update
        # pick frag_pad_pre
        frag_pad_pre = glv.ref.pick_refseq(
            chrom,
            abs_frag_pad_pre_stt,
            abs_frag_pad_pre_end).upper()

        # これは変わらない
        seq_target_ref = self.seq_target_ref

        # pick frag_pad_aft
        frag_pad_aft = glv.ref.pick_refseq(
            chrom,
            abs_frag_pad_aft_stt,
            abs_frag_pad_aft_end).upper()

        seq_template_ref = "{}{}{}".format(
            frag_pad_pre,
            seq_target_ref,
            frag_pad_aft)

        return frag_pad_pre, frag_pad_aft, seq_template_ref


    def adjust_seq_temlate_ref_by_enzyme(self):
        '''
        '''

        enzyme_cnt_per_variant = 0

        l_marker_info = list()
        l_seq_t_abs_pos = list()
        l_seq_t_rel_pos = list()
        l_seq_t_ref = list()
        l_seq_target = list()

        if self.var_type == glv.INDEL:

            enzyme_cnt_per_variant = 1
            marker_info = self._make_marker_info()

            l_marker_info = [marker_info]
            l_seq_t_abs_pos = [self.seq_template_ref_abs_pos]
            l_seq_t_rel_pos = [self.seq_template_ref_rel_pos]
            l_seq_t_ref = [self.seq_template_ref]
            l_seq_target = [self.SEQUENCE_TARGET]

        else:
            # duplicate line information by enzyme
            enzyme_cnt_per_variant = len(self.caps_dict)
            self._divide_information_by_enzyme(
                enzyme_cnt_per_variant,
                l_marker_info,
                l_seq_t_abs_pos,
                l_seq_t_rel_pos,
                l_seq_t_ref,
                l_seq_target)

        line_for_each_enzyme = list()        

        for num in range(enzyme_cnt_per_variant):
            enzyme_cnt = "{}/{}".format(num+1, enzyme_cnt_per_variant)
            set_enz_cnt = "{}-{}".format(self.set_n_str, enzyme_cnt)
            vseq_lens_ano_str = \
                "{}".format(','.join(map(str, self.vseq_lens_ano)))

            l_list = list()

            # out to marker out file
            l_list += [self.chrom]
            l_list += [self.pos]
            l_list += [self.targ_grp]
            l_list += [self.gts_segr_lens]

            l_list += [set_enz_cnt]
            l_list += [self.targ_ano]
            l_list += [self.var_type]
            l_list += [l_marker_info[num]]

            # vseq_lens_ano
            l_list += [vseq_lens_ano_str]

            # 1) g0_seq_target_len
            l_list += [self.gr_prod[0].seq_target_len]
            # 2) g0_seq_target
            l_list += [self.gr_prod[0].seq_target]

            # 1) g1_seq_target_len
            l_list += [self.gr_prod[1].seq_target_len]
            # 2) g1_seq_target
            l_list += [self.gr_prod[1].seq_target]

            # 1) seq_template_ref_len
            l_list += [len(l_seq_t_ref[num])]
            # 2) seq_template_ref_abs_pos
            l_list += [l_seq_t_abs_pos[num]]
            # 3) seq_template_ref_rel_pos
            l_list += [l_seq_t_rel_pos[num]]
            # 4) SEQUENCE_TARGET
            l_list += [l_seq_target[num]]
            # 5) seq_template_ref
            l_list += [l_seq_t_ref[num]]

            line_for_each_enzyme.append('\t'.join(map(str, l_list)))

        self.line = '\n'.join(map(str, line_for_each_enzyme))


    def _divide_information_by_enzyme(
        self,
        enzyme_cnt_per_variant,
        l_marker_info,
        l_seq_t_abs_pos,
        l_seq_t_rel_pos,
        l_seq_t_ref,
        l_seq_target):

        #log.debug("start")
        #log.debug("{}".format(self.caps_dict))

        for enzyme in self.caps_dict:

            # enzymeごとに、digest_posは、grに組み込む
            marker_info = self._make_marker_info(enzyme)

            l_marker_info.append(marker_info)

            # 何もなければすでにrefとしてセットされた値を
            seq_t_abs_pos = self.seq_template_ref_abs_pos
            seq_t_rel_pos = self.seq_template_ref_rel_pos
            seq_t_ref = self.seq_template_ref
            SEQUENCE_TARGET = self.SEQUENCE_TARGET

            # capsを検索する
            caps_result_dict, caps_result_dict_str = \
                EvalVariant.get_caps_result(
                    seq_t_ref, [enzyme])

            caps_result_cnt = len(caps_result_dict_str[enzyme])

            if caps_result_cnt != 0:
                # 切断されているなら、seq_templateを更新し、
                # 現在、
                # ref.
                # self.gr_prod[num]


                # 20200715
                # それぞれ、メソッドで使っているだけ
                seq_t_abs_pos = \
                    self._change_abs_pos_by_digest(
                        caps_result_dict_str[enzyme])

                # templateの相対pos string
                seq_t_rel_pos, SEQUENCE_TARGET = \
                    self._convert_to_rel_pos(seq_t_abs_pos)

                # まとめてpick
                frag_pad_pre, frag_pad_aft, seq_t_ref = \
                    self._get_seq_template_ref(
                        self.chrom, seq_t_abs_pos)

                # 確認用 capsを検索する
                #caps_result_dict, caps_result_dict_str = \
                #    Products.get_caps_result(
                #        seq_t_ref, [enzyme])
                #log.debug("rechecked {}".format(caps_result_dict_str))
                #log.debug("{}\n".format(seq_t_rel_pos))

            else:
                pass

            l_seq_t_abs_pos.append(seq_t_abs_pos)
            l_seq_t_rel_pos.append(seq_t_rel_pos)
            l_seq_t_ref.append(seq_t_ref)
            l_seq_target.append(SEQUENCE_TARGET)


    def serial_enz(self, enzyme):

        return  "{},{},{},{}".format(
            enzyme,
            self.caps_dict[enzyme]['diggno'],
            self.caps_dict[enzyme]['digpos'],
            self.caps_dict[enzyme]['pattern'])


    def _make_marker_info(self, enzyme=''):

        marker_info = ''

        if self.var_type == glv.INDEL:
            # comma, and only 2 info
            marker_info = "{},{}".format(self.long_gno, self.diff_len)

        else:
            # comma, and 3 info
            marker_info = "{},{},{},{}".format(
                enzyme,
                self.caps_dict[enzyme]['diggno'],
                self.caps_dict[enzyme]['digpos'],
                self.caps_dict[enzyme]['pattern'])

        return marker_info


    def _change_abs_pos_by_digest(self, caps_dig_pos):

        #log.debug("{}".format(caps_dig_pos))

        fixed_pre_stt = self.abs_frag_pad_pre_stt
        fixed_pre_end = self.abs_frag_pad_pre_end
        fixed_aft_stt = self.abs_frag_pad_aft_stt
        fixed_aft_end = self.abs_frag_pad_aft_end

        five_prime_biggest_pos = fixed_pre_stt
        three_prime_smallest_pos = fixed_aft_end

        # EcoRI': {'pattern': 'G^AATT_C'}
        # [17]
        # 123456789012345678901234567890
        # CTCTGTTCGGTGGAAGAATTCAGATTTCAGAGTCA
        #               G^AATT_C
        #                 /-> 17
        #                 AATTCAGATTTCAGAGTCA
        # 切断されるポイントは17。これは残る側。
        # 残る側に、切断ポイントを残していいのか。
        # 今は残している。

        # digest_positionごとに調査
        for rel_digest_pos in caps_dig_pos:
            # 絶対posに変換
            # 10001    100 -> 10001+100 - 10100
            abs_digest_pos = fixed_pre_stt + rel_digest_pos - 1

            #log.debug("rel={} abs={} pstt<{} pend>{} astt<{}".format(
            #    rel_digest_pos,
            #    abs_digest_pos,
            #    fixed_pre_stt,
            #    fixed_pre_end,
            #    fixed_aft_stt))

            # |fixed_pre_stt
            #                |fixed_pre_end
            #                               |fixed_aft_stt
            # <--------------><=============<-------------->

            # １つポジションをずらす。それにより認識サイトが
            # 壊れる
            if abs_digest_pos < fixed_pre_end:
                # 5'側では、一つ先で切る
                # always
                five_prime_biggest_pos = abs_digest_pos + 1

            elif fixed_aft_stt < abs_digest_pos:
                # only once
                # 3'側では、一つ手前で切る
                three_prime_smallest_pos = abs_digest_pos - 1
                break

        seq_template_ref_abs_pos = \
            self._set_seq_template_ref_abs_pos(
                five_prime_biggest_pos,
                three_prime_smallest_pos)

        #log.debug("{}".format(self.seq_template_ref_abs_pos))
        #log.debug("{}".format(seq_template_ref_abs_pos))

        return seq_template_ref_abs_pos


    def _get_seq_template_ref(self, chrom, seq_template_ref_abs_pos):

        abs_frag_pad_pre_stt, abs_frag_pad_pre_end, \
        abs_around_seq_pre_stt, abs_around_seq_pre_end, \
        abs_pos, \
        abs_around_seq_aft_stt, abs_around_seq_aft_end, \
        abs_frag_pad_aft_stt, abs_frag_pad_aft_end = \
            self._separate_pos_str(seq_template_ref_abs_pos)

        # update
        # pick frag_pad_pre
        frag_pad_pre = glv.ref.pick_refseq(
            chrom,
            abs_frag_pad_pre_stt,
            abs_frag_pad_pre_end).upper()

        # これは変わらない
        seq_target_ref = self.seq_target_ref

        # pick frag_pad_aft
        frag_pad_aft = glv.ref.pick_refseq(
            chrom,
            abs_frag_pad_aft_stt,
            abs_frag_pad_aft_end).upper()

        seq_template_ref = "{}{}{}".format(
            frag_pad_pre,
            seq_target_ref,
            frag_pad_aft)

        return frag_pad_pre, frag_pad_aft, seq_template_ref


    def _set_seq_template_ref_abs_pos(
        self,
        abs_frag_pad_pre_stt,
        abs_frag_pad_aft_end):

        return "{}/{}/{}/{}/{}/{}/{}/{}/{}".format(
            abs_frag_pad_pre_stt,
            self.abs_frag_pad_pre_end,
            self.abs_around_seq_pre_stt,
            self.abs_around_seq_pre_end,
            self.pos,
            self.abs_around_seq_aft_stt,
            self.abs_around_seq_aft_end,
            self.abs_frag_pad_aft_stt,
            abs_frag_pad_aft_end)


