# -*- coding: utf-8 -*-

import sys
import os
import errno
import time

# global configuration
import vprimer.glv as glv
import vprimer.utils as utl

from vprimer.logging_config import LogConf


class OutList(object):

    def __init__(self):

        self.outf_prefix = {
            'prepare'       : {'no':  0,  'fn': ''},
            'variant'       : {'no': 10,  'fn': '010_variant'},
            'marker'        : {'no': 20,  'fn': '020_marker'},
            'primer'        : {'no': 30,  'fn': '030_primer'},
            'formfail'      : {'no': 40,  'fn': '040_formatF'},
            'formsafe'      : {'no': 50,  'fn': '050_formatS'},
        }

        # list of full path of out text
        self.distin_files = list()

    def open_log(self):

        global log
        log = LogConf.open_log(__name__)


    def prepare_distin_files(self):
        ''' access distinct_group information with outfile info
        '''

        for distin_dict in glv.conf.distin_g_list:
        # distin_dict={
        #    0: 'gAkenohoshi',
        #    1: 'gAkitakomachi',
        #    'region': ['rg0'],
        #    'pick_mode': 'all'}.

            for reg in distin_dict['region']:

                # reconstruct cause region is list
                distin_out_dict = {
                    0: distin_dict[0],
                    1: distin_dict[1],
                    'region': reg,
                    'pick_mode': distin_dict['pick_mode'],
                }

                out_file_dict = dict()

                for key, no_fn in self.outf_prefix.items():

                    # file name to write out result to text
                    out_file_path, base_name = self._make_distin_fname(
                        no_fn['fn'],
                        distin_dict[0],             # g0
                        distin_dict[1],             # g1
                        reg,                        # region
                        distin_dict['pick_mode'])   # pick_mode

                    hdr_text, hdr_list, hdr_dict = self.make_header(key)

                    out_file_dict[key] = {
                        'out_path': out_file_path,
                        'base_nam': base_name,
                        'hdr_text': hdr_text,
                        'hdr_list': hdr_list,
                        'hdr_dict': hdr_dict,
                    }

                # add dictionary item
                distin_out_dict.update(out_file_dict)
            # make list of dictionary
            self.distin_files.append(distin_out_dict)

        #log.debug("{}".format(self.distin_files))

# [
# 	{
# 		0: 'gAkenohoshi', 
# 		1: 'gAkitakomachi', 
# 		'region': 'rg0', 
# 		'pick_mode': 'all', 
# 		'variant': 
# 			{
# 
# 				'path': '/lustre7/home/lustre4/ibrcuser/software/vprimer_v1/vprimer/out_vprimer1/01_variant~distin~gAkenohoshi~gAkitakomachi~rg0~all~20-200.txt', 
# 
# 				'hdr_txt': 'CHROM\tPOS\tvar_type\tg0_ano\tg1_ano\tset_n\tALT_vc\tlens\tg0_name\tg0_altyp\tg0_alint\tg0_gtlst\tg0_lens\tg0_longest\tg1_name\tg1_altyp\tg1_alint\tg1_gtlst\tg1_lens\tg1_longest\tsegr_ptn\tlong_side\tdiff_len\tano_vseq_str', 
# 
# 				'hdr_lst': ['CHROM', 'POS', 'var_type', 'g0_ano', 'g1_ano', 'set_n', 'ALT_vc', 'lens', 'g0_name', 'g0_altyp', 'g0_alint', 'g0_gtlst', 'g0_lens', 'g0_longest', 'g1_name', 'g1_altyp', 'g1_alint', 'g1_gtlst', 'g1_lens', 'g1_longest', 'segr_ptn', 'long_side', 'diff_len', 'ano_vseq_str'], 
# 
# 				'hdr_dct': {'CHROM': 1, 'POS': 2, 'var_type': 3, 'g0_ano': 4, 'g1_ano': 5, 'set_n': 6, 'ALT_vc': 7, 'lens': 8, 'g0_name': 9, 'g0_altyp': 10, 'g0_alint': 11, 'g0_gtlst': 12, 'g0_lens': 13, 'g0_longest': 14, 'g1_name': 15, 'g1_altyp': 16, 'g1_alint': 17, 'g1_gtlst': 18, 'g1_lens': 19, 'g1_longest': 20, 'segr_ptn': 21, 'long_side': 22, 'diff_len': 23, 'ano_vseq_str': 24}
# 
# 			},
# 
# 		'fragment': 
# 			{
# 
# 				'path': '/lustre7/home/lustre4/ibrcuser/software/vprimer_v1/vprimer/out_vprimer1/02_fragment~distin~gAkenohoshi~gAkitakomachi~rg0~all~20-200.txt', 
# 
# 				'hdr_txt': 'CHROM\tPOS\tvar_type\tg0_ano\tg1_ano\tset_n\tALT_vc\tlens\tg0_name\tg0_altyp\tg0_alint\tg0_gtlst\tg0_lens\tg0_longest\tg1_name\tg1_altyp\tg1_alint\tg1_gtlst\tg1_lens\tg1_longest\tsegr_ptn\tlong_side\tdiff_len\tano_vseq_str', 
# 
# 				'hdr_lst': ['CHROM', 'POS', 'var_type', 'g0_ano', 'g1_ano', 'set_n', 'ALT_vc', 'lens', 'g0_name', 'g0_altyp', 'g0_alint', 'g0_gtlst', 'g0_lens', 'g0_longest', 'g1_name', 'g1_altyp', 'g1_alint', 'g1_gtlst', 'g1_lens', 'g1_longest', 'segr_ptn', 'long_side', 'diff_len', 'ano_vseq_str'], 
# 
# 				'hdr_dct': {'CHROM': 1, 'POS': 2, 'var_type': 3, 'g0_ano': 4, 'g1_ano': 5, 'set_n': 6, 'ALT_vc': 7, 'lens': 8, 'g0_name': 9, 'g0_altyp': 10, 'g0_alint': 11, 'g0_gtlst': 12, 'g0_lens': 13, 'g0_longest': 14, 'g1_name': 15, 'g1_altyp': 16, 'g1_alint': 17, 'g1_gtlst': 18, 'g1_lens': 19, 'g1_longest': 20, 'segr_ptn': 21, 'long_side': 22, 'diff_len': 23, 'ano_vseq_str': 24
# 
# 			}
# 		}
# 	}
# ]


    def _make_distin_fname(
        self, outf_pref, distin_0, distin_1, rg, pick_mode):
        """
        """

        #distin~gHitomebore~gKaluheenati~rg0~all~50-200.txt
        base_name = "{}~{}~{}~{}~{}~i{}-{}~p{}-{}".format(
                outf_pref,
                distin_0,
                distin_1,
                rg,
                pick_mode,
                glv.conf.min_indel_len,
                glv.conf.max_indel_len,
                glv.conf.min_product_size,
                glv.conf.max_product_size)

        out_file_path = "{}/{}.txt".format(glv.conf.out_dir, base_name)

        return out_file_path, base_name


    def make_header(self, type):

        hd_l = list()
        hd_d = dict()
        i = 1   # 0 is index

        if type == 'variant':
            # variant_id(6): chrom, pos, var_type, g0_ano, g1_ano, set_n
            hd_l, hd_d, i = self._mkmap('chrom',         hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('pos',           hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('targ_grp',      hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('gts_segr_lens', hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('set_n',         hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('targ_ano',      hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('len_g0g1_dif_long', hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('var_type',      hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('vseq_ano_str',  hd_l, hd_d, i)

        elif type == 'marker':
            hd_l, hd_d, i = self._mkmap('chrom',         hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('pos',           hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('targ_grp',      hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('gts_segr_lens', hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('set_enz_cnt',   hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('targ_ano',      hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('var_type',      hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('marker_info',   hd_l, hd_d, i)

            hd_l, hd_d, i = self._mkmap('vseq_lens_ano_str', hd_l, hd_d, i)

            # g0
            hd_l, hd_d, i = self._mkmap('g0_seq_target_len',
                 hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('g0_seq_target',
                 hd_l, hd_d, i)
            # g1
            hd_l, hd_d, i = self._mkmap('g1_seq_target_len',
                 hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('g1_seq_target',
                 hd_l, hd_d, i)

            # seq_template_ref
            hd_l, hd_d, i = self._mkmap('seq_template_ref_len',
                 hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('seq_template_ref_abs_pos',
                 hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('seq_template_ref_rel_pos',
                 hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('SEQUENCE_TARGET',
                 hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('seq_template_ref',
                 hd_l, hd_d, i)

        elif type == 'primer':  # there is no header

            hd_l, hd_d, i = self._mkmap('marker_id',   hd_l, hd_d, i)

            hd_l, hd_d, i = self._mkmap('chrom',       hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('pos',         hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('targ_grp',    hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('gts_segr_lens', hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('targ_ano',  hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('set_enz_cnt', hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('var_type',    hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('marker_info', hd_l, hd_d, i)

            hd_l, hd_d, i = self._mkmap('vseq_lens_ano_str', hd_l, hd_d, i)

            hd_l, hd_d, i = self._mkmap('try_cnt',     hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('complete',    hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('blast_check', hd_l, hd_d, i)

            hd_l, hd_d, i = self._mkmap('target_gno',  hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('target_len',  hd_l, hd_d, i)

            hd_l, hd_d, i = self._mkmap('g0_seq_target_len',
                hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('g0_seq_target',
                hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('g1_seq_target_len',
                hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('g1_seq_target',
                hd_l, hd_d, i)

            hd_l, hd_d, i = self._mkmap('seq_template_ref_len',
                hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('seq_template_ref_abs_pos',
                hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('seq_template_ref_rel_pos',
                hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('PRIMER_PAIR_0_PRODUCT_SIZE',
                hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('PRIMER_LEFT_0',
                hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('left_primer_id',
                hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('PRIMER_LEFT_0_SEQUENCE',
                hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('PRIMER_RIGHT_0',
                hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('right_primer_id',
                hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('PRIMER_RIGHT_0_SEQUENCE',
                hd_l, hd_d, i)

            hd_l, hd_d, i = self._mkmap('SEQUENCE_TARGET',
                hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('SEQUENCE_EXCLUDED_REGION',
                hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('SEQUENCE_TEMPLATE',
                hd_l, hd_d, i)

        elif type == 'formsafe':

            hd_l, hd_d, i = self._mkmap('chrom',            hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('pos',              hd_l, hd_d, i)

            hd_l, hd_d, i = self._mkmap('try_cnt',          hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('complete',         hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('comment',          hd_l, hd_d, i)

            hd_l, hd_d, i = self._mkmap('var_type',         hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('gts_segr_lens',    hd_l, hd_d, i)

            hd_l, hd_d, i = self._mkmap('enzyme',           hd_l, hd_d, i)

            hd_l, hd_d, i = self._mkmap('g0_name',          hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('g0_product_size',  hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('g0_digested_size', hd_l, hd_d, i)

            hd_l, hd_d, i = self._mkmap('diff_length',      hd_l, hd_d, i)

            hd_l, hd_d, i = self._mkmap('g1_digested_size', hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('g1_product_size',  hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('g1_name',          hd_l, hd_d, i)

            hd_l, hd_d, i = self._mkmap('left_primer_id',   hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('PRIMER_LEFT_0_SEQUENCE',
                hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('right_primer_id',  hd_l, hd_d, i)
            hd_l, hd_d, i = self._mkmap('PRIMER_RIGHT_0_SEQUENCE',
                hd_l, hd_d, i)


        hd_str = '\t'.join(map(str, hd_l))
        return hd_str, hd_l, hd_d


    def _mkmap(self, key, header_list, header_dict, idx):

        #log.debug("{} {} {} {}".format(
        #    key, header, header_dict, idx,
        #    ))

        header_list += [key]
        add_dict = {key: idx}
        header_dict.update(add_dict)
        idx += 1

        return header_list, header_dict, idx


    def complete_list(self):
        pass

