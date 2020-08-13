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
from joblib import Parallel, delayed

from vprimer.enzyme import Enzyme
from vprimer.eval_variant import EvalVariant


class Marker(object):

    def __init__(self):

        self.enzyme = Enzyme()

    def design_marker(self):

        # progress check
        if utl.progress_check('marker') == False:
            log.info("progress={} so skip variant.".format(
                glv.conf.progress))
            return
        log.info("Start processing {}".format('marker'))

        # primer3用フラグメントを作成する
        # for each distinguish_groups
        for distin_dict in glv.outlist.distin_files:

            # read variant file 
            variant_file = distin_dict['variant']['out_path']
            log.info("variant_file {}".format(variant_file))

            df_distin = pd.read_csv(
                variant_file, sep='\t', header=0, index_col=None)

            # Bio.Restriction.Restriction_Dictionary
            self.enzyme.read_enzyme_file()

            # file name to write out result to text
            out_txt_file = distin_dict['marker']['out_path']
            utl.save_to_tmpfile(out_txt_file)

            start = time.time()
            with open(out_txt_file, mode='a') as f:

                # write header
                #f.write("{}\n".format(distin_dict['marker']['hdr_text']))

                if glv.conf.parallel == True:
                    log.info("do Parallel cpu {} parallel {}".format(
                        glv.conf.thread,
                        glv.conf.parallele_full_thread))

                    Parallel(
                        n_jobs=glv.conf.parallele_full_thread,
                        backend="threading")(
                        [
                            delayed(self._loop_evaluate_for_marker)
                                (distin_dict, variant_df_row, f) \
                                for variant_df_row in df_distin.itertuples()
                        ]
                    )

                else:
                    log.info("do Serial cpu 1")

                    # each variant
                    for variant_df_row in df_distin.itertuples():
                        # バリアントがマーカーとして使えるかどうか、判断する。
                        # マーカー化可能なものはprimer3用の情報を準備する。
                        self._loop_evaluate_for_marker(
                            distin_dict, variant_df_row, f)

            utl.sort_file(
                'marker', distin_dict, out_txt_file,
                'chrom', 'pos', 'marker_info', 'string')

            log.info("marker {} {}".format(
                utl.elapsed_time(time.time(), start),
                distin_dict['marker']['base_nam']))


    def _loop_evaluate_for_marker(self, distin_dict, variant_df_row, f):

        evalv = EvalVariant()
        evalv.evaluate_for_marker(
            variant_df_row, distin_dict, self.enzyme.enzyme_list)
        
        #sys.exit(1)

        if evalv.marker_available == True:

            # マーカー化ならseq_template_ref を作成
            evalv.make_seq_template_ref()

            # enzymeごとにdigest posを調整したseq_templateを作る
            evalv.adjust_seq_temlate_ref_by_enzyme()

            # 書き出す
            f.write("{}\n".format(evalv.line))
            evalv.line = ''

