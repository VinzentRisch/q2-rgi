import os
import shutil
from unittest.mock import MagicMock, patch

from qiime2.plugin.testing import TestPluginBase

from q2_amr.card.kmer import kmer_query_mags_card, run_rgi_kmer_query
from q2_amr.types import (
    CARDAnnotationDirectoryFormat,
    CARDDatabaseDirectoryFormat,
    CARDKmerDatabaseDirectoryFormat,
)
from q2_amr.types._format import CARDMAGsKmerAnalysisDirectoryFormat


class CARDDKmerDatabaseFormat:
    pass


class TestKmer(TestPluginBase):
    package = "q2_amr.tests"

    def copy_analysis_file(
        self, tmp, input_file, input_type, kmer_size, minimum, threads
    ):
        des_path = os.path.join(tmp, "output_61mer_analysis_rgi_summary.txt")
        src_path = self.get_data_path("kmer_analysis_rgi_summary.txt")
        shutil.copy(src_path, des_path)

    def test_kmer_query_mags_card(self):
        mock_run_rgi_kmer_query = MagicMock(side_effect=self.copy_analysis_file)
        with patch(
            "q2_amr.card.kmer.run_rgi_kmer_query", side_effect=mock_run_rgi_kmer_query
        ), patch("q2_amr.card.kmer.load_card_db", return_value="61"):
            amr_annotations = CARDAnnotationDirectoryFormat()
            card_db = CARDDatabaseDirectoryFormat()
            kmer_db = CARDKmerDatabaseDirectoryFormat()
            annotation_dir = os.path.join(str(amr_annotations), "sample1", "bin1")
            os.makedirs(annotation_dir)
            des_path = os.path.join(annotation_dir, "amr_annotation.json")
            src_path = self.get_data_path("rgi_output.json")
            shutil.copy(src_path, des_path)
            result = kmer_query_mags_card(
                amr_annotations=amr_annotations, card_db=card_db, kmer_db=kmer_db
            )
            self.assertIsInstance(result, CARDMAGsKmerAnalysisDirectoryFormat)
            self.assertTrue(
                os.path.exists(
                    os.path.join(
                        str(result), "sample1", "bin1", "61mer_analysis_mags.txt"
                    )
                )
            )

    def test_run_rgi_kmer_query(self):
        with patch("q2_amr.card.kmer.run_command") as mock_run_command:
            run_rgi_kmer_query(
                tmp="path_tmp",
                input_file="path_inout_file",
                input_type="rgi",
                kmer_size="61",
                minimum="10",
                threads="4",
            )
            mock_run_command.assert_called_once_with(
                [
                    "rgi",
                    "kmer_query",
                    "--input",
                    "path_inout_file",
                    "--rgi",
                    "--kmer_size",
                    "61",
                    "--minimum",
                    "10",
                    "--threads",
                    "4",
                    "--output",
                    "output",
                    "--local",
                ],
                "path_tmp",
                verbose=True,
            )
