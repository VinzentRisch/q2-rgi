import os
import shutil
from unittest.mock import MagicMock, patch

from qiime2.plugin.testing import TestPluginBase

from q2_amr.card.kmer import (
    kmer_build_card,
    kmer_query_mags_card,
    kmer_query_reads_card,
    run_rgi_kmer_build,
    run_rgi_kmer_query,
)
from q2_amr.types import (
    CARDAnnotationDirectoryFormat,
    CARDDatabaseDirectoryFormat,
    CARDKmerDatabaseDirectoryFormat,
)
from q2_amr.types._format import (
    CARDGeneAnnotationDirectoryFormat,
    CARDMAGsKmerAnalysisDirectoryFormat,
    CARDReadsKmerAnalysisDirectoryFormat,
)


class TestKmer(TestPluginBase):
    package = "q2_amr.tests"

    def copy_analysis_file(
        self, tmp, input_file, input_type, kmer_size, minimum, threads
    ):
        des_path = os.path.join(tmp, f"output_61mer_analysis_{input_type}_summary.txt")
        src_path = self.get_data_path(f"kmer_analysis_{input_type}_summary.txt")
        shutil.copy(src_path, des_path)

    def _run_kmer_query_test(self, annotation_format, output_format, query_function):
        mock_run_rgi_kmer_query = MagicMock(side_effect=self.copy_analysis_file)
        with patch(
            "q2_amr.card.kmer.run_rgi_kmer_query", side_effect=mock_run_rgi_kmer_query
        ), patch("q2_amr.card.kmer.load_card_db", return_value="61"):
            amr_annotations = annotation_format()
            card_db = CARDDatabaseDirectoryFormat()
            kmer_db = CARDKmerDatabaseDirectoryFormat()
            if query_function == kmer_query_reads_card:
                annotation_dir = os.path.join(str(amr_annotations), "sample1")
                des_path = os.path.join(annotation_dir, "sorted.length_100.bam")
                src_path = self.get_data_path("output.sorted.length_100.bam")
                analysis_path = os.path.join("sample1", "61mer_analysis_reads.txt")
            else:
                annotation_dir = os.path.join(str(amr_annotations), "sample1", "bin1")
                des_path = os.path.join(annotation_dir, "amr_annotation.json")
                src_path = self.get_data_path("rgi_output.json")
                analysis_path = os.path.join(
                    "sample1", "bin1", "61mer_analysis_mags.txt"
                )

            os.makedirs(annotation_dir)
            shutil.copy(src_path, des_path)
            result = query_function(
                amr_annotations=amr_annotations, card_db=card_db, kmer_db=kmer_db
            )
            self.assertIsInstance(result, output_format)
            self.assertTrue(os.path.exists(os.path.join(str(result), analysis_path)))

    def test_kmer_query_mags_card(self):
        self._run_kmer_query_test(
            annotation_format=CARDAnnotationDirectoryFormat,
            output_format=CARDMAGsKmerAnalysisDirectoryFormat,
            query_function=kmer_query_mags_card,
        )

    def test_kmer_query_reads_card(self):
        self._run_kmer_query_test(
            annotation_format=CARDGeneAnnotationDirectoryFormat,
            output_format=CARDReadsKmerAnalysisDirectoryFormat,
            query_function=kmer_query_reads_card,
        )

    def test_run_rgi_kmer_query(self):
        with patch("q2_amr.card.kmer.run_command") as mock_run_command:
            run_rgi_kmer_query(
                tmp="path_tmp",
                input_file="path_input_file",
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
                    "path_input_file",
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

    def copy_kmer_files(
        self, tmp, input_directory, card_fasta, kmer_size, threads, batch_size
    ):
        src_des_list = [
            ("kmer_json_test.json", f"{kmer_size}_kmer_db.json"),
            ("kmer_txt_test.txt", f"all_amr_{kmer_size}mers.txt"),
        ]
        for scr_file, des_file in src_des_list:
            shutil.copy(self.get_data_path(scr_file), os.path.join(tmp, des_file))

    def test_kmer_build_card(self):
        mock_run_rgi_kmer_build = MagicMock(side_effect=self.copy_kmer_files)
        with patch(
            "q2_amr.card.kmer.run_rgi_kmer_build", side_effect=mock_run_rgi_kmer_build
        ), patch("q2_amr.card.kmer.load_card_db"), patch("glob.glob"):
            card_db = CARDDatabaseDirectoryFormat()
            result = kmer_build_card(card_db=card_db, kmer_size=32)
            self.assertIsInstance(result, CARDKmerDatabaseDirectoryFormat)
            for file in ["32_kmer_db.json", "all_amr_32mers.txt"]:
                self.assertTrue(os.path.exists(os.path.join(str(result), file)))

    def test_run_rgi_kmer_build(self):
        with patch("q2_amr.card.kmer.run_command") as mock_run_command:
            run_rgi_kmer_build(
                tmp="path_tmp",
                input_directory="path_directory",
                card_fasta="path_fasta",
                kmer_size="61",
                threads="10",
                batch_size="1000000",
            )
            mock_run_command.assert_called_once_with(
                [
                    "rgi",
                    "kmer_build",
                    "--input_directory",
                    "path_directory",
                    "--card",
                    "path_fasta",
                    "-k",
                    "61",
                    "--threads",
                    "10",
                    "--batch_size",
                    "1000000",
                ],
                "path_tmp",
                verbose=True,
            )
