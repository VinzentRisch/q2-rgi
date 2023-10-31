import subprocess

from q2_amr.card.utils import run_command


def run_rgi_kmer_query(tmp, input_file, input_type, kmer_size, minimum, threads):
    cmd = [
        "rgi",
        "kmer_query",
        "--input",
        input_file,
        f"--{input_type}",
        "--kmer_size",
        kmer_size,
        "--minimum",
        minimum,
        "--threads",
        threads,
        "--output",
        "--local",
    ]

    try:
        run_command(cmd, tmp, verbose=True)
    except subprocess.CalledProcessError as e:
        raise Exception(
            "An error was encountered while running rgi, "
            f"(return code {e.returncode}), please inspect "
            "stdout and stderr to learn more."
        )
