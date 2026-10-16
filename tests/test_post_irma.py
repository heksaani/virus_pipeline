"""CHANGE (2026-10-06): NEW. Post-IRMA parsing on IRMA-shaped output."""
from pathlib import Path

import pandas as pd

from post_irma import (
    _append_df_to_tsv_locked,
    get_flu_n_percentages,
    get_genome_n_percent,
)


def test_flu_n_percent_uses_numbered_files_and_ignores_padding(tmp_path: Path):
    (tmp_path / "s_4.fa").write_text(">s_4\nACGTNACGTN\n")              # 20% N
    (tmp_path / "s_4.pad.fa").write_text(">s_4\nNNNNNNNNNNACGTNACGTN\n")  # padded copy
    result = get_flu_n_percentages(tmp_path)
    assert result == {"HA": (2, 20.0)}


def test_genome_n_percent_empty_sequence_does_not_crash(tmp_path: Path):
    (tmp_path / "s.fa").write_text(">s\n")
    assert get_genome_n_percent(tmp_path) == (None, None)


def test_nextclade_tsv_keeps_na_label_and_commas(tmp_path: Path):
    dest = tmp_path / "nextclade_results.tsv"
    _append_df_to_tsv_locked(pd.DataFrame({"sample_name": ["s1"], "Segmentti": ["NA"],
                                           "substitutions": ["A12G,C34T"]}), dest)
    _append_df_to_tsv_locked(pd.DataFrame({"sample_name": ["s2"], "Segmentti": ["HA"],
                                           "substitutions": ["G5A;x"]}), dest)
    df = pd.read_csv(dest, sep="\t", dtype=str, keep_default_na=False)
    assert df["Segmentti"].tolist() == ["NA", "HA"]
    assert df["substitutions"].tolist() == ["A12G,C34T", "G5A;x"]
