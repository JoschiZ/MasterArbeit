#!/usr/bin/python
import collections.abc
import os
from pathlib import Path
from typing import TextIO, TypedDict, Literal
import click
import pandas as pd
from ecliptools.util.read_peak_json import read_peak_json


class Gene(TypedDict):
    peaks: set[str]
    strand: Literal["+", "-"]


def read_gene_peaks(in_file: TextIO | Path):
    if isinstance(in_file, Path):
        in_file = in_file.open("r")

    json_text = in_file.read()
    peaks = read_peak_json(json_text)

    genes_peaks: dict[str, Gene] = {}
    for peak in peaks:
        for annotation in peak.annotations:
            accession = annotation.name
            if annotation.info and annotation.info.parent:
                accession = annotation.info.parent

            if accession in genes_peaks:
                genes_peaks[accession]["peaks"].add(peak.position_string)
                genes_peaks[accession]["strand"] = annotation.strand
            else:
                genes_peaks[accession] = {"peaks": set()}
                genes_peaks[accession]["peaks"] = {peak.position_string}
                genes_peaks[accession]["strand"] = annotation.strand

    flattened: dict[str, dict] = {}
    for accession, gene_peak in genes_peaks.items():
        peaks_string = "{\"" + "\",\"".join(gene_peak["peaks"]) + "\"}"
        flattened[accession] = {"peaks": peaks_string, "strand": gene_peak["strand"]}

    df = pd.DataFrame.from_dict(flattened, orient="index")
    df["peaks"] = df["peaks"].apply(eval)

    df.rename(columns={"peaks": f"peaks-{os.path.basename(in_file.name)}"}, inplace=True)
    return df


def create_peak_union(*input_files: TextIO | Path):
    all_dfs = []
    if len(input_files) < 1:
        raise ValueError("You need to provide at least one file")

    for input_file in input_files:
        df = read_gene_peaks(input_file)
        all_dfs.append(df)

    df = all_dfs.pop()
    for left_df in all_dfs:
        df = df.merge(left_df.iloc[:, 0], how="outer", left_index=True, right_index=True, copy=False)
        df = df.combine_first(left_df)

    df["peak-union"] = pd.Series(dtype=object)

    for i, row in df.iterrows():
        overall_union = set()
        for column, value in row.items():
            if "peaks-" in column and not pd.isna(value):
                overall_union = overall_union | value
        df["peak-union"][i] = overall_union

    return df


def safe_df(df: pd.DataFrame, out_file: TextIO):
    df.index.name = "ID"
    df.to_csv(out_file, sep="\t", lineterminator="\n")


def save_len(iterable: collections.abc.Sized):
    if not pd.isna(iterable):
        return len(iterable)
    return 0


def count_gene_peaks(df: pd.DataFrame, out_file: TextIO | Path | None = None):
    count_columns: list[str] = []
    union_count_column: str | None = None
    for column in df.columns:
        if "peak" not in column:
            continue
        new_column_name = column.replace("peak", "count")
        if "union" not in new_column_name:
            count_columns.append(new_column_name)
        else:
            union_count_column = new_column_name

        df[new_column_name] = df[column].apply(save_len)

    df["union-max-diff"] = df[union_count_column] - (df[count_columns].max(axis=1) - df[count_columns].min(axis=1))

    df.sort_values(by=["union-max-diff", "count-union"], inplace=True, ascending=False)

    # early exist if it should not be saved to file
    if not out_file:
        return df

    if isinstance(out_file, Path):
        out_file = out_file.open("w")

    safe_df(df, out_file)
    return df


@click.command("count-gene-peaks")
@click.argument("out_file", type=click.File("w"))
@click.argument("in_files", type=click.File("r"), nargs=-1)
def count_gene_peaks_command(out_file: TextIO, in_files: TextIO):
    union = create_peak_union(*in_files)
    count_gene_peaks(union, out_file)


if __name__ == "__main__":
    test = True
    if test:
        path1 = Path("../../data/json/ENCFF663QIZ.mapped.json")
        path2 = Path("../../data/json/ENCFF128AKC.mapped.json")
        out = Path("../../data/counted_new.tsv")
        union = create_peak_union(path1, path2)
        counted = count_gene_peaks(union, out)
    else:
        count_gene_peaks_command()
