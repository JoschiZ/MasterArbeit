#!/usr/bin/python
import json
import sys
from pathlib import Path
from typing import TextIO, Hashable

import click
import pandas as pd
import requests

REQUEST_CACHE: dict[Hashable, str] = {}


def get_references(gene_id: Hashable) -> str:
    if pd.isna(gene_id):
        return ""

    server = "https://rest.ensembl.org"
    ext = f"/xrefs/id/{gene_id}?"
    full_url = server+ext

    if gene_id in REQUEST_CACHE:
        return REQUEST_CACHE[gene_id]

    r = requests.get(server + ext, headers={"Content-Type": "application/json"})

    if not r.ok:
        click.echo(f"Request Error for {full_url} {r.status_code}: {r.reason}", file=sys.stderr)
        return "\"[]\""

    decoded = r.json()

    synonym_set = set()
    for ref in decoded:
        if ref["synonyms"]:
            synonym_set.update(ref["synonyms"])

    synonym_set_string = "[\"" + "\",\"".join(synonym_set) + "\"]"
    if len(synonym_set) == 0:
        synonym_set_string = "[]"
    REQUEST_CACHE[gene_id] = synonym_set_string
    return synonym_set_string


def get_ref_df(df: pd.DataFrame) -> pd.DataFrame:
    reference_dict: dict[Hashable, str] = {}

    with click.progressbar(df.iterrows(), length=df.shape[0], file=sys.stderr) as bar:
        for i, row in bar:
            reference_dict[i] = get_references(i)

    return pd.DataFrame.from_dict(reference_dict, orient="index", columns=["aliases"])


def append_references(in_file: TextIO | Path, out_file: TextIO | Path | None, lookup_file: TextIO | Path) -> pd.DataFrame:

    main_df = pd.read_csv(in_file, sep="\t", header=0, index_col=0)

    ref_df = get_ref_df(main_df)
    joined = main_df.join(ref_df)

    if isinstance(lookup_file, Path):
        lookup_file = lookup_file.open("r")

    lookup_json = json.load(lookup_file)
    lookup_df = pd.DataFrame.from_dict(lookup_json, orient="index")
    lookup_df.drop("ID", axis=1, inplace=True)
    joined = joined.join(lookup_df)

    if not out_file:
        return joined

    if isinstance(out_file, Path):
        out_file = out_file.open("w")

    joined.to_csv(out_file, sep="\t", encoding="utf8", lineterminator="\n")
    return joined


@click.command("append-refs")
@click.option("-o", "--out-file", type=click.File("w"), default=sys.stdout)
@click.option("-i", "--in-file", type=click.File("r"), default=sys.stdin)
@click.argument("lookup_file", type=click.File("r"))
def append_references_command(out_file: TextIO, in_file: TextIO, lookup_file: TextIO):
    """
    This will take in annotated table and adds information from a JSON formatted lookup file to it.
    It will also search for different names in the ENSEMBLE database and append them as well.
    :param out_file:
    :param in_file: The base table to enrich with reference data
    :param lookup_file: A JSON file obtained by get-lookup-file
    :return:
    """
    append_references(in_file, out_file, lookup_file)


if __name__ == "__main__":
    test = True
    if test:
        in_path = Path("../../data/counted_new.tsv")
        annot_path = Path("../../data/json/full.lookup.json")
        out_path = Path("../../data/counted.referenced.tsv")
        result = append_references(in_path, out_path, annot_path)
    else:
        append_references_command()
