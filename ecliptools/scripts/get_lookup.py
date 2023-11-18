#!/usr/bin/python
import json
import sys
from pathlib import Path
from typing import TextIO

import click
import regex
import requests


PARENT_CACHE = {}


def get_lookup_data(ensembl_id: str, deep=True) -> dict:
    server = "https://rest.ensembl.org"
    head = ensembl_id.split(".", 1)[0]

    ext = f"/lookup/id/{head}?"

    r = requests.get(server + ext, headers={"Content-Type": "application/json", "species": "homo_sapiens"})

    if not r.ok:
        stderr = click.get_text_stream("stderr")
        stderr.write(f"Request Error for {server + ext} {r.status_code}: {r.reason}")
        return {}

    decoded = r.json()
    description = ""

    if "description" in decoded:
        description = decoded["description"]

    display_name = ""
    if "display_name" in decoded:
        display_name = decoded["display_name"]

    biotype = ""
    if "biotype" in decoded:
        biotype = decoded["biotype"]

    object_type = ""
    if "object_type" in decoded:
        object_type = decoded["object_type"]

    seq_region_name = ""
    if "seq_region_name" in decoded:
        seq_region_name = decoded["seq_region_name"]

    parentID = ""
    if "Parent" in decoded:
        parentID = decoded["Parent"]

    is_canonical = ""
    if "is_canonical" in decoded:
        is_canonical = decoded["is_canonical"]

    if deep and "Parent" in decoded:
        if decoded["Parent"] in PARENT_CACHE:
            parent = PARENT_CACHE[decoded["Parent"]]
        else:
            parent = get_lookup_data(decoded["Parent"], False)
            PARENT_CACHE[decoded["Parent"]] = parent

        if description == "" and "description" in parent:
            description = parent["description"]

        if display_name == "" and "display_name" in parent:
            display_name = parent["display_name"]

    return {
        "ID": ensembl_id,
        "object_type": object_type,
        "display_name": display_name,
        "biotype": biotype,
        "description": description,
        "seq_region_name": seq_region_name,
        "parent": parentID,
        "is_canonical": is_canonical
    }


@click.command("get-lookup-file")
@click.option("-o", "--out-file", type=click.File("w"), default=sys.stdout)
@click.argument("in_files", type=click.File("r"), nargs=-1)
def get_lookup_file(out_file: TextIO, in_files: tuple[TextIO]):
    """
    Takes in a arbitrary number of files and extracts all ENSENMBLE identifiers from them.
    Using those it calls the ENSEMBLE /lookup/id/ endpoint.
    The returned data is then written to JSON with the Accession as key.
    :param out_file: Defaults to stdout
    :param in_files:
    :return:
    """
    accession_set = set()
    for file in in_files:
        text = file.read()
        all_accessions = regex.findall("ENS.\d{11}", text)
        accession_set.update(all_accessions)

    accession_lookup = {}
    with click.progressbar(accession_set) as bar:
        for accession in bar:
            accession_lookup[accession] = get_lookup_data(accession)

    accession_lookup |= PARENT_CACHE
    json.dump(accession_lookup, out_file, indent=4)


if __name__ == "__main__":
    get_lookup_file()
