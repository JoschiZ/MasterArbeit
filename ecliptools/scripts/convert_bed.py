#!/usr/bin/python
import json
from pathlib import Path
from typing import TextIO
import click
import jsonpickle
import pandas as pd
from ecliptools.classes.Peak import Annotation, Info, BasePeak
from ecliptools.util.create_dynamic_column_names import create_dynamic_column_names


@click.group()
def convert():
    pass


def convert_df_to_json(df: pd.DataFrame, enst_lookup_json: dict[str, dict[str, any]]):
    out = []
    for _, row in df.iterrows():
        annotations: list[Annotation] = []

        if row[10] == "UNKNOWN":
            continue
        for i, column in row.items():
            if pd.isna(column):
                continue

            column_num = int(i)
            if column_num <= 8 or (column_num - 10) % 7 != 0:
                continue

            if row[column_num+3] != row[5]:
                continue

            info: Info | None = None
            unversioned_name = row[column_num+4].split(".")[0]
            if unversioned_name in enst_lookup_json:
                info_dict = enst_lookup_json[unversioned_name]
                info = Info(info_dict=info_dict)

            annotation = Annotation(
                chrom=row[column_num],
                start=int(row[column_num+1]),
                end=int(row[column_num+2]),
                strand=row[column_num+3],
                name=row[column_num+4],
                type=row[column_num+5],
                sub_type=row[column_num+6],
                info=info
            )
            annotations.append(annotation)

        base_peak = BasePeak(
            chrom=row[0],
            start=row[1],
            end=row[2],
            name=row[3],
            strand=row[5],
            annotations=annotations
        )

        out.append(base_peak)

    return out


def tojson(in_path: Path, lookup_file: TextIO | Path, out_file: TextIO | Path, clear_text=False):
    """
    Takes in a annotated bed file and returns a json file. Appends further information from a lookup file.
    Remove annotations that did not respect strandedness
    :param in_path: Path to the annotated bed file that is to be converted
    :param lookup_file: A JSON file containing additional information created by get_lookup_file
    :param out_file:
    :param clear_text:
    :return:
    """
    if isinstance(lookup_file, Path):
        lookup_file = lookup_file.open("r")

    if isinstance(out_file, Path):
        out_file = out_file.open("w")

    df = pd.read_csv(in_path, sep="\t", header=None, index_col=None, names=create_dynamic_column_names(in_path),
                     engine="python")

    lookup = json.load(lookup_file)
    json_df = convert_df_to_json(df, lookup)
    out_json = jsonpickle.encode(json_df, unpicklable=not clear_text, indent=4, include_properties=True)
    out_file.write(out_json)


@convert.command("to-json")
@click.argument("in_path", type=click.Path(exists=True, path_type=Path))
@click.argument("lookup_file", type=click.File("r"))
@click.argument("out_file", type=click.File("w"))
@click.option("--clear/--no-clear", default=False,
              help="If True this will output humanly readable json. False will allow further use with this tools.")
def tojson_command(in_path: Path, lookup_file: TextIO, out_file: TextIO, clear=False):
    """
    Takes in one or more annotated bed files and returns a json file.
    Appends further information from a lookup file.
    Remove annotations that did not respect strandedness.
    :param in_path: Path to the annotated bed file that is to be converted
    :param lookup_file: A JSON file containing additional information created by get_lookup_file
    :param out_file:
    :param clear: If True this outputs more humanly readable json, but further processing with ecliptools is not possible
    :return:
    """
    tojson(in_path, lookup_file, out_file, clear)


if __name__ == "__main__":
    tojson(Path("../../data/ENCFF128AKC.mapped.bed"), Path("../../data/json/full.lookup.json"), Path("../../data/json/ENCFF128AKC.mapped.json"))




