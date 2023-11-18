from __future__ import annotations

import os
from dataclasses import dataclass
from glob import iglob
from os import listdir
from os.path import isfile, join
from pathlib import Path
from typing import TextIO

import numpy as np
import openpyxl
import pandas as pd

"""
This script is not a dropin solution and requires you to follow specific naming conventions
"""

@dataclass()
class QuantificationTable:
    df: pd.DataFrame
    repeat: str
    antibody_number: str
    antibody: str

    def __init__(self, file: Path | TextIO, index: None | list[str] = None):
        df = pd.read_excel(file, header=1, skiprows=0)
        wb = openpyxl.load_workbook(file)
        sheet = wb.active
        name_cell = sheet.cell(1, 1).value
        name = name_cell.split(" ")[-1]
        details = name.split("_")
        if len(details) == 3:
            repeat, antibody_number, antibody = details
        elif len(details) == 4:
            repeat, antibody_number, _, antibody = details
        else:
            repeat, antibody_number, _, antibody, _ = details

        print(file)
        print(antibody)

        if repeat == "IP1":  # The sample sequence in this run was different
            index = ["WT", "GFP", "H562R", "R785C", "R493H"]
        if index is None:
            index = ["WT", "GFP", "R493H", "H562R", "R785C"]
        df = df.set_index(pd.Index(index))

        self.df = df
        self.repeat = repeat
        self.antibody_number = antibody_number
        self.antibody = antibody
        print("Created: " + name)

    @property
    def adj_total_band(self):
        adj_total_band = self.df["Adj. Total Band Vol. (Int)"]
        adj_total_band.name = self.repeat + "_" + self.antibody
        return adj_total_band

    @property
    def normalised_to_overall(self) -> pd.Series:
        except_gfp_control = self.adj_total_band[self.adj_total_band.index != "GFP"]
        except_gfp_control = except_gfp_control / except_gfp_control.sum()
        return except_gfp_control

    def get_ip_by_ip(self, dhx_30_quant: QuantificationTable):
        per_dhx30_gfp = round(
            self.normalised_to_overall /
            dhx_30_quant.normalised_to_overall
            , 2)
        per_dhx30_gfp.name = self.repeat + "_" + self.antibody
        return per_dhx30_gfp

    def get_relative_quantification(self, dhx_30_quant: QuantificationTable):
        per_gfp = self.get_ip_by_ip(dhx_30_quant)
        rel_to_wt = round(per_gfp / per_gfp["WT"], 2)
        rel_to_wt.name = self.repeat + "_" + self.antibody
        return rel_to_wt


def read_all_tables(files: Path | list[Path] | list[TextIO]):
    if isinstance(files, Path):
        files = [Path(join(files, f)) for f in listdir(files) if isfile(join(files, f))]

    all_quant_table = []
    for file in files:

        if "Quanti" in file.name or "~$" in file.name:
            continue
        print(file)
        all_quant_table.append(QuantificationTable(file))

    return all_quant_table


dhx_30_cache: dict[str, QuantificationTable] = {}


def get_dhx30_quantification(quant_tables: list[QuantificationTable], repeat: str) -> QuantificationTable:
    """
    Finds the DHX30 measurement for a given repeat in a list of Quantification Tables.
    Results are cached for fast access
    :param quant_tables:
    :param repeat:
    :return:
    """
    if repeat in dhx_30_cache:
        return dhx_30_cache[repeat]

    for quant_table in quant_tables:
        if "GFP" in quant_table.antibody and quant_table.repeat == repeat:
            dhx_30_cache[repeat] = quant_table
            return quant_table


if __name__ == "__main__":

    rootdir_glob = r'PATH_TO_FOLDER_CONTAINING_QUANTIFIED_EXCELS'
    # This will return absolute paths
    dir_list = [f for f in iglob(rootdir_glob, recursive=True) if os.path.isdir(f)]
    for dir in dir_list:
        try:
            result = read_all_tables(list(Path(dir).rglob("*.xlsx")))
            all_relative = []
            all_ip_by_ip = []
            all_normalised_to_overall = []
            for quant in result:
                dhx30 = get_dhx30_quantification(result, quant.repeat)
                all_ip_by_ip.append(quant.get_ip_by_ip(dhx30))
                all_relative.append(quant.get_relative_quantification(dhx30))
                all_normalised_to_overall.append(quant.normalised_to_overall)

            normed = pd.DataFrame(all_normalised_to_overall).transpose().fillna(0).replace([np.inf, -np.inf], 0)
            ip_by_ip = pd.DataFrame(all_ip_by_ip).transpose().fillna(0).replace([np.inf, -np.inf], 0)
            end_result = pd.DataFrame(all_relative).transpose().fillna(0).replace([np.inf, -np.inf], 0)

            ip_number = Path(dir).name
            ex_writer = pd.ExcelWriter(Path(f"data/result_{ip_number}.xlsx"))
            end_result.to_excel(ex_writer, sheet_name="NormedToWT")
            ip_by_ip.to_excel(ex_writer, sheet_name="IPByIP")
            normed.to_excel(ex_writer, sheet_name="NormedToOverall")
            ex_writer.close()
        except Exception as e:
            print(e)
            continue
