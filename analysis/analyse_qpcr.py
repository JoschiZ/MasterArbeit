import math
from pathlib import Path
import pandas as pd


def read_qpcr_result(path: Path):
    df = pd.read_csv(path, encoding="windows-1252", skiprows=14)
    df = df[df["Omit"] != True]
    df = df[df["Amp Status"].str.contains("No Amp") == False]
    df["Ct"] = df["Ct"].astype("float")
    return df


def get_average_by_repeat(df: pd.DataFrame):
    grouped = df.groupby(["Experiment Name", "Target Name", "Sample Name"])["Ct"]
    sum = grouped.sum()
    return sum.div(grouped.count(), axis=0)


def calculate_delta_ct(df: pd.DataFrame):
    dict = {}
    for k, v in df.items():
        gapdh_key = list(k)
        gapdh_key[1] = "GAPDH"
        gapdh_key = tuple(gapdh_key)

        row = list(k)
        row.append(v - df[gapdh_key])

        dict[len(dict)] = row

    new_df = pd.DataFrame.from_dict(dict, orient="index", columns=["Experiment Name", "Target Name", "Sample Name", "deltaCt"])
    return new_df


def calculate_delta_delta_ct(df: pd.DataFrame):
    wt = df[df["Sample Name"] == "WT"]
    kd = df[df["Sample Name"] == "KD"]

    wt_ct = wt.set_index(["Experiment Name", "Target Name"])["deltaCt"]
    kd_ct = kd.set_index(["Experiment Name", "Target Name"])["deltaCt"]

    delta_delta_ct = kd_ct - wt_ct
    delta_delta_ct = delta_delta_ct.reset_index()
    delta_delta_ct["Sample Name"] = "KD"
    delta_delta_ct.set_index(["Experiment Name", "Target Name", "Sample Name"], inplace=True)
    delta_delta_ct.rename(columns={"deltaCt": "deltaDeltaCt"}, inplace=True)

    df.set_index(["Experiment Name", "Target Name", "Sample Name"], inplace=True)

    result = pd.concat([df, delta_delta_ct], axis=1)
    result.fillna(0, inplace=True)
    return result


def calculate_rq(df: pd.DataFrame):
    df["RQ"] = 2**(-df["deltaDeltaCt"])
    return df


if __name__ == "__main__":
    data_path = Path("EXPRESSION_SUITE_RESULT_EXPORT.csv")
    df = read_qpcr_result(data_path)
    repeat_average = get_average_by_repeat(df)
    delta_ct = calculate_delta_ct(repeat_average)
    delta_delta_ct = calculate_delta_delta_ct(delta_ct)
    rq = calculate_rq(delta_delta_ct)

    rq.reset_index(inplace=True)
    rq.drop("Experiment Name", axis=1, inplace=True)
    grouped = rq.groupby(["Sample Name", "Target Name"])

    sum = grouped.sum()
    mean = sum.div(grouped.count(), axis=0)

    description = grouped.describe()
    description.to_excel(Path("DESCRIPTION.xlsx"))

    delta_ct.reset_index(inplace=True)
    delta_ct.sort_values(["Target Name", "Sample Name"], inplace=True)
    delta_ct.to_excel("RESULTS.xlsx", index=False)

