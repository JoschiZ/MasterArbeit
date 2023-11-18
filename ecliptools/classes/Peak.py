from dataclasses import dataclass
from typing import Literal


class Info(object):
    ID: str
    object_type: str
    display_name: str
    biotype: str
    description: str
    seq_region_name: str
    is_canonical: bool
    parent: str | None = None

    def __init__(self, info_dict: dict):
        self.ID = info_dict["ID"]
        self.object_type = info_dict["object_type"]
        self.display_name = info_dict["display_name"]
        self.biotype = info_dict["biotype"]
        self.description = info_dict["description"]
        self.seq_region_name = info_dict["seq_region_name"]
        self.parent = info_dict["parent"]
        self.is_canonical = bool(info_dict["is_canonical"])

    @property
    def base_display_name(self):
        splitted = self.display_name.split("-")
        if len(splitted[-1]) and splitted[-1].isdigit():
            return "-".join(splitted[0:-1])
        return self.display_name

    def to_dict(self):
        return {
            "display_name": self.base_display_name,
            "description": self.description,
        }


@dataclass
class Peak(object):
    chrom: str
    start: int
    end: int
    name: str
    strand: Literal["+", "-"]

    @property
    def unversioned_name(self):
        return self.name.split(".")[0]

    @property
    def position_string(self):
        return f"{self.chrom}:{self.start}-{self.end}"


@dataclass()
class Annotation(Peak):
    type: str
    sub_type: str
    info: Info | None


@dataclass()
class BasePeak(Peak):
    annotations: list[Annotation]
