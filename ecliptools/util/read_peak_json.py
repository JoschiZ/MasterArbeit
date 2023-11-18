#!/usr/bin/python
import jsonpickle

from ecliptools.classes.Peak import BasePeak, Peak, Annotation, Info


def read_peak_json(json_text: str) -> list[BasePeak]:
    return jsonpickle.decode(json_text, classes=[Peak, BasePeak, Annotation, Info])


if __name__ == "__main__":
    with open("../../data/json/ENCFF128AKC.mapped.json", "r") as f:
        result = read_peak_json(f.read())
