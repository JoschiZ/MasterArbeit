import math
import sys
from typing import TypedDict, Literal

import click
import requests

from classes.BedPeak import BedPeak


class SequenceResponse(TypedDict):
    query: str
    molecule: str
    version: int
    desc: str
    id: str
    seq: str


SEQUENCE_CACHE: dict[str, SequenceResponse] = {}


def get_ensemble_sequence(ensemble_accession: str, cds=False) -> SequenceResponse:
    if ensemble_accession in SEQUENCE_CACHE:
        return SEQUENCE_CACHE[ensemble_accession]

    server = "https://rest.ensembl.org"
    resource = f"/sequence/id/{ensemble_accession}?"
    request_url = server + resource

    headers = {
        "Content-Type": "application/json",
        "Type": cds.__str__()
    }

    response = requests.get(request_url, headers=headers)

    if not response.ok:
        click.echo(f"Request Error for {request_url} {response.status_code}: {response.reason}", file=sys.stderr)
        return SequenceResponse(query=ensemble_accession, molecule="FAILED", version=-1, desc="", id="", seq="")

    decoded_response = response.json()
    SEQUENCE_CACHE[ensemble_accession] = decoded_response
    return decoded_response


class Position(TypedDict):
    chrom: str
    start: int
    end: int
    strand: Literal[1, -1]


def parse_ensemble_position(desc: str) -> Position:
    """
    example: "chromosome:GRCh38:3:193593144:193697811:1"
    """
    split = desc.split(":")
    chrom = split[2]
    start = int(split[3])
    end = int(split[4])

    strand: Literal[1, -1]
    if split[5] == "1" or split[5] == "+":
        strand = 1
    elif split[5] == "-1" or split[5] == "-":
        strand = -1
    else:
        raise NotImplementedError("Strand has to be 1/+ or -1/-")

    return Position(chrom=chrom, start=start, end=end, strand=strand)


def get_binding_sequence(gene_accession: str, start: int, end: int) -> str:
    seq_response = get_ensemble_sequence(gene_accession)
    response_position = parse_ensemble_position(seq_response["desc"])

    start_offset = start - response_position["start"]
    end_offset = start_offset + end - start

    if start_offset < 0:
        click.echo(f"WARN: Sequence got clipped by {abs(start_offset)} at the start", sys.stderr)
        start_offset = 0
    if end_offset > len(seq_response["seq"]):
        clipped = abs(end_offset - len(seq_response['seq']))
        click.echo(f"WARN: Sequence got clipped by {clipped} at the end", sys.stderr)
        return seq_response["seq"][start_offset:]

    return seq_response["seq"][start_offset:end_offset]


def get_reading_frame(transcript_accession: str, peak_sequence: str):
    transcript_res = get_ensemble_sequence(transcript_accession, cds=True)
    transcript_seq = transcript_res["seq"]

    peak_start = transcript_seq.find(peak_sequence)
    peak_end = peak_start + len(peak_sequence)

    if peak_start == -1:
        return

    codon_counts = {}
    current_pos = 0
    while True:
        codon = transcript_seq[current_pos:current_pos + 3]
        current_pos += 3

        if current_pos -3 <= peak_start:
            continue
        if current_pos > peak_end:
            break

        if codon in codon_counts:
            codon_counts[codon] += 1
        else:
            codon_counts[codon] = 1

    return codon_counts



if __name__ == "__main__":
    peak_sequence = get_binding_sequence("ENSG00000198836", 193614721, 193695408)
    result = get_reading_frame("ENST00000361510", peak_sequence)
