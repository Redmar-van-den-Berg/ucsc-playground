#!/usr/bin/env python3

import os
import argparse
import urllib.request
from urllib.error import HTTPError
from filecache import filecache
import json

from dataclasses import dataclass, field

@dataclass
class Region:
    """Class to hold a region"""
    chrom: str
    start: int
    end: int
    size: int = field(init=False)

    def __post_init__(self):
        self.size = self.end - self.start

    def __repr__(self):
        return f"{self.chrom}:{self.start}-{self.end} ({self.size}bp)"

header = ['#name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds', 'proteinID', 'alignID']

@filecache(24 * 7 * 60 * 60)
def fetch(url):
    try:
        response = urllib.request.urlopen(url)
    except HTTPError as e:
        msg = f"Fetching '{url}' returned {e}"
        raise RuntimeError(msg)
    return json.loads(response.read())

def fetch_gene(gene, genome):
    url = f"https://api.genome.ucsc.edu/search?search={gene}&genome={genome}"
    return fetch(url)


def fetch_files(genome):
    url = f"https://api.genome.ucsc.edu/list/files?genome={genome}"
    return fetch(url)


def find_track(data, trackname):
    for track in data["positionMatches"]:
        if track["trackName"] == trackname:
            return track

def csv_to_int(string):
    return list(map(int, (x for x in string.split(',') if x)))

def get_region(transcript):
    """ Return the region for an Ensembl transcript"""
    url = f"https://rest.ensembl.org/lookup/id/{transcript}?content-type=application/json"
    js = fetch(url)

    if js["assembly_name"] != "GRCh38":
        msg = f"Unsupported reference {js['assembly_name']}"
        raise NotImplementedError(msg)

    chrom =  js["seq_region_name"]
    start = js["start"]
    end = js["end"]
    version = js["version"]
    return chrom, start, end, version

def uscs_to_tsv(ts):
    chromStarts = csv_to_int(ts["chromStarts"])
    blockSizes = csv_to_int(ts["blockSizes"])
    exonStarts = [ts["chromStart"] + offset for offset in chromStarts]
    exonEnds = [start + size for start, size in zip(exonStarts, blockSizes)]
    values = {
        "#name": ts["name"],
        "chrom": ts["chrom"],
        "strand": ts["strand"],
        "txStart": ts["chromStart"],
        "txEnd": ts["chromStart"] + chromStarts[-1] + blockSizes[-1],
        "cdsStart": ts["thickStart"],
        "cdsEnd": ts["thickEnd"],
        "exonCount": ts["blockCount"],
        "exonStarts": ','.join((str(x) for x in exonStarts)) + ',',
        "exonEnds": ','.join((str(x) for x in exonEnds)) + ',',
        "proteinID": "" if ts.get("geneName2") == "none" else ts["geneName2"],
        "alignID": ts["name2"],
    }
    print(*(values[key] for key in header), sep='\t')


def parse_transcript(ts):
    exons = exon_regions(ts)
    coding = coding_region(ts)

    print('-'*10, "EXONS", '-'*10)
    for e in exons: print(e)

    print('-'*10, "CODING", '-'*10)
    print(coding)
    print('-'*10, "END", '-'*10)
    print()

def coding_region(ts):
    return Region(ts["chrom"], ts["thickStart"], ts["thickEnd"])

def exon_regions(ts):
    """Parse transcript into usable format"""
    chrom = ts["chrom"]
    chrom_start = ts["chromStart"]
    exon_starts = [chrom_start + x for x in csv_to_int(ts["chromStarts"])]
    exon_sizes = csv_to_int(ts["blockSizes"])
    exon_ends = [start + size for start, size in zip(exon_starts, exon_sizes)]

    # Create a list of exon Region objects
    return [Region(chrom, start, end) for start, end in zip(exon_starts, exon_ends)]
    for region in exons:
        print(region)

    print(json.dumps(ts, indent=True))
    exit()
    pass

def main(transcript):
    chrom, start, end, version = get_region(transcript)

    # Add the version to the transcript name
    transcript = f"{transcript}.{version}"

    # DMD
    url = f"https://api.genome.ucsc.edu/getData/track?genome=hg38;track=knownGene;chrom={chrom};start={start};end={end}"
    data = fetch(url)

    for ts in data["knownGene"]:
        if ts["name"] == transcript:
            parse_transcript(ts)
            #uscs_to_tsv(ts)
            print(json.dumps(ts, indent=True))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('transcript')

    args = parser.parse_args()

    main(args.transcript)

