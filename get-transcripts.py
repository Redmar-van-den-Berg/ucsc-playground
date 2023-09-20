#!/usr/bin/env python3

import os
import argparse
import urllib.request
from urllib.error import HTTPError
from filecache import filecache
import json

from dataclasses import dataclass, field
from draw_blocks import make_drawing, draw_regions

# Set up logging
import logging
logger = logging.getLogger(__name__)

logging.basicConfig(level=logging.INFO)

@dataclass
class Region:
    """Class to hold a region"""
    name: str
    chrom: str
    start: int
    end: int
    size: int = field(init=False)

    def __post_init__(self):
        self.size = self.end - self.start

    def __repr__(self):
        return f"{self.name} {self.chrom}:{self.start:,}-{self.end:,} ({self.size}bp)"

header = ['#name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds', 'proteinID', 'alignID']

def fetch(url):
    logger.info(f"Fetching {url}")
    return _fetch(url)

@filecache(24 * 7 * 60 * 60)
def _fetch(url):
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


def draw_genomic_regions(ts, genomic_region):
    """Draw the regions for a transcript"""
    offset = ts["chromStart"]
    scale = 0.1

    blocks = list()
    for regions in genomic_region.values():
        # Make into tuples
        tupl = [(r.start-offset, r.end-offset) for r in regions]
        # Scale
        tupl = [(start*scale, end*scale) for start, end in tupl]
        blocks.append(tupl)

    return draw_regions(genomic_region, offset)
    return make_drawing(60, *blocks)

def parse_transcript(ts):
    exons = exon_regions(ts)
    coding = [coding_region(ts)]
    
    return {
        "exons": exons,
        "coding": coding,
    }
    drawing = draw_regions(ts, exons, [coding])
    print(drawing)

    print('-'*10, "EXONS", '-'*10)
    for e in exons: print(e)

    print('-'*10, "CODING", '-'*10)
    print(coding)
    print('-'*10, "END", '-'*10)
    print()

def coding_region(ts):
    return Region("CDS", ts["chrom"], ts["thickStart"], ts["thickEnd"])


def print_genomic_region(genomic_region):
    for name, regions in genomic_region.items():
        print('-'*10, name, '-'*10)
        for i in regions:
            print(i)

def exon_regions(ts):
    """Parse transcript into usable format"""
    chrom = ts["chrom"]
    chrom_start = ts["chromStart"]
    exon_starts = [chrom_start + x for x in csv_to_int(ts["chromStarts"])]
    exon_sizes = csv_to_int(ts["blockSizes"])
    exon_ends = [start + size for start, size in zip(exon_starts, exon_sizes)]

    # Create a list of exon Region objects
    regions = list()
    for i, range in enumerate(zip(exon_starts, exon_ends), 1):
        # Ensure numbering is consistent with strand
        if ts["strand"] == '-':
            exon_nr = len(exon_starts) - i + 1
        else:
            exon_nr = i
        start, end = range

        regions.append(Region(f"{exon_nr}", chrom, start, end))
    return regions

def fetch_uniprot_track(track, chrom, start, end, uniprot_id):
    url = f"https://api.genome.ucsc.edu/getData/track?genome=hg38;track={track};chrom={chrom};start={start};end={end}"
    data = fetch(url)

    domain_regions = list()
    for domain in data[track]:
        # Ensure that we only include domains for the correct protein
        if domain["uniProtId"] == uniprot_id:
            # Do some sanity checks
            if domain["chromStart"] != domain["thickStart"] or domain["chromEnd"] != domain["thickEnd"]:
                raise NotImplementedError

            domain_regions.append(Region(domain["name"], domain["chrom"], domain["chromStart"], domain["chromEnd"]))
            #print(json.dumps(domain, indent=True))
    return domain_regions

def main(transcript, format):
# Get the location for the specified transcript
    chrom, start, end, version = get_region(transcript)

    # Add the version to the transcript name
    transcript = f"{transcript}.{version}"

    url = f"https://api.genome.ucsc.edu/getData/track?genome=hg38;track=knownGene;chrom={chrom};start={start};end={end}"
    transcript_data = fetch(url)

    for ts in transcript_data["knownGene"]:
        if ts["name"] == transcript:
            genomic_region = parse_transcript(ts)
            break
    else:
        raise RuntimeError(f"transcript {transcript} not found")

    # Next, we get some more tracks we are interested in
    uniprot_tracks = ["unipDomain", "unipStruct", "unipLocCytopl", "unipLocTransMemb", "unipLocExtra", "unipRepeat"]
    uniprot_id = ts["geneName2"]
    for track in uniprot_tracks:
        genomic_region[track] = fetch_uniprot_track(track, chrom, start, end, uniprot_id)
    #print(json.dumps(data, indent=True))

    output(ts, genomic_region, format)

def output(transcript, genomic_region, format):
    if format == "text":
        print_genomic_region(genomic_region)
    elif format=="json":
        print(genomic_region)
    elif format == "svg":
        print(draw_genomic_regions(transcript, genomic_region))
    else:
        raise NotImplementedError


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('transcript')
    parser.add_argument('--format', choices=['text', 'json', 'svg'])

    args = parser.parse_args()

    main(args.transcript, args.format)

