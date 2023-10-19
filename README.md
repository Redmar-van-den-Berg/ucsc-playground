# UCSC playground
Some code to play around with the UCSC API

## Usage
```
usage: get-transcripts.py [-h] [--format {text,json,svg}] transcript

positional arguments:
  transcript

options:
  -h, --help            show this help message and exit
  --format {text,json,svg}
```

## Examples
### DMD
```bash
python3 get-transcripts.py ENST00000357033 --format text
```

### FLT3
```bash
python3 get-transcripts.py ENST00000241453.12 --format text
```
