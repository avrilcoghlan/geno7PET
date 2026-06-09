# Geno7PET

A tool for classifying _Vibrio cholerae_ genome assembly genotypes.

## About

## Example output

## Dependencies

Geno7PET can be run using [python](#pythonpip), [uv](#uv), [pixi](#pixi) or [Docker](#Docker). It also requires BLAST to
be installed, unless you're using [pixi](#pixi) or [Docker](#Docker).

## Installation and running

First download this repository, and then either install as a package using [python](#pythonpip) or [uv](#uv), or use it
in a conda environment with [pixi](#pixi).
You can also run it directly using [python](#pythonpip), [uv](#uv), [pixi](#pixi) or [Docker](#Docker).

### uv

It's simple to run directly. Note that when running from the home directory of the repository, the --project option
should be excluded.

```bash
uv run --project /path/to/genopet/repository/ Geno7PET </path/to/your.fasta> <out_file.txt>
```

You can also install it into your system and run it anywhere.

```bash
uv pip install --system .
cd /some/other/dir/
Geno7PET -h
```

### Python/pip

In the Geno7PET home directory, either run it directly or install it and run anywhere

```bash
# Run directly
python -m geno7PET -h

# Install
pip install .
cd /to/another/location
Geno7PET -h
```

### Pixi

Pixi will automatically install BLAST. In the Geno7PET home directory:

```bash
pixi run Geno7PET -h
```

Or as conda environment:
```bash
pixi install
pixi shell
Geno7PET -h
```

### Docker

A [Dockerfile](Dockerfile) is provided for building a Docker image of Geno7PET.

```bash
docker build --pull --rm -t geno7pet .
docker run --rm -i -v "${PWD}":/workdir geno7pet:latest /workdir/my.fasta /workdir/out.txt
```
