"""Write an nf-core sample sheet whose inputs are the files fetch_from_archive.sh brought back, never old local paths.

    python3 samplesheet_from_manifest.py FASTQ_DIR [--study PRJNA1492510] [--out samplesheet.csv]

Each run becomes one row named by its library (library_map.csv's name for the U251 libraries, the run accession
for the public rat runs), with the paired files in mate order; a file the fetch has not landed stops the sheet,
so a sheet never names a file that is not there.
"""
import argparse
import csv
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent


def runs(manifest, study):
    """Run accession to (library, [files in mate order]) for one study."""
    found = {}
    with open(manifest, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["study"] == study:
                found.setdefault(row["run"], [row["library"], []])[1].append(row["file"])
    return found


def sheet_row(library, files, fastq_dir):
    paths = [fastq_dir / name for name in sorted(files)]
    missing = [str(path) for path in paths if not path.is_file()]
    if missing:
        sys.exit("not fetched yet: " + ", ".join(missing))
    second = str(paths[1]) if len(paths) > 1 else ""
    return [library, str(paths[0]), second, "auto"]


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("fastq_dir", type=pathlib.Path)
    parser.add_argument("--study", default="PRJNA1492510")
    parser.add_argument("--manifest", default=str(HERE / "manifest.tsv"))
    parser.add_argument("--out", default="samplesheet.csv")
    args = parser.parse_args()
    fastq_dir = args.fastq_dir.resolve()
    rows = [sheet_row(lib, files, fastq_dir) for _, (lib, files) in sorted(runs(args.manifest, args.study).items())]
    with open(args.out, "w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["sample", "fastq_1", "fastq_2", "strandedness"])
        writer.writerows(sorted(rows))
    print(f"{len(rows)} samples -> {args.out}")


if __name__ == "__main__":
    main()
