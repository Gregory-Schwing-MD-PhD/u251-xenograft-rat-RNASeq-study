"""Write the U251 accession manifest: every raw read file the work needs, as the public archive holds it.

The archive is ENA's mirror of SRA, which carries GEO's runs: GSE338105 is BioProject PRJNA1492510 (the ten
U251 libraries), and the rat references of the floor analysis are PRJNA627944 and PRJNA238328 (GSE53960).
ENA's file report gives each FASTQ's address, byte count and md5, which the fetch checks every file against.

    python3 build_manifest.py                 # from the snapshots under ena/ (read 2026-11-15)
    python3 build_manifest.py --from-ena      # read ENA's file reports again first

The rat runs kept are only those the floor analysis used (FLOOR_RUNS); the U251 libraries are named through
library_map.csv, the repository's one join of run accessions to library names.
"""
import argparse
import csv
import pathlib
import urllib.request

HERE = pathlib.Path(__file__).resolve().parent
ENA = HERE / "ena"
STUDIES = ("PRJNA1492510", "PRJNA627944", "PRJNA238328")
FIELDS = ("run_accession,study_accession,secondary_study_accession,sample_alias,sample_title,library_layout,"
          "read_count,fastq_bytes,fastq_md5,fastq_ftp,first_public")
REPORT = "https://www.ebi.ac.uk/ena/portal/api/filereport?accession={}&result=read_run&fields={}&format=tsv"
# The rat runs the floor analysis used: the 22 runs whose FASTQs sat in u251_science/spread/floor/fastq.
FLOOR_RUNS = tuple((HERE / "floor_runs.txt").read_text().split())
COLUMNS = ("study", "run", "library", "layout", "file", "url", "bytes", "md5", "read_count", "geo_sample",
           "first_public")


def refresh_reports():
    """Read ENA's file report of every study again into ena/."""
    ENA.mkdir(exist_ok=True)
    for study in STUDIES:
        with urllib.request.urlopen(REPORT.format(study, FIELDS), timeout=120) as answer:
            (ENA / f"ena_{study}.tsv").write_bytes(answer.read())


def library_names():
    """Run accession to the project's library name, from the repository's library_map.csv."""
    with open((HERE / "library_map.csv" if (HERE / "library_map.csv").is_file() else HERE.parent / "library_map.csv"), newline="") as handle:
        return {row["run_accession"]: row["library"] for row in csv.DictReader(handle)}


def wanted(row):
    return row["study_accession"] == "PRJNA1492510" or row["run_accession"] in FLOOR_RUNS


def file_rows(row, names):
    """One manifest row for each FASTQ of one run, with the archive's own bytes and md5."""
    urls = row["fastq_ftp"].split(";")
    sizes = row["fastq_bytes"].split(";")
    sums = row["fastq_md5"].split(";")
    library = names.get(row["run_accession"], row["sample_title"].replace(" ", "_"))
    for url, size, md5 in zip(urls, sizes, sums):
        yield (row["study_accession"], row["run_accession"], library, row["library_layout"],
               url.rsplit("/", 1)[1], "https://" + url, size, md5, row["read_count"], row["sample_alias"],
               row["first_public"])


def manifest_rows():
    names = library_names()
    for study in STUDIES:
        with open(ENA / f"ena_{study}.tsv", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if wanted(row):
                    yield from file_rows(row, names)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--from-ena", action="store_true", help="read ENA's file reports again first")
    parser.add_argument("--out", default=str(HERE / "manifest.tsv"))
    args = parser.parse_args()
    if args.from_ena:
        refresh_reports()
    rows = sorted(manifest_rows())
    with open(args.out, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(COLUMNS)
        writer.writerows(rows)
    print(f"{len(rows)} files, {sum(int(r[6]) for r in rows)} bytes -> {args.out}")


if __name__ == "__main__":
    main()
