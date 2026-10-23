#!/usr/bin/env python3
r"""Rerun only the GISAID naming for samples that the pipeline already processed.

Usage::

    python gisaid_naming_rerun.py --run-dir <run_dir> --config <config.json> \
        (--sample <sample_name> | --all) [--force]

Examples::

    python gisaid_naming_rerun.py --run-dir /path/to/RESULTS_20261016 \
        --config /path/to/config.json --sample 26-040001
    python gisaid_naming_rerun.py --run-dir /path/to/RESULTS_20261016 \
        --config /path/to/config.json --all

``--all`` skips samples that already have a GISAID name
in the central list, so it only names the ones still missing a name.
``--sample`` names the chosen sample even if it has a name, keeping its
running number (e.g. to recreate its files); ``--force`` gives it a new one. MinION samples (``barcode01``, ...) are never GISAID-named: ``--all``
skips them and ``--sample barcodeNN`` stops with an error, so barcode names
never end up in the central name list.
"""
import argparse
import json
import logging
from collections.abc import Iterator
from pathlib import Path

from gisaid_naming import process_sample_folder

logger = logging.getLogger(__name__)


def is_minion_barcode(sample_name: str) -> bool:
    """Check whether a sample is a MinION barcode (named after its barcode folder).

    Args:
        sample_name: Sample name.

    Returns:
        ``True`` for names such as ``barcode01`` (case-insensitive).
    """
    return sample_name.lower().startswith("barcode")


def find_irma_sample_dirs(run_dir: Path) -> Iterator[tuple[str, Path]]:
    """Find the samples in a run folder that IRMA produced a consensus for.

    Looks for ``<run_dir>/<sample>/irma_results/<sample>/amended_consensus``.
    MinION barcode samples are skipped.

    Args:
        run_dir: Results folder of the run.

    Yields:
        ``(sample_name, irma_sample_path)`` for each sample, in name order.
    """
    for sample_dir in sorted(run_dir.iterdir()):
        if not sample_dir.is_dir():
            continue
        if is_minion_barcode(sample_dir.name):
            logger.info("Ohitetaan MinION-näyte %s (ei GISAID-nimeä)", sample_dir.name)
            continue
        irma_sample_path = sample_dir / "irma_results" / sample_dir.name
        if (irma_sample_path / "amended_consensus").exists():
            yield sample_dir.name, irma_sample_path


def main() -> None:
    """Parse the arguments and GISAID-name the chosen samples.

    Raises:
        SystemExit: If the run folder, ``gisaid_naming_file`` in the config,
            or the chosen sample's consensus is missing, or a MinION sample
            is chosen.
    """
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser(description="Rerun only GISAID naming for sample(s)")
    parser.add_argument("--run-dir", required=True, type=Path,
                        help="Ajon tuloskansio, esim. .../RESULTS_20261016")
    parser.add_argument("--config", required=True, type=Path,
                        help="config.json, josta luetaan gisaid_naming_file")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--sample", help="Yksittäisen näytteen nimi (kansion nimi run_dir:in alla)")
    group.add_argument("--all", action="store_true", help="Aja kaikille näytteille run_dir:in alla")
    parser.add_argument("--force", action="store_true",
                        help="Ohittaa aiemman juoksevan numeron ja antaa uuden")
    args = parser.parse_args()

    run_dir = args.run_dir.resolve()
    if not run_dir.exists():
        raise SystemExit(f"run_dir ei ole olemassa: {run_dir}")

    with open(args.config) as f:
        config = json.load(f)

    gisaid_excel_str = config.get("gisaid_naming_file")
    if not gisaid_excel_str:
        raise SystemExit("config.json: 'gisaid_naming_file' puuttuu")
    gisaid_excel = Path(gisaid_excel_str)

    if args.all:
        targets = list(find_irma_sample_dirs(run_dir))
        if not targets:
            raise SystemExit(f"Yhtään GISAID-nimettävää amended_consensus-näytettä ei löytynyt kansiosta {run_dir}")
    else:
        if is_minion_barcode(args.sample):
            raise SystemExit(f"{args.sample} on MinION-näyte; MinION-näytteitä ei GISAID-nimetä.")
        irma_sample_path = run_dir / args.sample / "irma_results" / args.sample
        if not (irma_sample_path / "amended_consensus").exists():
            raise SystemExit(f"amended_consensus-kansiota ei löydy: {irma_sample_path}")
        targets = [(args.sample, irma_sample_path)]

    for sample_name, irma_sample_path in targets:
        logger.info("Nimetään näyte: %s", sample_name)
        try:
            new_name = process_sample_folder(
                sample_path=irma_sample_path,
                excel_file=gisaid_excel,
                run_dir=run_dir,
                force_reassign=args.force,
                rename_existing=bool(args.sample),
            )
        except RuntimeError:
            logger.exception("GISAID-nimeäminen epäonnistui näytteelle %s", sample_name)
            continue

        if new_name:
            logger.info("-> %s : %s", sample_name, new_name)
        else:
            logger.warning("-> %s : ohitettiin (jo nimetty, puuttuu HA/NA tms.; syy yllä)", sample_name)


if __name__ == "__main__":
    main()