r"""Watch an input folder and run the IRMA pipeline when new FASTQ files arrive.

Each new FASTQ file restarts a quiet-period timer (the debounce). When no new
file has appeared for that long, ``run_pipeline.py`` is started on the whole
input folder. If files arrive while a run is in progress, one more run is
started as soon as it finishes.

Usage::

    python watch_and_run.py -i /path/to/input -o /path/to/output \
        --config /path/to/config.json [--debounce 55] \
        [--pipeline-script /path/to/run_pipeline.py] [--dry-run]

Options:
    -i, --input          Folder to watch for FASTQ files (searched recursively).
    -o, --output         Base output folder for pipeline results.
    --config             JSON config file passed on to the pipeline.
    --debounce           Seconds to wait after the last new FASTQ before
                         starting a run (default 55).
    --pipeline-script    Pipeline script to run (default: run_pipeline.py
                         next to this file).
    --dry-run            Pass --dry-run on to the pipeline.
"""
from __future__ import annotations

import argparse
import logging
import signal
import subprocess
import sys
import threading
from collections.abc import Callable
from pathlib import Path

try:
    from watchdog.events import FileSystemEventHandler
    from watchdog.observers import Observer
except ImportError:
    raise SystemExit("The 'watchdog' package is required")

FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
DEFAULT_DEBOUNCE_SECONDS = 55.0

logger = logging.getLogger(__name__)

def is_fastq(path: str | Path) -> bool:
    """Check whether a path has a FASTQ file suffix (gzipped or not).

    Args:
        path: File name or path.

    Returns:
        ``True`` for ``.fastq``, ``.fq``, ``.fastq.gz`` and ``.fq.gz``.
    """
    return str(path).endswith(FASTQ_SUFFIXES)


class DebounceFastqHandler(FileSystemEventHandler):
    """File-system event handler that batches FASTQ arrivals.

    Every new FASTQ file restarts a timer. When the timer runs out,
    ``trigger_cb`` is called once with all FASTQ paths seen since the last
    call.
    """

    def __init__(self, debounce_seconds: float, trigger_cb: Callable[[list[str]], None]):
        """Create the handler.

        Args:
            debounce_seconds: Quiet period to wait after the last FASTQ event.
            trigger_cb: Called with the list of new FASTQ paths.
        """
        super().__init__()
        self.debounce_seconds = float(debounce_seconds)
        self.trigger_cb = trigger_cb
        self._lock = threading.Lock()
        self._pending: set[str] = set()
        self._timer: threading.Timer | None = None

    def on_created(self, event):
        """Queue a newly created file."""
        if not event.is_directory:
            self._queue(str(event.src_path))

    def on_moved(self, event):
        """Queue a file that was renamed into place.

        rsync, scp and most sync tools write to a temporary name
        (e.g. ``.x.fastq.gz.AbC123``) and rename it when the copy is
        complete. That arrives as a move event with the final name in
        ``dest_path``, not as a create event.
        """
        if not event.is_directory:
            self._queue(str(event.dest_path))

    def _queue(self, path: str) -> None:
        """Add a FASTQ path to the pending set and restart the timer."""
        if not is_fastq(path):
            logger.debug("Ignored file event: %s", path)
            return
        logger.info("Detected new FASTQ: %s", path)
        with self._lock:
            self._pending.add(path)
            if self._timer:
                self._timer.cancel()
            self._timer = threading.Timer(self.debounce_seconds, self._on_timeout)
            self._timer.daemon = True
            self._timer.start()

    def _on_timeout(self) -> None:
        """Hand the pending paths to ``trigger_cb`` and clear them."""
        with self._lock:
            pending = list(self._pending)
            self._pending.clear()
            self._timer = None
        logger.info("Debounce timeout, triggering pipeline for %d files", len(pending))
        try:
            self.trigger_cb(pending)
        except Exception:
            logger.exception("Error while triggering pipeline")


class WatchAndRun:
    """Watch a folder and run the pipeline script, one run at a time."""

    def __init__(self, script: Path, input_dir: Path, output_dir: Path,
                 debounce: float = DEFAULT_DEBOUNCE_SECONDS, extra_args: list[str] | None = None):
        """Set up the watcher; call ``start()`` to begin watching.

        Args:
            script: Pipeline script to run (``run_pipeline.py``).
            input_dir: Folder to watch; also passed to the script as ``-i``.
            output_dir: Passed to the script as ``-o``.
            debounce: Seconds to wait after the last new FASTQ before a run.
            extra_args: Extra command-line arguments for the script.
        """
        self.script = script
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.debounce = debounce
        self.extra_args = extra_args or []
        self._observer = Observer()
        self._run_lock = threading.Lock()
        self._run_in_progress = False
        self._pending_during_run = False

    def _invoke_pipeline(self) -> None:
        """Run the pipeline script once and log how it ended."""
        cmd = [sys.executable, str(self.script), "-i", str(self.input_dir),
               "-o", str(self.output_dir), *self.extra_args]
        logger.info("Invoking pipeline: %s", " ".join(cmd))
        try:
            result = subprocess.run(cmd, check=False)
        except OSError:
            logger.exception("Could not start the pipeline")
            return
        if result.returncode == 0:
            logger.info("Pipeline finished successfully")
        else:
            logger.error("Pipeline finished with returncode=%s", result.returncode)

    def request_run(self) -> None:
        """Start a pipeline run in the background, or queue one.

        If a run is already in progress, one more run is started right after
        it finishes. Any number of requests during a run add up to a single
        follow-up run.
        """
        with self._run_lock:
            if self._run_in_progress:
                logger.info("A pipeline run is in progress; will run again when it finishes")
                self._pending_during_run = True
                return
            self._run_in_progress = True
        threading.Thread(target=self._run_until_idle, daemon=True).start()
    def _run_until_idle(self) -> None:
        """Run the pipeline repeatedly until no new run has been requested."""
        while True:
            try:
                self._invoke_pipeline()
            except Exception:
                logger.exception("Unexpected error while running the pipeline")
            with self._run_lock:
                if not self._pending_during_run:
                    self._run_in_progress = False
                    return
                self._pending_during_run = False
            logger.info("New files arrived during the run; starting another run")

    def start(self) -> None:
        """Start watching ``input_dir`` (recursively) in a background thread."""
        handler = DebounceFastqHandler(self.debounce, lambda _files: self.request_run())
        self._observer.schedule(handler, str(self.input_dir), recursive=True)
        logger.info("Starting observer on %s (debounce=%ss)", self.input_dir, self.debounce)
        self._observer.start()

    def stop(self) -> None:
        """Stop watching for new files."""
        logger.info("Stopping observer")
        try:
            self._observer.stop()
            self._observer.join(timeout=5)
        except Exception:
            logger.exception("Error stopping observer")


def main() -> None:
    """Parse arguments, run once for files already present, then watch.

    Runs until stopped with Ctrl+C or SIGTERM.

    Raises:
        SystemExit: If the input folder does not exist, or on shutdown.
    """
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
    logger.info("Starting WatchAndRun for IRMA pipeline")
    parser = argparse.ArgumentParser(description="Watch an input directory and run IRMA pipeline when FASTQ samples arrive")
    parser.add_argument("-i", "--input", required=True, type=Path, help="Input directory to watch for FASTQ files")
    parser.add_argument("-o", "--output", required=True, type=Path, help="Output directory for IRMA pipeline results")
    parser.add_argument("--debounce", type=float, default=DEFAULT_DEBOUNCE_SECONDS,
                        help="Seconds to wait after the last new FASTQ file before running the pipeline")
    parser.add_argument("--pipeline-script", type=Path, default=Path(__file__).parent / "run_pipeline.py",
                        help="Path to the IRMA pipeline script to invoke")
    parser.add_argument("--dry-run", action="store_true",
                        help="Pass --dry-run to the pipeline script")
    parser.add_argument("--config", required=True, type=Path,
                        help="Path to JSON config file for default parameters (e.g. conda envs)")

    args = parser.parse_args()

    input_dir = args.input.resolve()
    output_dir = args.output.resolve()
    if not input_dir.exists():
        raise SystemExit(f"Input directory not found: {input_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    extra = ["--config", str(args.config.resolve())]
    if args.dry_run:
        extra.append("--dry-run")

    watcher = WatchAndRun(script=args.pipeline_script, input_dir=input_dir,
                          output_dir=output_dir, debounce=args.debounce, extra_args=extra)

    def _term(signum, frame):
        logger.info("Received signal %s, shutting down", signum)
        watcher.stop()
        sys.exit(0)

    signal.signal(signal.SIGINT, _term)
    signal.signal(signal.SIGTERM, _term)

    existing_fastqs = [p for p in input_dir.rglob("*") if p.is_file() and is_fastq(p)]
    if existing_fastqs:
        logger.info("Existing FASTQ files found on startup (%d), triggering one initial run", len(existing_fastqs))
        watcher.request_run()

    watcher.start()
    while True:
        signal.pause()


if __name__ == "__main__":
    main()