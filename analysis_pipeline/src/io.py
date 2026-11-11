"""
IO helpers for the pairwise analysis pipeline.

Loads the curated QC table, checks required columns, builds globally unique
sample IDs (and verifies uniqueness across experiments), and provides file
writers that never overwrite, a small logger and terminal output helpers.
"""
import json
import logging
import os
import sys
import time
from pathlib import Path

import pandas as pd


def setup_logging(name, log_dir):
    """Full detail (INFO) goes to log_dir/<name>.log, which must not exist yet. The
    terminal shows only warnings and errors. Short summaries are printed to the
    terminal separately via console(), and the progress bar via progress(), so the
    run stays easy to read."""
    Path(log_dir).mkdir(parents=True, exist_ok=True)
    log_file = Path(log_dir) / f"{name}.log"
    _stop_if_exists(log_file)
    logger = logging.getLogger(f"analysis.{name}")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    logger.propagate = False
    fh = logging.FileHandler(log_file, mode="x")
    fh.setLevel(logging.INFO)
    fmt = logging.Formatter("%(asctime)s  %(levelname)s  %(message)s", "%Y-%m-%dT%H:%M:%SZ")
    fmt.converter = time.gmtime
    fh.setFormatter(fmt)
    logger.addHandler(fh)
    sh = logging.StreamHandler(sys.stdout)
    sh.setLevel(logging.WARNING)
    sh.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
    logger.addHandler(sh)
    log_uncaught_errors(logger, fh)
    return logger


def log_uncaught_errors(logger, file_handler):
    """Write an error that stops the program, with its traceback, to the log file, then let Python
    show it on the terminal as usual."""
    def hook(exc_type, exc, tb):
        record = logger.makeRecord(logger.name, logging.ERROR, "", 0, "run stopped: %s: %s",
                                   (exc_type.__name__, exc), (exc_type, exc, tb))
        file_handler.handle(record)
        sys.__excepthook__(exc_type, exc, tb)
    sys.excepthook = hook


def console(msg=""):
    """Print one clean line to the terminal (headlines / summaries), kept
    separate from the detailed file log. Uses tqdm.write when available so it
    never corrupts an active progress bar."""
    try:
        from tqdm import tqdm
    except ImportError:
        print(msg)
        return
    tqdm.write(str(msg))


def progress(iterable, **kwargs):
    """Wrap an iterable in a single-line tqdm progress bar if tqdm is installed. Without tqdm
    the iterable is returned unchanged."""
    try:
        from tqdm import tqdm
    except ImportError:
        return iterable
    return tqdm(iterable, **kwargs)


def load_table(cfg, logger=None):
    """Load the QC output, check required columns, add sample_unique."""
    df = pd.read_csv(cfg.input_table, sep="\t" if str(cfg.input_table).endswith(".tsv") else ",",
                     low_memory=False)
    required = [cfg.sample_col, cfg.unit_col, cfg.exp_col]
    if cfg.members_col:
        required.append(cfg.members_col)
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise KeyError(f"input table is missing required columns: {missing}")
    # A row with no experiment, sample or unit would be left out of every matrix.
    blank = [c for c in (cfg.exp_col, cfg.sample_col, cfg.unit_col) if df[c].isna().any()]
    if blank:
        raise ValueError(f"input table has rows with no value in {blank}")

    # globally unique sample id: exp + "__" + sample_col
    df["sample_unique"] = df[cfg.exp_col].astype(str) + "__" + df[cfg.sample_col].astype(str)

    if logger:
        logger.info("loaded %s: %d rows, %d columns", cfg.input_table, len(df), df.shape[1])
        logger.info("unique experiments=%d unique samples(%s)=%d unique units(%s)=%d",
                    df[cfg.exp_col].nunique(), cfg.sample_col, df[cfg.sample_col].nunique(),
                    cfg.unit_col, df[cfg.unit_col].nunique())
    check_sample_uniqueness(df, cfg, logger)
    return df


def check_sample_uniqueness(df, cfg, logger=None):
    """A sample_unique value must belong to exactly one experiment. A sample_col value shared
    by two experiments is only reported, since each experiment is analyzed on its own."""
    per = df.groupby("sample_unique")[cfg.exp_col].nunique()
    bad = int((per > 1).sum())
    per_raw = df.groupby(cfg.sample_col)[cfg.exp_col].nunique()
    bad_raw = int((per_raw > 1).sum())
    msg = (f"sample_unique values mapping to >1 experiment: {bad}, "
           f"{cfg.sample_col} values colliding across experiments: {bad_raw}")
    if bad:
        raise ValueError(msg)
    if logger:
        logger.info(msg)


def _stop_if_exists(path):
    if Path(path).exists():
        raise FileExistsError(f"{path} already exists. A run never overwrites a file: delete the "
                              "run folder or use a new run_name")


def write_new(path, write):
    """Write a file that must not exist yet. It is written under a temporary name in the same
    folder and renamed at the end, so a write that is stopped never leaves a cut-off file."""
    path = Path(path)
    _stop_if_exists(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        write(tmp)
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)
    return path


def write_table(df, path):
    """A table without row labels. It always keeps its column names, even with no rows."""
    if not len(df.columns):
        raise ValueError(f"{path}: a table with no columns cannot be read back, so it is not written")
    sep = "\t" if str(path).endswith(".tsv") else ","
    return write_new(path, lambda tmp: df.to_csv(tmp, sep=sep, index=False))


def write_matrix(df, path):
    """A matrix with its row labels, such as the unit x sample matrix or a correlation matrix.
    An empty matrix is written as an empty header, which reads back as an empty table."""
    return write_new(path, lambda tmp: df.to_csv(tmp))


def write_json(obj, path):
    def write(tmp):
        with open(tmp, "w") as f:
            json.dump(obj, f, indent=2)
    return write_new(path, write)
