"""
Shared helpers for the ALEdb data-QC pipeline.

Provides configuration loading, path resolution, logging, table IO,
checkpoint / removed-row writers and basic metric computation. Imported by every
step script via importlib (the file name starts with a digit, so it cannot be
imported by name).
"""
import os
import re
import sys
import json
import logging
import time
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

# The folder that holds every output of the run, output/qc/<run_name>/ under the folder the step
# is started in. Set by load_config.
_RUN_DIR = None
# A run name becomes a folder name, so it may not hold a path separator or start with a dot.
RUN_NAME = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*")


# Every setting of the pipeline, as dotted paths, with the kind of value it takes. The code
# holds no fallback values: each of these has to be written in the config with a value of
# its kind, and a config with any other key stops the run, so a misspelled, missing or
# wrongly written setting cannot silently change the result. A listed key is one setting,
# so a free-form mapping such as step01.deposit_overrides is one entry. A tuple lists the
# only words the setting accepts. A path is not looked up here: a path that does not exist
# stops the step that reads it.
KNOWN_KEYS = {
    "run_name": "run name",
    "mutil_dir": "text",
    "organism.species": "text", "organism.subspecies": "list",
    "organism.subspecies_table": "text or empty",
    "inputs.muts_all_dir": "text", "inputs.metadata_dir": "text",
    "inputs.species_normalization": "text", "inputs.condition_normalization": "text",
    "inputs.string_aliases": "text",
    "config_tables.starting_strain_mutation_table": "text",
    "config_tables.starting_strain_mutation_table_user": "text or empty",
    "config_tables.hypermutator_genes": "text",
    "config_tables.manual_experiment_exclusion": "text",
    "step01.deposit_overrides": "ExpIDs",
    "step02.remove_condition_groups": "list", "step02.match_cols": "list",
    "step03.strict_key_cols": "list", "step03.starting_strain_sample_label": "text",
    "step03.drop_reference_seqs": "list", "step03.gene_name_replacements": "mapping",
    "step04.high_mutation_count_threshold": "count",
    "step04.population_frequency_threshold": "fraction",
    "step04.hypermutator_frequency_token": "text", "step04.drop_gene_names": "list",
    "step04.high_prevalence_report_threshold": "fraction",
    "step04.high_prevalence_min_endpoints": "count",
    "step05.remove_projects": "list", "step05.string_valid_sources": "list",
    "step05.string_primary_sources": "list", "step05.clean_gene_strings": "switch",
    "step05.keep_unmapped_as_annotated": "switch", "step05.locus_tag_column": "text",
}
# The number kinds and the values they allow.
NUMBERS = {"count": "a whole number of at least 1", "fraction": "a number above 0 and at most 1",
           "significance level": "a number above 0 and below 1"}
# How to write "none" for the kinds that have an empty value.
EMPTY = {"text or empty": '""', "list": "[]", "mapping": "{}", "ExpIDs": "{}"}
# Names that are not settings of this pipeline, each with the text shown when a config has one.
REMOVED_KEYS = {}


class UniqueKeyLoader(yaml.SafeLoader):
    """YAML loader that stops on a key written twice in one mapping. The plain loader keeps
    the last one and says nothing."""

    def construct_mapping(self, node, deep=False):
        first = {}
        for key_node, _ in node.value:
            if key_node.tag == "tag:yaml.org,2002:merge":
                continue
            key = self.construct_object(key_node, deep=True)
            line = key_node.start_mark.line + 1
            try:
                if key in first:
                    raise ValueError(f"setting {key!r} is written twice in "
                                     f"{key_node.start_mark.name}: lines {first[key]} and {line}")
                first[key] = line
            except TypeError:       # an unhashable key, the plain loader reports it
                pass
        return super().construct_mapping(node, deep)


# The last part of every path in KNOWN_KEYS and REMOVED_KEYS ("drop_gene_names", ...).
SETTING_NAMES = {p.rsplit(".", 1)[-1] for p in KNOWN_KEYS.keys() | REMOVED_KEYS.keys()}


def unknown_keys(d, prefix=""):
    """(dotted path, problem) for every key of d that is not a setting where it stands.
    problem is "dot" for a key written with a dot ("step04.drop_gene_names: []" on one line,
    which no step reads), "inside" for the name of a setting found inside a free-form mapping
    (it was indented one level too deep and would be ignored), and "unknown" otherwise."""
    out = []
    for k, v in d.items():
        path = f"{prefix}{k}"
        if "." in str(k):
            out.append((path, "dot"))
        elif any(p.startswith(path + ".") for p in KNOWN_KEYS):     # a section
            if isinstance(v, dict):
                out.extend(unknown_keys(v, path + "."))
        elif path not in KNOWN_KEYS:
            out.append((path, "unknown"))
        elif isinstance(v, dict):
            out.extend((f"{path}.{k2}", "inside") for k2 in v if k2 in SETTING_NAMES)
    return out


def lookup(cfg, path):
    """(written, value) of the setting at a dotted path."""
    node = cfg
    for part in path.split("."):
        if not (isinstance(node, dict) and part in node):
            return False, None
        node = node[part]
    return True, node


def missing_keys(cfg):
    """(dotted path, problem) for every required setting that is absent ("missing") or written
    with nothing after the colon ("blank"). An empty list, mapping or string ([], {}, "")
    counts as written."""
    out = []
    for path in sorted(KNOWN_KEYS):
        written, value = lookup(cfg, path)
        if not written:
            out.append((path, "missing"))
        elif value is None:
            out.append((path, "blank"))
    return out


def _shown(v):
    if isinstance(v, bool):
        return f"{str(v).lower()} (YAML reads yes, no, on and off as true or false)"
    return f"the text {v!r}" if isinstance(v, str) else repr(v)


def _why_not_text(x):
    """What to say about x, found where text is expected."""
    if x is None:
        return "found an empty entry. Write it or remove it"
    if isinstance(x, bool):
        return (f"found {str(x).lower()}, which YAML makes of yes, no, on, off, true and false. "
                "Put the text in quotation marks")
    if isinstance(x, (int, float)):
        return (f"found the number {x!r}, which YAML makes of unquoted digits (1.00 becomes 1.0). "
                "Put the text in quotation marks, exactly as intended")
    return f"found {x!r}"


def kind_problem(v, kind):
    """What is wrong with the value v written for a setting of this kind, or None."""
    whole = isinstance(v, int) and not isinstance(v, bool)
    if isinstance(kind, tuple):
        return None if v in kind else f"must be one of {', '.join(kind)}, found {_shown(v)}"
    if kind == "run name":
        if not isinstance(v, str):
            return f"must be text, {_why_not_text(v)}"
        return None if RUN_NAME.fullmatch(v) else (
            "must be a folder name made of letters, digits, '_', '-' and '.', starting with a "
            f"letter or digit, found {v!r}")
    if kind == "switch":
        return None if isinstance(v, bool) else f"must be true or false, found {_shown(v)}"
    if kind in NUMBERS:
        if isinstance(v, float) and not np.isfinite(v):
            return f"must be a finite number, found {v!r}"
        if whole or (kind != "count" and isinstance(v, float)):
            in_range = {"count": v >= 1, "fraction": 0 < v <= 1, "significance level": 0 < v < 1}
            return None if in_range[kind] else f"must be {NUMBERS[kind]}, found {v!r}"
        if isinstance(v, str):
            hint = ". Write a plain number without quotation marks" + (
                ", such as 0.05 rather than 5e-2" if kind != "count" else "")
        elif isinstance(v, float) and v.is_integer() and v >= 1:
            hint = f". Write {int(v)} without the decimal point"
        else:
            hint = ""
        return f"must be {NUMBERS[kind]}, found {_shown(v)}{hint}"
    if kind == "list":
        if not isinstance(v, list):
            return f"must be a list such as [a, b], or [] for none, found {_shown(v)}"
        entries = v
    elif kind in ("mapping", "ExpIDs"):
        if not isinstance(v, dict):
            return (f"must be a mapping, one 'name: value' per line, or {{}} for none, "
                    f"found {_shown(v)}")
        if kind == "ExpIDs":
            bad = [x for x in v.values() if not (isinstance(x, int) and not isinstance(x, bool))]
            if bad:
                return f"every ExpID must be a whole number, found {_shown(bad[0])}"
        entries = [*v, *v.values()] if kind == "mapping" else list(v)
    elif not isinstance(v, str):
        return f"must be text, {_why_not_text(v)}"
    else:
        entries = []
    bad = [x for x in entries if not isinstance(x, str)]
    if bad:
        return f"every entry must be text, {_why_not_text(bad[0])}"
    if kind == "text" and not v:
        return "must not be empty"
    return None


def check_keys(cfg, config_path):
    if not isinstance(cfg, dict):
        raise ValueError(f"{config_path} is empty or is not a mapping of settings")
    lines = []
    for k, problem in unknown_keys(cfg):
        if k in REMOVED_KEYS:
            why = f"removed, {REMOVED_KEYS[k]}"
        elif problem == "dot":
            why = "write it as nested keys under its section, not as one key with a dot"
        elif problem == "inside":
            why = f"this is the name of a setting, but it is indented under {k.rsplit('.', 1)[0]}"
        else:
            why = ("not a setting of this pipeline, check the spelling and the section it is "
                   "nested under")
        lines.append(f"  {k}: {why}")
    for k, problem in missing_keys(cfg):
        if problem == "blank":
            kind = KNOWN_KEYS[k]
            why = "no value, write one" + (f", or {EMPTY[kind]} for none" if kind in EMPTY else "")
        else:
            why = "missing, every setting has to be written in the config"
        lines.append(f"  {k}: {why}")
    for k, kind in KNOWN_KEYS.items():
        _, value = lookup(cfg, k)
        problem = kind_problem(value, kind) if value is not None else None
        if problem:
            lines.append(f"  {k}: {problem}")
    if lines:
        raise ValueError(f"settings to fix in {config_path}:\n" + "\n".join(lines))


def load_config(config_path):
    """Load the YAML config given with --config. Every output of the run goes to
    output/qc/<run_name>/ under the folder the step is started in."""
    global _RUN_DIR
    if not config_path:
        raise ValueError("no config file given, pass one with --config")
    config_path = Path(config_path)
    if not config_path.is_absolute():
        config_path = (Path.cwd() / config_path).resolve()
    with open(config_path) as f:
        cfg = yaml.load(f, Loader=UniqueKeyLoader)
    check_keys(cfg, config_path)
    cfg["_config_path"] = str(config_path)
    _RUN_DIR = Path.cwd() / "output" / "qc" / cfg["run_name"]
    return cfg


def resolve(rel_or_abs):
    """Resolve a config path. A relative path is read from the folder the step is started in,
    like a path typed on the command line."""
    p = Path(rel_or_abs)
    if p.is_absolute():
        return p
    return (Path.cwd() / p).resolve()


# Every output of a run sits in its own folder, output/qc/<run_name>/.
def out_dir():
    return _RUN_DIR


def removed_dir():
    return out_dir() / "removed"


def qc_dir():
    return out_dir() / "qc"


def checkpoints_dir():
    return out_dir() / "checkpoints"


def logs_dir():
    return out_dir() / "logs"


def organism_input_warning(cfg, logger):
    """Warn loudly if a non-E.coli organism is selected while the config still
    points at E. coli-specific reference files / tables."""
    species = cfg["organism"]["species"]
    if species == "Escherichia coli":
        return
    flags = []
    for key in ("string_aliases",):
        v = str(cfg["inputs"].get(key, ""))
        if "511145" in v:
            flags.append(f"inputs.{key} = {v}")
    if "CP009273" in str(cfg["step03"]["drop_reference_seqs"]):
        flags.append("step03.drop_reference_seqs still contains CP009273 (E. coli)")
    if cfg["step05"]["locus_tag_column"] == "b_number":
        flags.append("step05.locus_tag_column is 'b_number' (E. coli), set a generic name e.g. 'locus_tag'")
    curated = resolve(cfg["config_tables"]["starting_strain_mutation_table"])
    if sum(1 for _ in open(curated)) > 1:
        flags.append(f"curated starting strain mutation table is non-empty: {curated}")
    if flags:
        logger.warning("ORGANISM CHECK: species='%s' but the config still points at E. coli-specific inputs:", species)
        for f in flags:
            logger.warning("  - %s", f)
        logger.warning("  Update these to your organism's reference files/tables before trusting results.")


def setup_logging(step_id, step_name):
    """Create a logger that writes to logs/<step_id>_<step_name>.log and stdout, with times in UTC.
    The log is the first file a step writes, so a step that already ran in this run folder stops
    here. An error that stops the step is written to the log, with where it happened."""
    logs_dir().mkdir(parents=True, exist_ok=True)
    log_file = logs_dir() / f"{step_id}_{step_name}.log"
    _stop_if_exists(log_file)
    logger = logging.getLogger(f"qc.{step_id}")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    fmt = logging.Formatter("%(asctime)s  %(levelname)s  %(message)s", "%Y-%m-%dT%H:%M:%SZ")
    fmt.converter = time.gmtime
    fh = logging.FileHandler(log_file, mode="x")
    fh.setFormatter(fmt)
    logger.addHandler(fh)
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    logger.addHandler(sh)
    log_uncaught_errors(logger, fh)
    logger.info("step %s (%s) started", step_id, step_name)
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


def add_mutil_to_path(cfg):
    """Put mutil_dir, the folder that holds the mutil package, on sys.path."""
    mutil_dir = resolve(cfg["mutil_dir"])
    if str(mutil_dir) not in sys.path:
        sys.path.insert(0, str(mutil_dir))
    return mutil_dir


# ---------------------------------------------------------------------
# Table IO
# ---------------------------------------------------------------------
def _sep_for(path):
    return "," if str(path).endswith(".csv") else "\t"


def read_table(path, **kw):
    path = Path(path)
    return pd.read_csv(path, sep=_sep_for(path), low_memory=False, **kw)


def _stop_if_exists(path):
    if Path(path).exists():
        raise FileExistsError(f"{path} already exists. A run never overwrites a file: delete the "
                              "run folder or use a new run_name")


def _write_new(path, write):
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
    if not len(df.columns):
        raise ValueError(f"{path}: a table with no columns cannot be read back, so it is not written")
    return _write_new(path, lambda tmp: df.to_csv(tmp, sep=_sep_for(path), index=False))


def save_removed(df_removed, path, reason):
    """Write a removed-rows table, annotating it with a removal reason column."""
    df_removed = df_removed.copy()
    if "removal_reason" not in df_removed.columns:
        df_removed["removal_reason"] = reason
    return write_table(df_removed, path)


def write_json(obj, path):
    def write(tmp):
        with open(tmp, "w") as f:
            json.dump(obj, f, indent=2, default=_json_default)
    return _write_new(path, write)


def write_checkpoint(metrics, path):
    return write_json(metrics, path)


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, (set,)):
        return sorted(o)
    return str(o)


# ---------------------------------------------------------------------
# Metric helpers
# ---------------------------------------------------------------------
def sample_cols(df):
    """The columns that identify a sample: sample_endpoint, which holds exp, or else exp and
    sample, since step 03 rewrites sample without exp."""
    return ["sample_endpoint"] if "sample_endpoint" in df.columns else (
        ["exp", "sample"] if "sample" in df.columns else None)


def gene_col(df):
    return "Gene_name_normalized" if "Gene_name_normalized" in df.columns else (
        "Gene (Scrollable)" if "Gene (Scrollable)" in df.columns else None)


def n_unique(df, col):
    return int(df[col].nunique()) if col and col in df.columns else None


def basic_metrics(df):
    sc = sample_cols(df)
    gc = gene_col(df)
    return {
        "n_rows": int(len(df)),
        "n_columns": int(df.shape[1]),
        "n_unique_experiments": n_unique(df, "exp"),
        "n_unique_samples": int(df.groupby(sc).ngroups) if sc else None,
        "sample_id_columns": sc,
        "n_unique_genes": n_unique(df, gc),
        "gene_column": gc,
    }


# ---------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------
def validate_common(df, logger, required_cols=None):
    problems = []
    if df.columns.duplicated().any():
        dups = df.columns[df.columns.duplicated()].tolist()
        problems.append(f"duplicated column names: {dups}")
    if required_cols:
        missing = [c for c in required_cols if c not in df.columns]
        if missing:
            problems.append(f"missing required columns: {missing}")
    if not len(df):
        problems.append("the table has no rows")
    for c in ["Position", "exp"]:
        if c in df.columns and len(df) and df[c].isna().all():
            problems.append(f"key column '{c}' is entirely empty")
    for p in problems:
        logger.error("VALIDATION: %s", p)
    if problems:
        raise ValueError("validation failed: " + " | ".join(problems))


def check_conservation(n_before, n_after, n_removed, logger):
    """removed rows + retained rows == input rows. Stops the run on a mismatch,
    since every removal is recorded and a mismatch therefore means rows were lost or duplicated."""
    if n_after + n_removed != n_before:
        msg = (f"row conservation mismatch: before={n_before} after={n_after} removed={n_removed} "
               f"(after+removed={n_after + n_removed})")
        logger.error("VALIDATION: %s", msg)
        raise ValueError(msg)
    return True


def finish_step(step_id, step_name, logger):
    logger.info("step %s (%s) finished", step_id, step_name)
