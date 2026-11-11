"""
Configuration loading for the pairwise analysis pipeline.

Reads relative paths from the folder the analysis is started in and exposes a small
Config object. Every output of a run goes to output/analysis/<run_name>/ under that
folder. The test unit is defined only by unit_col, and unit_name is its tag in the
outputs. There is no gene/operon switch anywhere in the pipeline.
"""
import io
import math
import re
from pathlib import Path

import yaml

# A run name becomes a folder name, so it may not hold a path separator or start with a dot.
RUN_NAME = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*")


def resolve(rel_or_abs):
    """A relative path is read from the folder the analysis is started in."""
    p = Path(rel_or_abs)
    return p if p.is_absolute() else (Path.cwd() / p).resolve()


# Every setting of the pipeline, with the kind of value it takes. The code holds no fallback
# values: each of these has to be written in the config with a value of its kind, and a config
# with any other key stops the run, so a misspelled, missing or wrongly written setting cannot
# silently change the result. A tuple lists the only words the setting accepts. A path is not
# looked up here: an input_table that does not exist stops the run when it is read.
KNOWN_KEYS = {
    "run_name": "run name", "input_table": "text",
    "sample_col": "text", "unit_col": "text", "unit_name": "text", "exp_col": "text",
    "members_col": "text or empty",
    "min_endpoints": "count", "fixed_min_count": "count",
    "significance_threshold": "significance level",
    "discover_fdr_method": ("BH", "DBH"),
}
# The number kinds and the values they allow.
NUMBERS = {"count": "a whole number of at least 1", "fraction": "a number above 0 and at most 1",
           "significance level": "a number above 0 and below 1"}
# How to write "none" for the kinds that have an empty value.
EMPTY = {"text or empty": '""'}
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


def check_keys(d, path):
    if not isinstance(d, dict):
        raise ValueError(f"{path} is empty or is not a mapping of settings")
    lines = [f"  {k}: removed, {REMOVED_KEYS[k]}" if k in REMOVED_KEYS
             else f"  {k}: not a setting of this pipeline, check the spelling"
             for k in d if k not in KNOWN_KEYS]
    # A setting indented under another one is read as part of that block and ignored.
    for name, value in d.items():
        if name in KNOWN_KEYS and isinstance(value, dict):
            lines += [f"  {name}.{k}: removed, {REMOVED_KEYS[k]}" if k in REMOVED_KEYS
                      else f"  {name}.{k}: this is the name of a setting, move it to the left margin"
                      for k in value if k in KNOWN_KEYS or k in REMOVED_KEYS]
    # An empty list, mapping or string ([], {}, "") counts as written, a blank value does not.
    for k in sorted(KNOWN_KEYS):
        if k not in d:
            lines.append(f"  {k}: missing, every setting has to be written in the config")
        elif d[k] is None:
            none = f", or {EMPTY[KNOWN_KEYS[k]]} for none" if KNOWN_KEYS[k] in EMPTY else ""
            lines.append(f"  {k}: no value, write one{none}")
    for k, kind in KNOWN_KEYS.items():
        problem = kind_problem(d[k], kind) if d.get(k) is not None else None
        if problem:
            lines.append(f"  {k}: {problem}")
    if lines:
        raise ValueError(f"settings to fix in {path}:\n" + "\n".join(lines))


def _shown(v):
    if isinstance(v, bool):
        return f"{str(v).lower()} (YAML reads yes, no, on and off as true or false)"
    return f"the text {v!r}" if isinstance(v, str) else repr(v)


def _why_not_text(x):
    """What to say about x, found where text is expected."""
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
    if kind in NUMBERS:
        if isinstance(v, float) and not math.isfinite(v):
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
    if not isinstance(v, str):
        return f"must be text, {_why_not_text(v)}"
    if kind == "text" and not v:
        return "must not be empty"
    return None


class Config:
    def __init__(self, d, path, raw):
        check_keys(d, path)
        self.path = path
        self.raw = raw              # the config file's bytes, as read once and parsed
        self.run_name = d["run_name"]
        self.input_table = resolve(d["input_table"])
        # Every output of the run sits in its own folder, the matrices in its matrices/ folder.
        self.output_dir = Path.cwd() / "output" / "analysis" / self.run_name
        self.matrices_dir = self.output_dir / "matrices"
        self.sample_col = d["sample_col"]
        self.unit_col = d["unit_col"]
        self.unit_name = d["unit_name"]
        self.exp_col = d["exp_col"]
        self.members_col = d["members_col"] or None
        self.min_endpoints = int(d["min_endpoints"])
        self.fixed_min_count = int(d["fixed_min_count"])
        self.significance_threshold = float(d["significance_threshold"])
        # DISCOVER FDR method: 'BH' (standard Benjamini-Hochberg) or 'DBH' (DISCOVER's discrete
        # Benjamini-Hochberg, less conservative for discrete tests). Fisher always uses BH.
        self.discover_fdr_method = d["discover_fdr_method"]

    @property
    def unit_tag(self):
        return str(self.unit_name)


def load_config(path):
    """Load the YAML config given with --config. The file is read once, and the run saves those
    same bytes as its copy of the config."""
    if not path:
        raise ValueError("no config file given, pass one with --config")
    path = Path(path)
    if not path.is_absolute():
        path = (Path.cwd() / path).resolve()
    raw = path.read_bytes()
    stream = io.BytesIO(raw)
    stream.name = str(path)         # so a message about the YAML names the file
    d = yaml.load(stream, Loader=UniqueKeyLoader)
    return Config(d, str(path), raw)
