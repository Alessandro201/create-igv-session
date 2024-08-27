#!/usr/bin/env python
from collections.abc import Iterable
import json
from pathlib import Path
import argparse
import re
from typing import Union


DEFAULT_STYLESHEET = Path(__file__).parent / "stylesheet.json"
DEFAULT_SESSION_TEMPLATE = Path(__file__).parent / "session_template.json"

BASE_FILETYPES_ALLOWED = [
    ".bed",
    ".bed.gz",
    ".bigWig",
    ".bigWig.gz",
    ".gtf",
    ".gtf.gz",
]


BASE_ORDER_VALUE = 10
BASE_TRACK_FEATURES = {
    ".bed": {
        "format": "bed",
        "type": "annotation",
        "height": 25,
    },
    ".bed.gz": {
        "format": "bed",
        "type": "annotation",
        "height": 25,
    },
    ".bigWig": {
        "format": "bigWig",
        "type": "wig",
    },
    ".bigWig.gz": {
        "format": "bigWig",
        "type": "wig",
    },
    ".gtf": {},
}


def parse_args():
    parser = argparse.ArgumentParser(
        prog="create-igv-session.py",
        description="Search in the given directories for track to display on IGV and prepare an IGV session",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        help="Where to save the igv session file",
        default="./igv-session.json",
    )
    parser.add_argument(
        "-p",
        "--prefix",
        type=str,
        help="Prefix the localhost URL to the (relative) paths of the tracks found",
        default="http://localhost:8001/",
    )
    parser.add_argument(
        "-s",
        "--stylesheet",
        type=str,
        help="json style sheet to specify the style of the tracks",
        default=DEFAULT_STYLESHEET,
    )
    parser.add_argument(
        "-t",
        "--template",
        type=str,
        help="json base template of the session",
        default=DEFAULT_SESSION_TEMPLATE,
    )
    parser.add_argument(
        "-f",
        "--force",
        action="store_true",
        help="Overwrite the output if it exists.",
    )
    parser.add_argument("sources", help="Directories to search for tracks", nargs="+", type=str)

    args = parser.parse_args()

    args.output = Path(args.output)
    args.stylesheet = Path(args.stylesheet)
    args.template = Path(args.template)
    args.sources = [Path(s) for s in args.sources]
    args.prefix = args.prefix if args.prefix.endswith("/") else args.prefix + "/"

    if args.output.exists() and not args.force:
        print("The output already exists. If you want to overwrite it use the `-f` or `--force` flag")
        exit(1)

    if not args.stylesheet.exists():
        print("The given stylesheet file does not exists")
        exit(1)

    for src in args.sources:
        if not src.is_dir():
            print(f"The given source is not a valid directory: {src}")
            exit(1)

    return args


def read_json(stylesheet: Path) -> dict:
    text = Path(stylesheet).read_text()
    return json.loads(text)


def set_base_track_featues(file_path: str, order: int, prefix: str = "") -> dict:
    track = {}
    track["name"] = Path(file_path).name
    track["url"] = prefix + file_path

    for extension, track_features in BASE_TRACK_FEATURES.items():
        if file_path.endswith(extension):
            track.update(track_features)

    track["order"] = order
    return track


def update_track_features(track: "dict[str, str | int]", track_features: "dict[str, dict]"):
    for pattern, feature in track_features.items():
        track_name = track.get("name")

        # The track name should be set by set_base_track_features()
        assert track_name is not None

        if not isinstance(track_name, str):
            print(f"Track name: {track_name} is not a string: type {type(track_name)}")
            exit(1)

        if re.match(pattern, track_name):
            track.update(feature)
    return track


def add_to_template(template_path: str, tracks: "list[dict[str, str | int]]") -> dict:
    template = json.loads(Path(template_path).read_text())

    if template.get("tracks") is None:
        print(
            f"In the template there must be a section called 'tracks'. Check out the base template"
            f"to see an example: {DEFAULT_SESSION_TEMPLATE}"
        )
        exit(1)

    template["tracks"].extend(tracks)
    return template


def multikey_sort(string: str, patterns: "list[list[str]]") -> "list[tuple[Union[int, str]]]":
    """
    Match a string against the patterns, and return a "key" which can be used for sorting.

    The "patterns" parameter is a list of list of strings. It can be thought of as groups of patterns,
    and each group will try to match its patterns in order. The first to match will add its index to the key,
    and if no one does the a number just higher that the last index will be added.

    This way allows to sort the items in a tree-like fashion:

                  _ 0 -> (0, 0)
            _ 0 _/
           /     \\_ 1 -> (0, 1)
    tree _/
          \\       _ 0 -> (1, 0)
           \\_ 1 _/
                 \\_ 1 -> (1, 1)


    The patterns may also contain regex capture groups, however the patterns of each group need to have
    the same number of capture groups, otherwise they can't be compared and you will get an error.


    Example:
        patterns = [[".*abc.*", ".*def.*"], [".*txt"]]
        str1 = "abc.txt"  # multikey_sort -> [(0,), (0,)]
        str2 = "def.txt"  # multikey_sort -> [(1,), (0,)]
        str3 = "abc.gz"   # multikey_sort -> [(0,), (1,)]

        sorted((str1, str2, str3), key=lambda x: multikey_sort(x, patterns))  # str1, str3, str2
    """

    if not patterns:
        return [(string,)]

    multikey = []
    for pattern_group in patterns:
        matched = False
        for i, pattern in enumerate(pattern_group):
            if match_obj := re.match(pattern, string):
                matched = True
                if groups := match_obj.groups():
                    multikey.append(groups)
                else:
                    multikey.append((i,))
                break

        if not matched:
            multikey.append((len(pattern_group),))

    return multikey


def is_excluded(string: str, exclusion_patterns: "list[str]") -> bool:
    """
    Checks whether the string matches any exclusion pattern
    """
    return any((re.match(pattern, string) for pattern in exclusion_patterns))


def filetype_is_allowed(path: str, extensions: Iterable[str]) -> bool:
    """
    Checks whether the file extension is in the extensions allowed
    """
    return path.endswith(tuple(extensions))


def sort_paths(paths: "Iterable[str]", style: dict, remove_rep_number: bool = True) -> "list[str]":
    """ """

    # The default is to match the part before the first dot
    pattern: str = style.get("group_by_regex", r"^([^.]*).*$")

    groups = {}
    for path in paths:
        name = Path(path).name
        if remove_rep_number:
            name = re.sub(r"_R\d+\.", ".", name)

        match = re.search(pattern, name)
        group_key = match.group(1) if match else ""

        group_paths = groups.get(group_key, list())
        group_paths.append(path)
        groups[group_key] = group_paths

    sort_patterns: Union[list[list[str]], None] = style.get("set_track_order_within_group", None)

    paths = list()
    for group_key in sorted(groups.keys(), key=lambda p: Path(p).name):
        group_paths = groups[group_key]
        paths.extend(
            sorted(
                group_paths,
                key=lambda path: multikey_sort(Path(path).name, sort_patterns),
            )
        )

    return paths


def main() -> int:
    args = parse_args()

    paths = [str(p) for src in args.sources for p in Path(src).rglob("*")]
    style = read_json(args.stylesheet)
    exclude_patterns: list[str] = style.get("exclude", [])
    filetypes_allowed: list[str] = style.get("filetypes_allowed", BASE_FILETYPES_ALLOWED)
    paths = filter(lambda p: filetype_is_allowed(p, filetypes_allowed), paths)
    paths = filter(lambda p: not is_excluded(p, exclude_patterns), paths)
    paths = sort_paths(paths, style)

    tracks = [set_base_track_featues(path, prefix=args.prefix, order=i + 10) for i, path in enumerate(paths)]
    track_features = style.get("track_features", {})
    tracks = [update_track_features(t, track_features) for t in tracks]

    json_session = add_to_template(args.template, tracks)
    Path(args.output).write_text(json.dumps(json_session, indent=2))
    return 0
