#!/usr/bin/env python
import argparse
import glob
import json
import re
import textwrap
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path
from pprint import pprint
from typing import Union

STYLES_FOLDER = Path(__file__).parent / "styles/"
AVAILABLE_STYLESHEETS = [file for file in STYLES_FOLDER.rglob("*.json") if not file.name.endswith(".template.json")]
AVAILABLE_TEMPLATES = [style for style in STYLES_FOLDER.rglob("*.template.json")]
AVAILABLE_PRESETS = [style.stem for style in AVAILABLE_STYLESHEETS]
HELP_STYLES = Path(STYLES_FOLDER, "help_styles.txt").read_text()


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
    ".narrowPeak": {
        "format": "bed",
        "type": "annotation",
        "height": 25,
    },
    ".broadPeak": {
        "format": "bed",
        "type": "annotation",
        "height": 25,
    },
    ".narrowPeak.gz": {
        "format": "narrowPark",
        "type": "annotation",
        "height": 25,
    },
    ".broadPeak.gz": {
        "format": "broadPeak",
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
    ".gtf": {
        "format": "gtf",
        "type": "annotation",
    },
    ".gff": {
        "format": "gff",
        "type": "annotation",
    },
}
BASE_OVERLAY_TRACK_FEATURES = {"type": "merged", "autoscale": False, "alpha": 0.5, "height": 50}


def parse_args():
    parser = argparse.ArgumentParser(
        prog="create-igv-session.py",
        description="Search in the given directories for track to display on IGV and prepare an IGV session",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--help-styles", help="Print styles descriptions and how to use them", action="store_true")
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        help="Where to save the igv session file",
        default="./igv-session.json",
    )
    parser.add_argument(
        "-u",
        "--url",
        type=str,
        help="localhost URL to prefix to the (relative) paths of the tracks found",
        default="http://localhost:8001/",
    )
    parser.add_argument(
        "-p",
        "--preset",
        type=str,
        help="Preset of stylesheet and associated template. If you don't want any preset use the `--no-preset` flag",
        choices=AVAILABLE_PRESETS,
        default="default",
    )
    parser.add_argument(
        "--no-preset",
        action="store_true",
        help="Remove any preset. If you use this flag you must then give your own stylesheet and template",
    )

    parser.add_argument(
        "-s",
        "--stylesheet",
        type=str,
        help="Define a json stylesheet of the tracks that will be applied AFTER the preset. "
        "If you don't want the preset use the flag `--no-preset`",
    )
    parser.add_argument(
        "-t",
        "--template",
        type=str,
        help="Define a json session template that will replace the one used in the preset. "
        "It's required with `--no-preset`",
    )
    parser.add_argument(
        "-f",
        "--force",
        action="store_true",
        help="Overwrite the output if it exists.",
    )
    parser.add_argument(
        "-d",
        "--debug",
        action="store_true",
        help="Print debuggin information.",
    )
    parser.add_argument("sources", help="Directories to search for tracks", nargs="*", type=str)

    args = parser.parse_args()

    if args.help_styles:
        print(HELP_STYLES)
        exit(0)

    args.output = Path(args.output)

    if args.no_preset and not args.stylesheet:
        print("With `--no-preset` you must give a stylesheet with `--stylesheet`")
        exit(1)
    elif args.no_preset and not args.template:
        print("With `--no-preset` you must give a template with `--template`")
        exit(1)

    preset_stylesheet = None
    preset_template = None
    if args.preset and args.no_preset is False:
        preset_stylesheet: Path = STYLES_FOLDER / (args.preset + ".json")
        preset_template: Path = STYLES_FOLDER / (args.preset + ".template.json")

    if preset_stylesheet and args.stylesheet:
        args.stylesheet = [preset_stylesheet, Path(args.stylesheet)]
    else:
        args.stylesheet = [preset_stylesheet or Path(args.stylesheet)]
    args.template = Path(args.template or preset_template)

    args.url = args.url if args.url.endswith("/") else args.url + "/"

    if args.output.exists() and not args.force:
        print("The output already exists. If you want to overwrite it use the `-f` or `--force` flag")
        exit(1)

    for stylesheet in args.stylesheet:
        if not stylesheet.exists():
            print(f"The given stylesheet file does not exists: {stylesheet}")
            exit(1)

    return args


@dataclass
class Track:
    path: Path
    index: Union[Path, None] = None


@dataclass
class RNABigWig(Track):
    reverse_strand: Union[Path, None] = None


def read_json(stylesheet: Path) -> dict:
    text = Path(stylesheet).read_text()
    return json.loads(text)


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
            if match_obj := re.match(pattern, string, flags=re.IGNORECASE):
                matched = True
                if groups := match_obj.groups():
                    multikey.append(groups)
                else:
                    multikey.append((i,))
                break

        if not matched:
            multikey.append((len(pattern_group),))

    return multikey


def match_any_pattern(trackpath: Track, patterns: "list[str]") -> bool:
    """
    Checks whether the string matches any pattern
    """
    return any((re.match(pattern, trackpath.path.as_posix()) for pattern in patterns))


def sort_paths(trackpaths: Iterable[Track], sort_patterns: "list[list[str]]") -> "list[Track]":
    trackpaths = list(trackpaths)
    trackpaths_sortkeys = [multikey_sort(trackpath.path.as_posix(), sort_patterns) for trackpath in trackpaths]

    # Check to see if any pattern is wrong and produced not sortable results
    for i_pattern, pattern in enumerate(sort_patterns):
        for i_trackpath in range(len(trackpaths_sortkeys)):
            try:
                _ = trackpaths_sortkeys[0][i_pattern] > trackpaths_sortkeys[i_trackpath][i_pattern]
            except TypeError:
                print(
                    textwrap.dedent(f"""
                    ERROR: When the regex pattern 
                      [
                        {
                        "\n\
                        ".join(pattern)
                    }
                      ]
                    is applied to the paths it produced results that cauldn't be compared and sorted.

                    For instance on '{trackpaths[0].path.as_posix()}' it produced {trackpaths_sortkeys[0][i_pattern]}
                    while on '{trackpaths[i_trackpath].path.as_posix()}' it produced a '{
                        trackpaths_sortkeys[i_trackpath][i_pattern]
                    }'

                    This is likely due to the fact that you used a regex capture group capturing it's contents when you didn't want them.
                    For example the pattern '.*\\.bed(\\.gz)?' when applied on 'text.bed' will return (None,) because the regex pattern matched the string, but the capture group did not capture anything.
                    When you want to match a pattern but not capture anything you need to prefix the capture group with '?:' as in '.*\\.bed(?:\\.gz)?'"
                    """)
                )
                exit(1)

    # Sorting
    trackpaths = sorted(
        trackpaths,
        key=lambda trackpath: multikey_sort(trackpath.path.as_posix(), sort_patterns),
    )

    return trackpaths


def group_tracks_with_indexes(paths: "list[Path]") -> "list[Track]":
    paths = sorted(paths, key=lambda p: p.as_posix())

    tracks = [Track(paths[0])]

    for i, path in enumerate(paths[1:]):
        if path.name.endswith(".tbi") and paths[i].as_posix() == path.as_posix()[:-4]:
            tracks[-1].index = path
        else:
            tracks.append(Track(path))

    return tracks


def group_rna_strands(tracks: "list[Track]", rna_track_patterns: "list[str]") -> "list[RNABigWig]":
    rna_tracks: "list[Track | RNABigWig]" = list()
    other_tracks = list()

    for track in tracks:
        matched = False
        for pattern in rna_track_patterns:
            if re.match(pattern, track.path.as_posix()):
                matched = True
                break

        if matched:
            rna_tracks.append(track)
        else:
            other_tracks.append(track)

    paths_dict: dict[str, list[Track]] = {}

    # Finding pairs of forward-reverse tracks
    # All files which are not a pair are coppied to other_tracks
    rna_tracks = sorted(rna_tracks, key=lambda x: x.path.name)
    for rna_track in rna_tracks:
        if not rna_track.path.name.endswith("bigWig"):
            other_tracks.append(rna_track)
            continue

        if "reverse" not in rna_track.path.name and "forward" not in rna_track.path.name:
            other_tracks.append(rna_track)
            continue

        # Using full path to be sure that I select only forward-reverse strand in the same folder
        name = rna_track.path.as_posix().replace("forward", "").replace("reverse", "")
        path_list = paths_dict.get(name, [])
        path_list.append(rna_track)
        paths_dict[name] = path_list

    # Group together the strands in a single RNABigWig object
    rna_tracks: list[RNABigWig] = list()
    for name, path_list in paths_dict.items():
        forward_in_current_group = False

        for rna_track in path_list:
            if "forward" in rna_track.path.name:
                rna_tracks.append(RNABigWig(rna_track.path))
                forward_in_current_group = True

            elif "reverse" in rna_track.path.name:
                if not forward_in_current_group:
                    print(
                        "Found a rnaseq bigwig file without it's forward strand. "
                        "The files should be named exactly the same save for 'forward' and 'reverse'. "
                        "Not considering this file as paired."
                    )
                    break

                rna_tracks[-1].reverse_strand = rna_track.path

    return rna_tracks + other_tracks


def set_base_track_features(trackpath: Track, order: int, prefix: str = "") -> dict:
    track = {}

    if isinstance(trackpath, RNABigWig):
        assert trackpath.reverse_strand is not None

        track.update(BASE_OVERLAY_TRACK_FEATURES)
        track["name"] = f"Overlay {trackpath.path.name.replace('forward', '')}"
        track["type"] = "merged"
        track["tracks"] = [
            set_base_track_features(Track(trackpath.path), order, prefix),
            set_base_track_features(Track(trackpath.reverse_strand), order, prefix),
        ]

        return track

    track["name"] = trackpath.path.stem
    track["filename"] = trackpath.path.name
    track["track_path"] = trackpath.path.as_posix()
    track["order"] = order
    track["url"] = prefix + trackpath.path.as_posix()

    if trackpath.index:
        track["indexUrl"] = prefix + trackpath.index.as_posix()

    for extension, track_features in BASE_TRACK_FEATURES.items():
        if trackpath.path.name.endswith(extension):
            track.update(track_features)

    return track


def update_track_features(track: "dict[str, str | int]", track_features: "dict[str, dict]"):
    is_overlay = track.get("type") == "merged"
    if is_overlay:
        # The track name should be set by set_base_track_features()
        name = track.get("name")
        assert name is not None
        assert isinstance(name, str)

        match_string = name

        # set_base_track_features() should have set the parameter "tracks" which containes the track to be overlayed
        tracks = track.get("tracks")
        try:
            assert tracks is not None
            assert isinstance(tracks, list)  # list[dict[str, int | str]]
        except AssertionError:
            from pprint import pformat

            print(f"Error found, printing data usefule for debug: \ntrack: \n{pformat(track)}")
            raise
        track["tracks"] = [update_track_features(t, track_features) for t in tracks]

    else:
        # The track path should be set by set_base_track_features()
        track_path = track.get("track_path")
        try:
            assert track_path is not None
            assert isinstance(track_path, str)
        except AssertionError:
            from pprint import pformat

            print(f"Error found, printing data usefule for debug: \ntrack: \n{pformat(track)}")
            raise

        match_string = track_path

    for pattern, feature in track_features.items():
        if re.match(pattern, match_string):
            track.update(feature)

    return track


def add_to_template(template_path: str, tracks: "list[dict[str, str | int]]") -> dict:
    template = json.loads(Path(template_path).read_text())

    if template.get("tracks") is None:
        print(
            f"In the template there must be a section called 'tracks'. Check out the base template"
            f"to see an example: {STYLES_FOLDER / 'default.template.json'}"
        )
        exit(1)

    template["tracks"].extend(tracks)
    return template


def main() -> int:
    args = parse_args()

    styles = list()
    for stylesheet in args.stylesheet:
        style = read_json(stylesheet)
        styles.append(style)

    exclude_patterns: list[str] = []
    include_patterns: list[str] = []
    rna_track_patterns: list[str] = []
    overlay_rna_strands: bool = False
    filetypes_allowed: list[str] = []
    sort_patterns: list[list[str]] = []
    paths: list[str] = []
    # Find all files in the directories given as SOURCE
    paths.extend([p for src in args.sources for p in Path(src).rglob("*")])
    for style in styles:
        style_sources = style.get("sources", [])
        paths.extend([Path(p) for src in style_sources for p in glob.glob(src)])
        exclude_patterns.extend(style.get("exclude", []))
        include_patterns.extend(style.get("include", []))
        rna_track_patterns.extend(style.get("rna_tracks", []))
        overlay_rna_strands = style.get("overlay_rna_strands", overlay_rna_strands)
        filetypes_allowed.extend(style.get("filetypes_allowed", []))
        sort_patterns = style.get("tracks_order") or sort_patterns

    exclude_patterns = exclude_patterns or []
    include_patterns = include_patterns or [".*"]
    rna_track_patterns = rna_track_patterns or []
    filetypes_allowed = filetypes_allowed or [".*"]
    filetypes_allowed = map(lambda ft: ft if ft.startswith(".*") else f".*{ft}", filetypes_allowed)
    filetypes_allowed = map(lambda ft: ft if ft.endswith("$") else f"{ft}$", filetypes_allowed)
    filetypes_allowed = list(filetypes_allowed)
    if args.debug:
        print(f"Patterns: ")
        print("sources: ")
        pprint(paths)
        print("\n\ninclude patterns: ")
        pprint(include_patterns)
        print("\n\nexclude patterns: ")
        pprint(exclude_patterns)
        print("\n\nrna track patterns: ")
        pprint(rna_track_patterns)
        print("\n\noverlay rna strands")
        pprint(overlay_rna_strands)
        print("\n\nfiletypes allowed")
        pprint(filetypes_allowed)
        print("\n\nsort patterns")
        pprint(sort_patterns)

    if args.debug:
        print(f"\n\n\nFound {len(paths)} possible tracks")

    trackpaths: list[Track] = group_tracks_with_indexes(paths)
    if args.debug:
        print(f"\n\n\nMerged index files in {len(trackpaths)} possible tracks")

    trackpaths: filter[Track] = list(
        filter(lambda t: match_any_pattern(t, include_patterns + rna_track_patterns), trackpaths)
    )
    if args.debug:
        print(f'\n\n\nTracks which matched an "include" filter: ')
        pprint(trackpaths)
    trackpaths: filter[Track] = list(filter(lambda t: match_any_pattern(t, filetypes_allowed), trackpaths))
    if args.debug:
        print(f"\n\n\nTracks with the correct filetype: ")
        pprint(trackpaths)
    trackpaths: filter[Track] = list(filter(lambda t: not match_any_pattern(t, exclude_patterns), trackpaths))
    if args.debug:
        print(f"\n\n\nTracks after filtering out excluded paths: ")
        pprint(trackpaths)
    if overlay_rna_strands:
        trackpaths: "list[Track | RNABigWig]" = group_rna_strands(trackpaths, rna_track_patterns)
        if args.debug:
            print(f"\n\n\nGrouped RNA trands: ")
            pprint(trackpaths)
    trackpaths: "list[Track | RNABigWig]" = sort_paths(trackpaths, sort_patterns)
    if args.debug:
        print(f"\n\n\nTracks sorted: ")
        pprint(trackpaths)

    ##### Overlay RNA tracks

    tracks = [set_base_track_features(tp, prefix=args.url, order=i + 10) for i, tp in enumerate(trackpaths)]
    for style in styles:
        track_features = style.get("tracks_features", {})
        tracks = [update_track_features(t, track_features) for t in tracks]
    if args.debug:
        print(f"\n\n\nAddition of features to tracks: ")
        pprint(trackpaths)

    json_session = add_to_template(args.template, tracks)
    Path(args.output).write_text(json.dumps(json_session, indent=2))
    return 0
