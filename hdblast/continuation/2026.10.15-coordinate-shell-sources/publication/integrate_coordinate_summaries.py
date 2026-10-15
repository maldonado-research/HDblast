"""Plan, apply or review an explicitly approved coordinate-source summary config.

Default mode is read-only. No scientific result or public execution is assumed.
The only public files this tool may write are the nine fixed summary paths.
It never installs figures, invokes a network command or changes Git state.
"""

import argparse
from collections import Counter
import copy
import hashlib
import html
from html.parser import HTMLParser
import json
from pathlib import Path, PurePosixPath
import re
import subprocess


ROOT = Path(__file__).resolve().parent
TASK = ROOT.parent
MARK = "coordinate-shell-sources-20261015"
CAPSULE = "HDBLAST_COORDINATE_SHELL_SOURCES_20261015"
DELIVERY = "hdblast/continuation/2026.10.15-coordinate-shell-sources"
BASE = "https://github.com/maldonado-research/HDblast/blob/main/research/" + CAPSULE + "/"
TREE = "https://github.com/maldonado-research/HDblast/tree/main/research/" + CAPSULE
EVIDENCE = "https://github.com/maldonado-research/HDblast/blob/main/" + DELIVERY + "/PUBLICATION_RECEIPT.json"
VERSION = "2026.10.15-coordinate-shell-sources (GitHub checkpoint)"
OWNED = {
    "hdblast-worktree": (
        "README.md", "hdblast/README.md", "hdblast/CURRENT_STATUS.md",
        "hdblast/PROGRAM_MAP.md", "site/index.html", "site/es/index.html",
    ),
    "shared-worktree": ("README.md", "index.html", "projects/hdblast/index.html"),
}
REPOSITORIES = {
    "hdblast-worktree": "maldonado-research/HDblast",
    "shared-worktree": "maldonado-research/maldonado-research.github.io",
}
UNRESOLVED = {
    "en": "Full U2 remains quantitatively unresolved.",
    "es": "U2 completo sigue sin resolverse cuantitativamente.",
}
PENDING = {
    "en": "See the separate publication evidence for replay, CI and live-readback status.",
    "es": "Consulte la evidencia de publicación independiente para conocer el estado de la reproducción, la integración continua y la comprobación del sitio publicado.",
}


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def digest(value, label):
    require(isinstance(value, str) and re.fullmatch(r"[0-9a-f]{64}", value),
            "Full lowercase SHA-256 required: " + label)
    return value


def text(value, label):
    require(isinstance(value, str) and value.strip() == value and value
            and "\x00" not in value, "Nonempty, trimmed text required: " + label)
    return value


def canonical(value):
    value = text(value, "relative path")
    path = PurePosixPath(value)
    require(not path.is_absolute() and path.as_posix() == value
            and not any(part in ("", ".", "..") for part in path.parts)
            and re.fullmatch(r"[A-Za-z0-9_./-]+", value),
            "Canonical relative path required: " + value)
    return value


def local_file(root, name):
    path = root / canonical(name)
    require(not path.is_symlink() and path.resolve().is_relative_to(root.resolve())
            and path.is_file(), "Regular file inside the declared root required: " + str(path))
    # A symlink in an intermediate directory could redirect a permitted write.
    require(all(not parent.is_symlink() for parent in path.parents if parent != root.parent),
            "Symlinked file ancestry is not supported: " + str(path))
    return path


def git(wt, *args):
    """Read-only Git commands only; callers supply fixed subcommands."""
    return subprocess.check_output(["git", "--no-optional-locks", "-C", str(wt), *args])


def kind(folder, name):
    if name.endswith(".md"):
        return "markdown"
    if folder == "shared-worktree" and name == "index.html":
        return "card"
    return "detailed"


def language(name):
    return "es" if name == "site/es/index.html" else "en"


class Parser(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.meta, self.links, self.jsonld = {}, [], []
        self.script = None

    def handle_starttag(self, tag, attrs):
        values = dict(attrs)
        if tag == "meta":
            key = values.get("name") or values.get("property")
            self.meta.setdefault(key, []).append(values.get("content"))
        if tag == "a" and "href" in values:
            self.links.append(values["href"])
        if tag == "script" and values.get("type") == "application/ld+json":
            require(self.script is None, "Nested JSON-LD script")
            self.script = []

    def handle_endtag(self, tag):
        if tag == "script" and self.script is not None:
            self.jsonld.append(json.loads("".join(self.script)))
            self.script = None

    def handle_data(self, data):
        if self.script is not None:
            self.script.append(data)


def parsed(source):
    parser = Parser()
    parser.feed(source)
    parser.close()
    require(parser.script is None, "Unclosed JSON-LD script")
    return parser


def unique_object(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, "Duplicate JSON object key: " + key)
        result[key] = value
    return result


def research_node(documents):
    nodes = [node for document in documents for node in document.get("@graph", [])
             if node.get("@type") in ("Dataset", "CreativeWork") and "description" in node]
    require(len(nodes) == 1, "Exactly one described Dataset/CreativeWork JSON-LD node required")
    return nodes[0]


def description_values(parser):
    result = []
    for key in ("description", "og:description", "twitter:description"):
        require(len(parser.meta.get(key, [])) == 1, "Unique description meta field required: " + key)
        result.extend(parser.meta[key])
    result.append(research_node(parser.jsonld)["description"])
    return result


def replace_json_field(source, field, before, after):
    """Replace one pinned JSON string value without reserializing other bytes."""
    candidates = []
    for match in re.finditer(r'"' + re.escape(field) + r'"\s*:\s*("(?:[^"\\]|\\.)*")', source):
        if json.loads(match.group(1)) == before:
            candidates.append(match)
    require(len(candidates) == 1, "Unique pinned JSON-LD field required: " + field)
    match = candidates[0]
    return source[:match.start(1)] + json.dumps(after, ensure_ascii=False) + source[match.end(1):]


def metadata_transform(old, metadata, words, shared):
    require(isinstance(metadata, dict), "Detailed page metadata pins required")
    required = {"old_description", "old_version", "current_line", "old_is_based_on"}
    require(set(metadata) == required, "Detailed metadata keys must be exactly " + repr(sorted(required)))
    old_description = text(metadata["old_description"], "old_description")
    old_version = text(metadata["old_version"], "old_version")
    prior = parsed(old)
    require(description_values(prior) == [old_description] * 4, "Four pinned baseline descriptions required")
    prior_node = research_node(prior.jsonld)
    require(prior_node.get("version") == old_version, "Pinned baseline JSON-LD version differs")
    require(old_version != VERSION, "Baseline version already names this checkpoint")
    heads = list(re.finditer(r"<head\b[^>]*>.*?</head>", old, re.S))
    require(len(heads) == 1, "Unique HTML head required")
    head = heads[0].group()
    updated_head = head
    for key in ("description", "og:description", "twitter:description"):
        candidates = []
        for match in re.finditer(r"<meta\b[^>]*>", updated_head):
            attrs = re.search(r'(?:name|property)="' + re.escape(key) + r'"', match.group())
            if attrs:
                candidates.append(match)
        require(len(candidates) == 1, "Unique literal meta tag required: " + key)
        match = candidates[0]
        content = re.search(r'content="([^"]*)"', match.group())
        require(content and html.unescape(content.group(1)) == old_description,
                "Pinned literal meta content differs: " + key)
        start, end = match.start() + content.start(1), match.start() + content.end(1)
        updated_head = updated_head[:start] + html.escape(words["description"], quote=True) + updated_head[end:]
    updated_head = replace_json_field(updated_head, "description", old_description, words["description"])
    updated_head = replace_json_field(updated_head, "version", old_version, VERSION)
    expected_json = copy.deepcopy(prior.jsonld)
    node = research_node(expected_json)
    node.update({"description": words["description"], "version": VERSION})
    old_source = metadata["old_is_based_on"]
    if old_source is not None:
        text(old_source, "old_is_based_on")
        require(prior_node.get("isBasedOn") == old_source, "Pinned current provenance differs")
        updated_head = replace_json_field(updated_head, "isBasedOn", old_source, TREE)
        node["isBasedOn"] = TREE
    new = old[:heads[0].start()] + updated_head + old[heads[0].end():]
    old_line, new_line = None, None
    if shared:
        old_line = text(metadata["current_line"], "shared current_line")
        match = re.fullmatch(r'(<p\b[^>]*class="version-line"[^>]*>)(.*?)(</p>)', old_line, re.S)
        require(match and old_line.count(old_version) == 1 and old.count(old_line) == 1,
                "Unique full shared current-version paragraph required")
        require(MARK not in old_line and "data-hdb-" in match.group(1),
                "Pinned existing HDBLAST current-version paragraph required")
        new_line = (match.group(1) + html.escape(VERSION + ": " + words["current_summary"]
                    + " " + UNRESOLVED["en"]) + match.group(3))
        new = new.replace(old_line, new_line, 1)
    else:
        require(metadata["current_line"] is None, "Only the shared detailed page has a current_line replacement")
    after = parsed(new)
    require(description_values(after) == [words["description"]] * 4,
            "Four matching replacement descriptions required")
    require(after.jsonld == expected_json, "JSON-LD changed outside the declared current fields")
    outside_head = new[heads[0].start() + len(updated_head):]
    if shared:
        require(outside_head.count(new_line) == 1, "Shared current-version paragraph must be outside the head")
        outside_head = outside_head.replace(new_line, old_line, 1)
    require(outside_head == old[heads[0].end():],
            "Historical HTML body changed outside the declared current-version paragraph")
    return new, {"old_head": head, "new_head": updated_head,
                 "old_current_line": old_line, "new_current_line": new_line}


def links(words, note, mode):
    if mode == "markdown":
        return ("[" + words["proof_label"] + "](" + BASE + note + ") · ["
                + words["checks_label"] + "](" + BASE + "README.md).")
    e = html.escape
    return ('<a href="' + BASE + note + '">' + e(words["proof_label"]) + '</a> · <a href="'
            + BASE + 'README.md">' + e(words["checks_label"]) + "</a>.")


def addition(config, folder, name):
    mode, lang = kind(folder, name), language(name)
    words, note = config["words"][lang], config["result_note"]["path"]
    limits = words["limits"] + " " + UNRESOLVED[lang]
    if mode == "markdown":
        return ("<!-- " + MARK + " -->\n**" + words["title"] + ":** " + words["result"]
                + "\n\n" + limits + "\n\n" + links(words, note, mode) + " " + PENDING[lang]
                + " [" + words["evidence_label"] + "](" + EVIDENCE + ").\n\n")
    e = html.escape
    if mode == "card":
        return ('<p data-hdb-coordinate-sources="' + MARK + '"><strong>' + e(words["title"])
                + ":</strong> " + e(words["card_result"]) + " "
                + e(words["card_limits"] + " " + UNRESOLVED[lang]) + " " + links(words, note, mode)
                + " " + e(PENDING[lang]) + ' <a href="' + EVIDENCE + '">'
                + e(words["evidence_label"]) + "</a>.</p>")
    result = ('<div data-hdb-coordinate-sources="' + MARK + '"><p><strong>' + e(words["title"])
              + ":</strong> " + e(words["result"]) + "</p><p>" + e(limits) + "</p><p>"
              + links(words, note, mode) + "</p><p>" + e(PENDING[lang]) + ' <a href="'
              + EVIDENCE + '">' + e(words["evidence_label"]) + "</a>.</p>")
    if config["figure"] is not None:
        figure = config["figure"]
        src = ("../assets/" if lang == "es" else "assets/") + figure["asset_name"]
        result += ('<figure><a href="' + BASE + figure["path"] + '"><img src="' + src
                   + '" alt="' + e(words["alt"], quote=True)
                   + '" loading="lazy" style="max-width:100%;height:auto"></a><figcaption>'
                   + e(words["caption"]) + "</figcaption></figure>")
    return result + "</div>\n\n"


def transform(config, row, old):
    folder, name = row["repository"], row["path"]
    mode = kind(folder, name)
    require(MARK not in old, "New checkpoint marker already occurs in baseline: " + name)
    anchor = text(row["anchor"], "exact insertion anchor")
    require(old.count(anchor) == 1, "Insertion anchor must occur exactly once: " + name)
    require("<" in anchor and MARK not in anchor, "Existing literal checkpoint anchor required: " + name)
    if mode == "markdown":
        require(re.fullmatch(r"<!-- [a-z0-9-]+ -->", anchor) is not None,
                "Exact historical Markdown marker required")
    else:
        tag = "p" if mode == "card" else "div"
        require(re.fullmatch(r"<" + tag + r' data-hdb-[a-z-]+="[a-z0-9-]+">', anchor) is not None,
                "Exact historical HDBLAST block opening required: " + name)
        body_start = old.find("<body")
        require(body_start >= 0 and old.index(anchor) > body_start, "Insertion anchor must be in the body")
    detail = {}
    if mode == "detailed":
        changed, detail = metadata_transform(old, row["metadata"], config["words"][language(name)],
                                             folder == "shared-worktree")
    else:
        require(row["metadata"] is None, "Markdown and card entries cannot replace metadata")
        changed = old
    added = addition(config, folder, name)
    require(added.count(MARK) == 1 and anchor not in added, "Rendered block contains an ambiguous marker or anchor")
    new = changed.replace(anchor, added + anchor, 1)
    detail["addition"] = added
    require(new.count(MARK) == 1, "Exactly one new checkpoint marker required: " + name)
    require(reverse(new, detail) == old, "Exact baseline-byte reversal failed: " + name)
    if mode != "markdown":
        before, after = parsed(old), parsed(new)
        require(not (Counter(before.links) - Counter(after.links)), "Historical link occurrence removed: " + name)
        require(re.findall(r"<nav\b[^>]*>.*?</nav>", old, re.S)
                == re.findall(r"<nav\b[^>]*>.*?</nav>", new, re.S), "Navigation changed: " + name)
        if mode == "card":
            pattern = r'<article\b[^>]*class="project-card"[^>]*>.*?</article>'
            prior_cards, new_cards = re.findall(pattern, old, re.S), re.findall(pattern, new, re.S)
            require(prior_cards and len(prior_cards) == len(new_cards), "Project-card count changed")
            target = [i for i, card in enumerate(prior_cards) if anchor in card
                      and re.search(r'<h3><a href="\./projects/hdblast/">HDBLAST</a></h3>', card)]
            require(len(target) == 1, "Unique HDBLAST project card containing the anchor required")
            require(all(a == b if i != target[0] else b.replace(added, "", 1) == a
                        for i, (a, b) in enumerate(zip(prior_cards, new_cards))),
                    "An existing project card changed outside the HDBLAST insertion")
            detail["other_project_cards_preserved"] = len(prior_cards) - 1
    return new, detail


def reverse(new, detail):
    added = detail["addition"]
    require(new.count(added) == 1, "Unique rendered addition required for reversal")
    restored = new.replace(added, "", 1)
    if "new_head" in detail:
        require(restored.count(detail["new_head"]) == 1, "Unique updated HTML head required for reversal")
        restored = restored.replace(detail["new_head"], detail["old_head"], 1)
        if detail["new_current_line"] is not None:
            require(restored.count(detail["new_current_line"]) == 1,
                    "Unique updated current line required for reversal")
            restored = restored.replace(detail["new_current_line"], detail["old_current_line"], 1)
    return restored


def load_config(path, approved_sha256):
    raw = path.read_bytes()
    require(sha(raw) == digest(approved_sha256, "--approved-config-sha256"),
            "Config bytes differ from the explicitly supplied approved SHA-256")
    config = json.loads(raw, object_pairs_hook=unique_object)
    expected_keys = {"schema_version", "checkpoint", "capsule", "delivery", "approval", "scientific_scope",
                     "repositories", "result_note", "result_inputs", "words", "figure", "summaries"}
    require(set(config) == expected_keys and type(config["schema_version"]) is int
            and config["schema_version"] == 1, "Unsupported or incomplete config schema")
    require((config["checkpoint"], config["capsule"], config["delivery"]) == (MARK, CAPSULE, DELIVERY),
            "Config checkpoint/capsule/delivery differs from this helper")
    approval = config["approval"]
    require(isinstance(approval, dict) and set(approval) == {"approved", "approved_by", "basis"}
            and approval["approved"] is True, "An explicit approved config is required")
    text(approval["approved_by"], "approval.approved_by")
    text(approval["basis"], "approval.basis")
    required_scope = {
        "boundary_trace_and_proper_time_derivative": True,
        "one_copy_coordinate_scalar_and_isotropic_israel_sources": True,
        "scalar_degrees": [0, 2, 4],
        "full_u2_quantitatively_unresolved": True,
    }
    require(json.dumps(config["scientific_scope"], sort_keys=True) == json.dumps(required_scope, sort_keys=True),
            "Approved scientific scope must explicitly retain the unresolved full-U2 limit")
    require(set(config["repositories"]) == set(OWNED), "Exactly the two repository slots are required")
    for folder, pin in config["repositories"].items():
        require(set(pin) == {"repository", "worktree", "baseline_commit"}, "Unexpected repository config keys")
        require(pin["repository"] == REPOSITORIES[folder], "Repository identity differs: " + folder)
        require(pin["worktree"] == str(ROOT / folder), "Worktree must be this task's isolated publication slot: " + folder)
        require(isinstance(pin["baseline_commit"], str)
                and re.fullmatch(r"[0-9a-f]{40}", pin["baseline_commit"]), "Full immutable baseline commit required")
    require(isinstance(config["result_inputs"], list), "result_inputs must be a list")
    inputs = [config["result_note"]] + config["result_inputs"]
    names = set()
    input_pins = []
    for item in inputs:
        require(isinstance(item, dict) and set(item) == {"path", "sha256"}, "Each result input needs path/sha256")
        name = canonical(item["path"])
        require(name not in names, "Duplicate result input: " + name)
        names.add(name)
        data = local_file(TASK, name).read_bytes()
        require(sha(data) == digest(item["sha256"], name), "Result input changed: " + name)
        input_pins.append({"path": name, "bytes": len(data), "sha256": sha(data)})
    require(config["result_note"]["path"].endswith(".md"), "Authoritative result note must be Markdown")
    require(set(config["words"]) == {"en", "es"}, "Approved English and Spanish wording required")
    common = {"title", "result", "limits", "description", "proof_label", "checks_label", "evidence_label"}
    for lang in ("en", "es"):
        required = common | ({"card_result", "card_limits", "current_summary"} if lang == "en" else set())
        if config["figure"] is not None:
            required |= {"alt", "caption"}
        require(set(config["words"][lang]) == required, "Wording keys differ for " + lang)
        for key, value in config["words"][lang].items():
            text(value, "words." + lang + "." + key)
            require(MARK not in value and re.search(r"</?[A-Za-z][^>]*>|<![^>]*>|<\?[^>]*>", value) is None
                    and "\n" not in value,
                    "Wording must be plain single-paragraph text without markers/HTML: " + key)
            if key.endswith("label") or key == "title":
                require(not any(char in value for char in "[]*"), "Markdown delimiters are not allowed in labels/title")
    rows = config["summaries"]
    require(isinstance(rows, list) and len(rows) == 9, "Exactly nine summary rows required")
    keys = []
    for row in rows:
        require(isinstance(row, dict) and set(row) == {"repository", "path", "baseline_sha256", "anchor", "metadata"},
                "Unexpected or incomplete summary row")
        keys.append((row["repository"], row["path"]))
        digest(row["baseline_sha256"], row["path"])
    require(len(set(keys)) == 9 and set(keys) == {(folder, name) for folder, names in OWNED.items() for name in names},
            "Summary rows must cover only the nine fixed paths, without duplicates")
    return config, {"path": str(path.resolve()), "bytes": len(raw), "sha256": sha(raw)}, input_pins


def figure_receipt(config):
    figure = config["figure"]
    if figure is None:
        return {"status": "PASS", "mode": "no_figure_assets", "source": None,
                "site_assets": [], "all_asset_bytes_equal_to_source": True}
    require(isinstance(figure, dict) and set(figure) == {"path", "sha256", "asset_name"},
            "Figure needs exactly path/sha256/asset_name")
    source = local_file(TASK, figure["path"]).read_bytes()
    require(sha(source) == digest(figure["sha256"], "figure source"), "Pinned authoritative figure changed")
    require(source.startswith(b"\x89PNG\r\n\x1a\n"), "Selected figure must be a PNG")
    name = figure["asset_name"]
    require(isinstance(name, str) and re.fullmatch(r"[a-z0-9_-]+\.png", name), "Safe PNG asset filename required")
    assets = []
    for folder, prefix, url in (
        ("hdblast-worktree", "site/assets/", "https://maldonado-research.github.io/HDblast/assets/"),
        ("shared-worktree", "projects/hdblast/assets/", "https://maldonado-research.github.io/projects/hdblast/assets/"),
    ):
        asset = local_file(Path(config["repositories"][folder]["worktree"]), prefix + name).read_bytes()
        require(asset == source, "Root must first install exact authoritative figure bytes: " + folder)
        assets.append({"repository": REPOSITORIES[folder], "path": prefix + name, "url": url + name,
                       "bytes": len(asset), "sha256": sha(asset)})
    return {"status": "PASS", "mode": "figure_assets", "source": {
        "path": figure["path"], "bytes": len(source), "sha256": sha(source)},
        "site_assets": assets, "all_asset_bytes_equal_to_source": True}


def repository_observation(folder, pin):
    wt, base = Path(pin["worktree"]), pin["baseline_commit"]
    require(wt.is_dir() and not wt.is_symlink(), "Isolated worktree does not exist: " + str(wt))
    require(git(wt, "rev-parse", "--show-toplevel").decode().strip() == str(wt), "Git root differs from worktree slot")
    require(git(wt, "rev-parse", base + "^{commit}").decode().strip() == base,
            "Baseline commit is unavailable or differs")
    require(subprocess.run(["git", "--no-optional-locks", "-C", str(wt), "merge-base", "--is-ancestor", base, "HEAD"],
                           check=False, stdout=subprocess.DEVNULL).returncode == 0,
            "Pinned baseline is not an ancestor of current HEAD")
    old_names = {name.decode() for name in git(wt, "ls-tree", "-rz", "--name-only", base).split(b"\0") if name}
    changed = {name.decode() for name in git(wt, "diff", "--no-ext-diff", "--name-only", "-z", base, "--").split(b"\0") if name}
    require(set(OWNED[folder]) <= old_names, "An owned summary is absent from the immutable baseline")
    require(not ((changed & old_names) - set(OWNED[folder])),
            "An existing tracked file changed outside the nine summaries: " + folder)
    return {"repository": REPOSITORIES[folder], "worktree": str(wt), "base": base,
            "observed_head": git(wt, "rev-parse", "HEAD").decode().strip(),
            "baseline_tree": git(wt, "rev-parse", base + "^{tree}").decode().strip(),
            "other_existing_files_preserved": len(old_names) - len(OWNED[folder]),
            "existing_file_check": "git diff pinned immutable baseline against current index and worktree",
            "new_files_outside_summaries": "not authenticated by this helper"}


def prepare(config, review_only=False):
    observations = [repository_observation(folder, config["repositories"][folder]) for folder in OWNED]
    figures = figure_receipt(config)
    planned, records = [], []
    for row in config["summaries"]:
        folder, name = row["repository"], row["path"]
        pin = config["repositories"][folder]
        wt = Path(pin["worktree"])
        entry = git(wt, "ls-tree", pin["baseline_commit"], "--", name).decode().split()
        require(entry and entry[0] in ("100644", "100755"), "Baseline summary must be a regular Git blob")
        old_bytes = git(wt, "show", pin["baseline_commit"] + ":" + name)
        require(sha(old_bytes) == row["baseline_sha256"], "Baseline summary hash differs: " + name)
        new, detail = transform(config, row, old_bytes.decode("utf-8"))
        new_bytes = new.encode("utf-8")
        path = local_file(wt, name)
        require(bool(path.stat().st_mode & 0o111) == (entry[0] == "100755"),
                "Existing summary executable mode differs from the baseline: " + name)
        indexed = git(wt, "ls-files", "--stage", "--", name).decode().split()
        require(indexed and indexed[0] == entry[0] and indexed[2] == "0",
                "Existing summary index mode or merge stage differs from the baseline: " + name)
        current = path.read_bytes()
        require(current in (old_bytes, new_bytes), "Independent existing summary edits must be preserved: " + name)
        require(not review_only or current == new_bytes, "Review requires exact expected integrated bytes: " + name)
        reverted = reverse(new, detail).encode("utf-8")
        require(reverted == old_bytes, "Exact byte reversal failed: " + name)
        planned.append((path, current, new_bytes))
        record = {"repository": REPOSITORIES[folder], "slot": folder, "path": name,
                  "baseline_commit": pin["baseline_commit"], "baseline_sha256": sha(old_bytes),
                  "observed_sha256": sha(current), "expected_updated_sha256": sha(new_bytes),
                  "updated_sha256": sha(new_bytes) if current == new_bytes else None,
                  "reversed_sha256": sha(reverted), "exact_byte_reversal_passed": True,
                  "insertion_anchor": row["anchor"], "new_checkpoint_markers": 1,
                  "all_prior_bytes_preserved": kind(folder, name) != "detailed",
                  "historical_body_preserved_except_declared_current_line": kind(folder, name) == "detailed"}
        if name.endswith(".html"):
            record["historical_link_occurrences_and_navigation_preserved"] = True
        if "other_project_cards_preserved" in detail:
            record["other_project_cards_preserved"] = detail["other_project_cards_preserved"]
            record["only_hdblast_card_modified"] = True
        if kind(folder, name) == "detailed":
            record.update({"description_fields_consistent": 4, "current_jsonld_version": VERSION,
                           "jsonld_changes_confined_to_declared_current_fields": True,
                           "shared_current_line_replaced": folder == "shared-worktree"})
        records.append(record)
    return planned, observations, figures, records


def write_receipt(name, value):
    path = ROOT / name
    require(not path.is_symlink(), "Refusing symlinked receipt output: " + name)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True, help="Root-reviewed JSON config")
    parser.add_argument("--approved-config-sha256", required=True, help="SHA-256 of the exact approved config bytes")
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--apply", action="store_true", help="Write the nine planned summaries and local receipts")
    modes.add_argument("--review-only", action="store_true", help="Require integrated bytes and write local receipts")
    args = parser.parse_args()
    config, config_pin, input_pins = load_config(args.config, args.approved_config_sha256)
    planned, observations, figures, records = prepare(config, args.review_only)
    mode = "review_only" if args.review_only else "apply" if args.apply else "plan_only"
    if args.apply:
        # Validate every target before the first write; never overwrite independent edits.
        for path, observed, _ in planned:
            require(path.read_bytes() == observed, "Summary changed after preflight: " + str(path))
        for path, _, expected in planned:
            path.write_bytes(expected)
        _, observations, figures, records = prepare(config, review_only=True)
    integrated = all(path.read_bytes() == expected for path, _, expected in planned)
    receipt = {"status": "PASS" if mode != "plan_only" else "PASS_READ_ONLY_PLAN",
               "checkpoint": MARK, "capsule": CAPSULE, "delivery": DELIVERY, "mode": mode,
               "config": config_pin, "approval": config["approval"], "scientific_scope": config["scientific_scope"],
               "integration_helper": {"path": str(Path(__file__).resolve()),
                                      "sha256": sha(Path(__file__).read_bytes())},
               "base_commits": {folder: pin["baseline_commit"] for folder, pin in config["repositories"].items()},
               "changed_owned_summaries": records, "preservation": observations,
               "result_input_pins": input_pins, "wording": config["words"],
               "all_nine_current_files_equal_expected_integration": integrated,
               "figure_assets_installed_by_this_tool": False,
               "figure_integration": figures, "publication_evidence_status": "PENDING_FINAL_ADDITIVE_RECEIPT",
               "observed_publication": {"network_accessed": False, "git_mutated": False,
                                        "pr_ci_merge_deployment_or_live_readback_checked": False},
               "scope": "Local immutable-baseline content integration and preservation only. This receipt does not establish scientific approval, source-capsule authenticity, CI, a merge, deployment or public readback beyond the recorded config assertion and file checks."}
    if mode != "plan_only":
        for name in ("SITE_INTEGRATION_REVIEW.json", "SUMMARY_INTEGRATION_PREFLIGHT.json"):
            write_receipt(name, receipt)
        write_receipt("FIGURE_INTEGRATION_RECEIPT.json", figures)
        write_receipt("SITE_BASELINES.json", receipt["base_commits"])
    print(json.dumps(receipt, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
