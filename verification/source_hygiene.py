#!/usr/bin/env python3
"""Static hygiene checks for authoritative LaTeX manuscripts.

This check detects accidental conversation/template material and broken
cross-references. It does not certify mathematical correctness, and it does
not treat truthful authorship or assistance disclosures as errors.
"""
import collections
import pathlib
import re
import sys

unsafe = {
    "role delimiters": re.compile(r"(?im)^\s*(?:SYSTEM|DEVELOPER|ASSISTANT|USER)\s*:\s"),
    "chat markup": re.compile(r"<\|(?:im_start|im_end|endoftext)\|>|\[INST\]|\[/INST\]"),
    "internal prompts": re.compile(r"(?i)\b(?:system prompt|developer instruction|hidden chain.of.thought|private scratchpad|prompt injection)\b"),
    "assistant filler": re.compile(r"(?i)\b(?:as an AI language model|here is the requested (?:proof|paper)|I have now checked your|copy and paste this prompt)\b"),
    "unfinished placeholders": re.compile(r"(?i)\b(?:TODO|FIXME|TBD|PLACEHOLDER)\b"),
}

def unescaped_tex_comments(src):
    lines = []
    for line in src.splitlines():
        for i, char in enumerate(line):
            if char == "%" and (i == 0 or line[i - 1] != "\\"):
                line = line[:i]
                break
        lines.append(line)
    return "\n".join(lines)

def inspect(path):
    raw = pathlib.Path(path).read_text(encoding="utf-8")
    errors = []
    for kind, pattern in unsafe.items():
        for m in pattern.finditer(raw):
            errors.append(f"{path}:{raw.count(chr(10), 0, m.start()) + 1}: {kind}")
    code = unescaped_tex_comments(raw)
    labels = re.findall(r"\\label(?:\[[^\]]*\])?\{([^}]*)\}", code)
    for key, count in collections.Counter(labels).items():
        if count > 1:
            errors.append(f"{path}: duplicate label {key}")
    bibkeys = set(re.findall(r"\\bibitem(?:\[[^\]]*\])?\{([^}]*)\}", code))
    for m in re.finditer(r"\\cite(?:\[[^\]]*\]){0,2}\{([^}]*)\}", code, flags=re.S):
        for key in m.group(1).split(","):
            if key.strip() not in bibkeys:
                errors.append(f"{path}: undefined bibliographic key {key.strip()}")
    for m in re.finditer(r"\\(?:eqref|cref|Cref|ref|autoref)\{([^}]*)\}", code, flags=re.S):
        for key in m.group(1).split(","):
            if key.strip() not in labels:
                errors.append(f"{path}: undefined cross-reference {key.strip()}")
    if "\\documentclass" not in raw or "\\end{document}" not in raw:
        errors.append(f"{path}: incomplete standalone LaTeX manuscript")
    print(f"checked {path}: {len(raw.splitlines())} lines, {len(labels)} labels, {len(bibkeys)} bibliography entries")
    return errors

if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise SystemExit("usage: source_hygiene.py FILE.tex [FILE.tex ...]")
    all_errors = []
    for name in sys.argv[1:]:
        all_errors.extend(inspect(name))
    if all_errors:
        print("\n".join(all_errors), file=sys.stderr)
        raise SystemExit(1)
    print("PASS: no matched conversation artefacts, placeholders, or broken source references")
