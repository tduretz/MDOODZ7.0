#!/usr/bin/env python3
"""Build a printable PDF of the defence document.

Invocation
----------

    python3 build_pdf.py [defence.md]

Output: ``defence.pdf`` (and a side-channel ``defence.html``) next to
the input markdown. Page budget target (D4): 25 pages at A4 11 pt.

Rendering pipeline
------------------

1. Read the markdown, convert to HTML with the stdlib-adjacent
   ``markdown`` package (fenced_code + tables + toc extensions).
2. Wrap in a primer-series-palette HTML shell with print-oriented CSS.
3. Invoke Google Chrome in headless mode to render HTML → PDF.

Chrome is used because it's already installed on the workstation (see
``design.md`` §D7 — "symlink to the reading-library builder, or a thin
copy"; this is the thin copy). The alternatives — ``pandoc`` and
``weasyprint`` — are not required to be present on every machine that
ships the defence, so relying on Chrome keeps the pipeline portable
across Roman's MacBook and the AWS-provisioned workstation.

Images are resolved relative to the markdown file. The ``figs/``
subdirectory (containing ``defence_*.png`` and ``vcycle_*.gif, .png``)
must sit next to the markdown or the paths in the rendered HTML will
404.
"""
from __future__ import annotations

import pathlib
import shutil
import subprocess
import sys

import markdown

CHROME_CANDIDATES = [
    "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
    "/Applications/Chromium.app/Contents/MacOS/Chromium",
    "/usr/bin/google-chrome",
    "/usr/bin/chromium",
    "/usr/bin/chromium-browser",
]


CSS = """
@page {
    size: A4;
    margin: 2.0cm 1.8cm 2.0cm 1.8cm;
}
html, body {
    background: #fbfaf6;
    color: #1f2a37;
    font-family: "Charter", "Georgia", "DejaVu Serif", "Liberation Serif", serif;
    font-size: 10.5pt;
    line-height: 1.38;
    margin: 0;
    padding: 0;
}
body {
    max-width: 100%;
}
h1, h2, h3, h4 {
    font-family: "Helvetica Neue", "Arial", sans-serif;
    color: #1f2a37;
    page-break-after: avoid;
    line-height: 1.22;
}
h1 {
    font-size: 22pt;
    margin-top: 0.2em;
    border-bottom: 2px solid #1f2a37;
    padding-bottom: 0.2em;
}
h2 {
    font-size: 15pt;
    margin-top: 1.4em;
    border-bottom: 1px solid #d9d3c5;
    padding-bottom: 0.15em;
}
h3 { font-size: 12pt; margin-top: 1.1em; }
h4 { font-size: 11pt; margin-top: 0.9em; color: #4b3b8f; }
p { text-align: justify; }
em, i { color: #4b3b8f; }
strong, b { color: #1f2a37; }
code, pre {
    font-family: "Menlo", "DejaVu Sans Mono", "Liberation Mono", monospace;
    font-size: 9.5pt;
    color: #1f2a37;
}
pre {
    background: #f2efe6;
    border-left: 3px solid #4b3b8f;
    padding: 0.6em 0.9em;
    overflow-x: hidden;
    page-break-inside: avoid;
}
code {
    background: #f2efe6;
    padding: 0 0.25em;
    border-radius: 2px;
}
table {
    border-collapse: collapse;
    width: 100%;
    margin: 0.8em 0;
    font-size: 9.5pt;
    page-break-inside: avoid;
}
th, td {
    border: 1px solid #d9d3c5;
    padding: 4px 7px;
    text-align: left;
    vertical-align: top;
}
th {
    background: #ece8dc;
    font-family: "Helvetica Neue", "Arial", sans-serif;
    font-weight: 600;
}
img {
    max-width: 85%;
    display: block;
    margin: 0.6em auto;
    page-break-inside: avoid;
}
blockquote {
    border-left: 3px solid #b04b3b;
    color: #3b3b3b;
    margin: 0.7em 0;
    padding: 0.3em 0.9em;
    background: #f7f2e9;
}
hr {
    border: 0;
    border-top: 1px solid #d9d3c5;
    margin: 1.2em 0;
}
ul, ol { margin: 0.4em 0 0.7em 1.3em; }
li { margin: 0.1em 0; }
a { color: #3b7d8a; text-decoration: none; }
.page-break { page-break-after: always; }
"""


HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{title}</title>
<style>{css}</style>
</head>
<body>
{body}
</body>
</html>
"""


def find_chrome() -> str:
    for path in CHROME_CANDIDATES:
        if pathlib.Path(path).exists():
            return path
    which = shutil.which("google-chrome") or shutil.which("chromium") \
        or shutil.which("chromium-browser")
    if which:
        return which
    raise RuntimeError(
        "Could not find a Chrome/Chromium binary — checked "
        f"{CHROME_CANDIDATES} and $PATH"
    )


def render_markdown(md_path: pathlib.Path) -> str:
    src = md_path.read_text(encoding="utf-8")
    html_body = markdown.markdown(
        src,
        extensions=["fenced_code", "tables", "toc", "sane_lists"],
        output_format="html5",
    )
    return HTML_TEMPLATE.format(title=md_path.stem, css=CSS, body=html_body)


def main(argv: list[str]) -> int:
    here = pathlib.Path(__file__).resolve().parent
    md_arg = argv[1] if len(argv) > 1 else "defence.md"
    md_path = (here / md_arg).resolve() if not pathlib.Path(md_arg).is_absolute() \
        else pathlib.Path(md_arg).resolve()
    if not md_path.exists():
        print(f"error: {md_path} does not exist", file=sys.stderr)
        return 2

    out_html = md_path.with_suffix(".html")
    out_pdf = md_path.with_suffix(".pdf")

    out_html.write_text(render_markdown(md_path), encoding="utf-8")
    print(f"wrote {out_html.relative_to(here)}")

    chrome = find_chrome()
    cmd = [
        chrome,
        "--headless=new",
        "--disable-gpu",
        "--no-sandbox",
        "--no-pdf-header-footer",
        f"--print-to-pdf={out_pdf}",
        out_html.as_uri(),
    ]
    print(" ".join(f'"{c}"' if " " in c else c for c in cmd))
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        print(proc.stdout)
        print(proc.stderr, file=sys.stderr)
        return proc.returncode
    print(f"wrote {out_pdf.relative_to(here)}")

    size_kb = out_pdf.stat().st_size // 1024
    print(f"defence.pdf: {size_kb} KB")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
