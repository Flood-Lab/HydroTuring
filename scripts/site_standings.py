#!/usr/bin/env python3
"""Write the per-law standings the site's two charts are drawn from, and the
README's copies of those charts.

The site shows the merged probes as a pie split by conservation law, and
each evaluated model as one bar split the same way: how many probes under
each law could score the model, and how many of those it passed. Both come
from the repository. A probe's law is the one its probe.yaml declares, and
a model's standing is the one tests/test_docs_in_sync.py checks the tables
against: the rows archived at the model's newest version in
models/result.csv, with N/A left out of both numbers.

The numbers live in a JSON block inside site/index.html, so the page stays
one self-contained file. The README cannot run the site's script, so the
same numbers are drawn into SVG files under res/, in the site's colours and
with its English labels, once for a light and once for a dark theme. The
README's flowchart is the site's "How it works" flowchart, copied out of
site/index.html with the theme's colours filled in. Run this after archiving
rows, merging a probe or editing the site's flowchart:

    python3 scripts/site_standings.py          # rewrite the block and the images
    python3 scripts/site_standings.py --check  # exit 1 if either is stale

Nothing here needs a dependency.
"""

from __future__ import annotations

import csv
import html as htmllib
import json
import math
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SITE = ROOT / "site" / "index.html"
LAWS = ("mass", "energy", "momentum")
BLOCK = re.compile(
    r'(<script type="application/json" id="ht-standings">)(.*?)(</script>)', re.DOTALL
)


LAW_LINE = re.compile(r"^law:\s*['\"]?(\w+)", re.MULTILINE)
ID_LINE = re.compile(r"^id:\s*['\"]?([\w./-]+)", re.MULTILINE)


def probe_laws() -> dict[str, str]:
    """Every merged probe's id, mapped to the law its probe.yaml declares.

    Read from the `law:` and `id:` lines rather than the directory, because a
    probe's id need not start with its law (`real/<slug>` is allowed).
    """
    laws = {}
    for path in ROOT.glob("probes/*/*/probe.yaml"):
        text = path.read_text(encoding="utf-8")
        law, probe_id = LAW_LINE.search(text), ID_LINE.search(text)
        if not law or not probe_id:
            raise SystemExit(f"{path} has no top-level law: or id: line")
        if law.group(1) not in LAWS:
            raise SystemExit(f"{path} declares law {law.group(1)!r}, not one of {LAWS}")
        laws[probe_id.group(1)] = law.group(1)
    return laws


def evaluated_models() -> list[str]:
    """The models the site lists: everything but the reference models and the template."""
    return sorted(
        path.parent.name
        for path in ROOT.glob("models/*/model.yaml")
        if not path.parent.name.startswith(("_", "reference_"))
    )


def standings() -> dict:
    laws = probe_laws()
    with open(ROOT / "models" / "result.csv", newline="", encoding="utf-8") as fh:
        rows = [row for row in csv.DictReader(fh) if row["probe"] in laws]
    version = {row["model"]: row["version"] for row in rows}
    latest = {
        (row["model"], row["probe"]): row["verdict"]
        for row in rows
        if row["version"] == version[row["model"]]
    }
    models = {}
    for model in evaluated_models():
        if model not in version:
            continue
        counts = {law: [0, 0] for law in LAWS}
        for (name, probe), verdict in latest.items():
            if name != model or verdict == "N/A":
                continue
            counts[laws[probe]][1] += 1
            counts[laws[probe]][0] += verdict == "PASS"
        models[model] = counts
    probes = {law: sum(1 for value in laws.values() if value == law) for law in LAWS}
    return {"probes": probes, "models": models}


def render(data: dict) -> str:
    """One model per line, so a re-archived model shows up as a one-line diff."""
    lines = [f'{{"probes":{json.dumps(data["probes"], separators=(",", ":"))},', '"models":{']
    items = list(data["models"].items())
    for i, (model, counts) in enumerate(items):
        comma = "," if i < len(items) - 1 else ""
        lines.append(f'"{model}":{json.dumps(counts, separators=(",", ":"))}{comma}')
    lines.append("}}")
    return "\n" + "\n".join(lines) + "\n"


def embedded(html: str) -> dict:
    match = BLOCK.search(html)
    if not match:
        raise SystemExit(f'{SITE} has no <script type="application/json" id="ht-standings"> block')
    return json.loads(match.group(2))


# ---------------------------------------------------------------------------
# The README's images. Everything they say is read from site/index.html: the
# colours from its :root tokens, the labels from its English strings, the model
# order from its English models table, and the flowchart from its markup.
# ---------------------------------------------------------------------------

THEMES = ("light", "dark")
SANS = '"IBM Plex Sans","Helvetica Neue",Arial,sans-serif'
CONDENSED = '"IBM Plex Sans Condensed","IBM Plex Sans","Helvetica Neue",Arial,sans-serif'
MONO = '"IBM Plex Mono",ui-monospace,SFMono-Regular,Menlo,monospace'


def image_paths() -> dict[str, Path]:
    """Each README image, keyed by chart and theme, e.g. "pie-dark"."""
    paths = {}
    for theme in THEMES:
        suffix = "" if theme == "light" else "-dark"
        paths[f"pie-{theme}"] = ROOT / "res" / f"probes-by-law{suffix}.svg"
        paths[f"bars-{theme}"] = ROOT / "res" / f"models-by-law{suffix}.svg"
        paths[f"infra-{theme}"] = ROOT / "res" / f"hydroturing-infra{suffix}.svg"
    return paths


def palette(html: str, theme: str) -> dict[str, str]:
    """The site's colour tokens for one theme, without the leading `--`."""
    pattern = r":root\{(.*?)\}" if theme == "light" else r':root\[data-theme="dark"\]\{(.*?)\}'
    match = re.search(pattern, html, re.DOTALL)
    if not match:
        raise SystemExit(f"{SITE} has no {theme} :root colour block")
    return {
        name: value.strip() for name, value in re.findall(r"--([\w-]+):\s*([^;]+);", match.group(1))
    }


def english(html: str, key: str) -> str:
    """One string from the site's English dictionary, the first of the three."""
    match = re.search(r'"' + re.escape(key) + r'":"((?:[^"\\]|\\.)*)"', html)
    if not match:
        raise SystemExit(f"{SITE} has no English string {key!r}")
    return json.loads('"' + match.group(1) + '"')


def chart_order(html: str, data: dict) -> list[str]:
    """The site's bar order: its English models table first, then anything it leaves out."""
    match = re.search(r'"models\.rows":`(.*?)`', html, re.DOTALL)
    listed = re.findall(r'<td class="mono">([^<]+)</td>', match.group(1)) if match else []
    listed = [name.strip() for name in listed]
    return [m for m in listed if m in data["models"]] + [
        m for m in data["models"] if m not in listed
    ]


def num(x: float) -> str:
    """A coordinate with at most two decimals, so the output is stable across runs."""
    text = f"{x:.2f}".rstrip("0").rstrip(".")
    return "0" if text == "-0" else text


def esc(text: str) -> str:
    return htmllib.escape(text, quote=True)


def ink_on(fill: str) -> str:
    """White or black, whichever clears the fill better; the site's inkOn()."""
    rgb = [int(fill.lstrip("#")[i : i + 2], 16) / 255 for i in (0, 2, 4)]
    lin = [c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4 for c in rgb]
    lum = 0.2126 * lin[0] + 0.7152 * lin[1] + 0.0722 * lin[2]
    return "#FFFFFF" if 1.05 / (lum + 0.05) >= (lum + 0.05) / 0.05 else "#000000"


def text_width(text: str, size: float, factor: float = 0.56) -> float:
    """A generous estimate of a proportional label's width, for spacing only."""
    return len(text) * size * factor


def card(width: float, height: float, c: dict[str, str]) -> str:
    """The site's <figure>: a surface with a hairline border and rounded corners."""
    return (
        f'<rect x="0.5" y="0.5" width="{num(width - 1)}" height="{num(height - 1)}" rx="10" '
        f'fill="{c["surface"]}" stroke="{c["line"]}"/>'
    )


def svg_open(width: float, height: float, title: str) -> str:
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{num(width)}" height="{num(height)}" '
        f'viewBox="0 0 {num(width)} {num(height)}" role="img">\n<title>{esc(title)}</title>'
    )


def percentages(counts: list[int]) -> list[int]:
    """Whole percentages by largest remainder, so the legend always sums to 100."""
    total = sum(counts)
    exact = [n / total * 100 for n in counts]
    pcts = [math.floor(x) for x in exact]
    order = sorted(range(len(counts)), key=lambda i: exact[i] - pcts[i], reverse=True)
    for i in order[: 100 - sum(pcts)]:
        pcts[i] += 1
    return pcts


def pie_svg(html: str, data: dict, theme: str) -> str:
    c = palette(html, theme)
    counts = [data["probes"][law] for law in LAWS]
    total = sum(counts)
    pcts = percentages(counts)
    names = [english(html, f"law.{law}") for law in LAWS]
    title_text = english(html, "chart.pie.title")
    width, pad_x, top = 470, 20, 22
    pie_top = top + 24 + 14
    height = pie_top + 200 + 14
    summary = ", ".join(
        f"{n} {count} ({p}%)" for n, count, p in zip(names, counts, pcts, strict=True)
    )
    out = [svg_open(width, height, f"{title_text}: {summary}"), card(width, height, c)]
    out.append(
        f'<text x="{pad_x}" y="{top + 17}" font-family=\'{CONDENSED}\' font-size="15" '
        f'font-weight="600" fill="{c["ink"]}">{esc(title_text)}</text>'
    )
    cx, cy, r = pad_x + 100, pie_top + 100, 96
    angle = -math.pi / 2
    labels = []
    for law, n in zip(LAWS, counts, strict=True):
        if not n:
            continue
        sweep = n / total * 2 * math.pi
        a0, a1 = angle, angle + sweep
        angle = a1
        fill = c[f"law-{law}"]
        if n == total:
            d = f"M{num(cx)} {num(cy - r)} A{r} {r} 0 1 1 {num(cx)} {num(cy + r)} A{r} {r} 0 1 1 {num(cx)} {num(cy - r)} Z"
        else:
            x0, y0 = cx + r * math.cos(a0), cy + r * math.sin(a0)
            x1, y1 = cx + r * math.cos(a1), cy + r * math.sin(a1)
            large = 1 if sweep > math.pi else 0
            d = f"M{num(cx)} {num(cy)} L{num(x0)} {num(y0)} A{r} {r} 0 {large} 1 {num(x1)} {num(y1)} Z"
        out.append(f'<path d="{d}" fill="{fill}" stroke="{c["surface"]}" stroke-width="2"/>')
        if sweep > 0.3:
            mid = (a0 + a1) / 2
            tx, ty = cx + r * 0.62 * math.cos(mid), cy + r * 0.62 * math.sin(mid)
            labels.append(
                f'<text x="{num(tx)}" y="{num(ty + 5)}" text-anchor="middle" font-family=\'{SANS}\' '
                f'font-size="14" font-weight="600" fill="{ink_on(fill)}">{n}</text>'
            )
    out.extend(labels)
    # The legend sits to the right of the pie, vertically centred on it.
    line, gap = 23.5, 9
    legend_x = pad_x + 200 + 40
    legend_top = cy - (len(LAWS) * line + (len(LAWS) - 1) * gap) / 2
    for i, (law, name, n, p) in enumerate(zip(LAWS, names, counts, pcts, strict=True)):
        mid = legend_top + i * (line + gap) + line / 2
        out.append(
            f'<rect x="{legend_x}" y="{num(mid - 6)}" width="12" height="12" rx="3" fill="{c[f"law-{law}"]}"/>'
        )
        out.append(
            f'<text x="{legend_x + 21}" y="{num(mid + 5)}" font-family=\'{SANS}\' font-size="14.5" '
            f'fill="{c["ink-2"]}">{esc(name)}<tspan dx="9" font-weight="600" fill="{c["ink"]}">{n}</tspan>'
            f'<tspan dx="9" fill="{c["muted"]}">{p}%</tspan></text>'
        )
    out.append("</svg>\n")
    return "\n".join(out)


def track_path(x: float, y: float, w: float, h: float, r: float) -> str:
    """A bar square on the left and rounded on the right, as the site's track is."""
    r = min(r, w / 2, h / 2)
    return (
        f"M{num(x)} {num(y)}H{num(x + w - r)}A{num(r)} {num(r)} 0 0 1 {num(x + w)} {num(y + r)}"
        f"V{num(y + h - r)}A{num(r)} {num(r)} 0 0 1 {num(x + w - r)} {num(y + h)}H{num(x)}Z"
    )


def bars_svg(html: str, data: dict, theme: str) -> str:
    c = palette(html, theme)
    models = chart_order(html, data)
    names = {law: english(html, f"law.{law}") for law in LAWS}
    probes_label = english(html, "chart.probes")
    title_text = english(html, "chart.bars.title")
    passed_key, not_passed, total_label = (
        english(html, "chart.passedkey"),
        english(html, "chart.notpassed"),
        english(html, "chart.total"),
    )
    width, pad_x, top, gap = 860, 20, 22, 18
    legend_top = top + 24 + 14
    grid_top = legend_top + 22 + 14
    head_h, row_h = 46, 22
    height = grid_top + head_h + row_h * len(models) + 14
    name_w = max(len(m) for m in models) * 12.3 * 0.62
    total_w = 52
    content_w = width - 2 * pad_x
    panel_w = (content_w - name_w - total_w - 4 * gap) / 3
    col_x = [pad_x + name_w + gap + i * (panel_w + gap) for i in range(3)]
    right = width - pad_x
    value_w = 3.2 * 11.5
    scale_w = panel_w - 6 - value_w

    out = [
        svg_open(width, height, f"{title_text}, for each evaluated model"),
        card(width, height, c),
    ]
    out.append(
        f'<text x="{pad_x}" y="{top + 17}" font-family=\'{CONDENSED}\' font-size="15" '
        f'font-weight="600" fill="{c["ink"]}">{esc(title_text)}</text>'
    )
    # Legend: the three law colours mean "passed", the track "not passed".
    mid = legend_top + 11
    for i, law in enumerate(LAWS):
        out.append(
            f'<rect x="{pad_x + i * 14}" y="{num(mid - 6)}" width="12" height="12" rx="3" fill="{c[f"law-{law}"]}"/>'
        )
    x = pad_x + 3 * 14 - 2 + 9
    out.append(
        f'<text x="{x}" y="{num(mid + 5)}" font-family=\'{SANS}\' font-size="13.5" fill="{c["ink-2"]}">{esc(passed_key)}</text>'
    )
    x += text_width(passed_key, 13.5) + 20
    out.append(
        f'<rect x="{num(x)}" y="{num(mid - 6)}" width="12" height="12" rx="3" fill="{c["track"]}"/>'
    )
    out.append(
        f'<text x="{num(x + 21)}" y="{num(mid + 5)}" font-family=\'{SANS}\' font-size="13.5" fill="{c["ink-2"]}">{esc(not_passed)}</text>'
    )
    # One panel per law, spanning the header and every row.
    grid_h = head_h + row_h * len(models)
    for i, law in enumerate(LAWS):
        out.append(
            f'<rect x="{num(col_x[i] - 7 + 0.5)}" y="{num(grid_top + 0.5)}" width="{num(panel_w + 14 - 1)}" '
            f'height="{num(grid_h - 1)}" rx="8" fill="none" stroke="{c["line-strong"]}"/>'
        )
        out.append(
            f'<rect x="{num(col_x[i])}" y="{num(grid_top + 10)}" width="12" height="12" rx="3" fill="{c[f"law-{law}"]}"/>'
        )
        out.append(
            f'<text x="{num(col_x[i] + 18)}" y="{num(grid_top + 20.5)}" font-family=\'{SANS}\' font-size="13" '
            f'font-weight="600" fill="{c["ink"]}">{esc(names[law])}</text>'
        )
        count = probes_label.replace("{n}", str(data["probes"][law]))
        out.append(
            f'<text x="{num(col_x[i])}" y="{num(grid_top + 35.5)}" font-family=\'{SANS}\' font-size="11.5" '
            f'fill="{c["muted"]}">{esc(count)}</text>'
        )
    out.append(
        f'<text x="{num(right)}" y="{num(grid_top + head_h - 9)}" text-anchor="end" font-family=\'{SANS}\' '
        f'font-size="11.5" fill="{c["muted"]}">{esc(total_label)}</text>'
    )
    for row, model in enumerate(models):
        counts = data["models"][model]
        mid = grid_top + head_h + row * row_h + row_h / 2
        out.append(
            f'<text x="{pad_x}" y="{num(mid + 4.3)}" font-family=\'{MONO}\' font-size="12.3" '
            f'fill="{c["ink"]}">{esc(model)}</text>'
        )
        for i, law in enumerate(LAWS):
            p, s = counts[law]
            all_ = data["probes"][law] or 1
            if s:
                track_w = max(2.0, s / all_ * scale_w)
                out.append(
                    f'<path d="{track_path(col_x[i], mid - 5, track_w, 10, 3)}" fill="{c["track"]}"/>'
                )
                fill_w = p / s * track_w
                if p == s:
                    out.append(
                        f'<path d="{track_path(col_x[i], mid - 5, track_w, 10, 3)}" fill="{c[f"law-{law}"]}"/>'
                    )
                elif p:
                    out.append(
                        f'<rect x="{num(col_x[i])}" y="{num(mid - 5)}" width="{num(fill_w)}" height="10" fill="{c[f"law-{law}"]}"/>'
                    )
            value, colour = (f"{p}/{s}", c["ink-2"]) if s else ("–", c["muted"])
            out.append(
                f'<text x="{num(col_x[i] + panel_w)}" y="{num(mid + 4)}" text-anchor="end" font-family=\'{SANS}\' '
                f'font-size="11.5" fill="{colour}">{value}</text>'
            )
        passed = sum(counts[law][0] for law in LAWS)
        scored = sum(counts[law][1] for law in LAWS)
        total = f"{passed} / {scored}" if scored else "–"
        out.append(
            f'<text x="{num(right)}" y="{num(mid + 4.3)}" text-anchor="end" font-family=\'{SANS}\' '
            f'font-size="12.5" fill="{c["ink-2"]}">{total}</text>'
        )
    out.append("</svg>\n")
    return "\n".join(out)


INFRA = re.compile(r'<svg class="infra" viewBox="0 0 (\d+) (\d+)"(.*?)</svg>', re.DOTALL)


def infra_svg(html: str, theme: str) -> str:
    """The site's "How it works" flowchart, standing alone: its English labels,
    the theme's colours in place of the CSS variables, and the figure's card
    behind it."""
    match = INFRA.search(html)
    if not match:
        raise SystemExit(f'{SITE} has no <svg class="infra"> flowchart')
    c = palette(html, theme)
    width, height, body = match.group(1), match.group(2), match.group(3)
    # The page replaces every labelled node's text with its English string on
    # load, so the string, not the markup's placeholder, is what the site shows.
    body = re.sub(
        r'(<(text|title)\b[^>]*?)\s+data-i18n="([^"]+)"([^>]*>)[^<]*(</\2>)',
        lambda m: m.group(1) + m.group(4) + esc(english(html, m.group(3))) + m.group(5),
        body,
    )
    if "data-i18n" in body:
        raise SystemExit("the flowchart has a translated node this script cannot fill in")

    def colour(m: re.Match) -> str:
        if m.group(1) not in c:
            raise SystemExit(f"the flowchart uses --{m.group(1)}, which the {theme} palette lacks")
        return c[m.group(1)]

    body = re.sub(r"var\(--([\w-]+)\)", colour, body)
    lines = [line.strip() for line in body.strip().splitlines()]
    # The first line still carries the attributes after viewBox, up to the tag's end.
    head = lines[0]
    rest = "\n".join(lines[1:])
    figure = (
        f'<rect x="0.5" y="0.5" width="{int(width) - 1}" height="{int(height) - 1}" rx="10" '
        f'fill="{c["surface"]}" stroke="{c["line"]}"/>'
    )
    # The card goes behind everything, after the title, style and defs.
    rest = rest.replace("</defs>", "</defs>\n" + figure, 1)
    return f'<svg class="infra" width="{width}" height="{height}" viewBox="0 0 {width} {height}" {head}\n{rest}\n</svg>\n'


def readme_images(html: str, data: dict) -> dict[Path, str]:
    paths = image_paths()
    images = {}
    for theme in THEMES:
        images[paths[f"pie-{theme}"]] = pie_svg(html, data, theme)
        images[paths[f"bars-{theme}"]] = bars_svg(html, data, theme)
        images[paths[f"infra-{theme}"]] = infra_svg(html, theme)
    return images


README = ROOT / "README.md"
NUMBERS = re.compile(r"(<!-- law-numbers:[^\n]*-->\n)(.*?)(<!-- /law-numbers -->)", re.DOTALL)


def numbers_block(html: str, data: dict) -> str:
    """The site's "Show the numbers" table under the bars, as Markdown, so the
    per-law results are text a screen reader can reach and not only an image."""

    def cell(passed: int, scored: int) -> str:
        return f"{passed} / {scored}" if scored else "–"

    header = (
        [english(html, "models.th1")]
        + [english(html, f"law.{law}") for law in LAWS]
        + [english(html, "chart.total")]
    )
    lines = [
        "<details>",
        f"<summary>{esc(english(html, 'chart.table'))}</summary>",
        "",
        "| " + " | ".join(header) + " |",
        "| --- |" + " ---: |" * (len(LAWS) + 1),
    ]
    for model in chart_order(html, data):
        counts = data["models"][model]
        passed = sum(counts[law][0] for law in LAWS)
        scored = sum(counts[law][1] for law in LAWS)
        cells = " | ".join(cell(*counts[law]) for law in LAWS)
        lines.append(f"| `{model}` | {cells} | {cell(passed, scored)} |")
    lines += ["", "</details>"]
    return "\n".join(lines) + "\n"


def readme_numbers(readme: str) -> str:
    match = NUMBERS.search(readme)
    if not match:
        raise SystemExit(f"{README} has no <!-- law-numbers: --> ... <!-- /law-numbers --> block")
    return match.group(2)


def main(argv: list[str]) -> int:
    html = SITE.read_text(encoding="utf-8")
    readme = README.read_text(encoding="utf-8")
    data = standings()
    if "--check" in argv:
        status = 0
        if embedded(html) != data:
            print("site/index.html standings are stale; run python3 scripts/site_standings.py")
            status = 1
        else:
            print("site/index.html standings match models/result.csv")
        stale = [
            path.relative_to(ROOT).as_posix()
            for path, text in readme_images(html, data).items()
            if not path.exists() or path.read_text(encoding="utf-8") != text
        ]
        if readme_numbers(readme) != numbers_block(html, data):
            stale.append("README.md numbers table")
        if stale:
            print(
                f"README charts are stale ({', '.join(stale)}); run python3 scripts/site_standings.py"
            )
            status = 1
        else:
            print("README images and numbers table match the site")
        return status
    embedded(html)  # fail loudly if the block is missing
    updated = BLOCK.sub(lambda m: m.group(1) + render(data) + m.group(3), html, count=1)
    SITE.write_text(updated, encoding="utf-8")
    images = readme_images(updated, data)
    for path, text in images.items():
        path.write_text(text, encoding="utf-8")
    readme_numbers(readme)  # fail loudly if the block is missing
    table = numbers_block(updated, data)
    README.write_text(
        NUMBERS.sub(lambda m: m.group(1) + table + m.group(3), readme, count=1), encoding="utf-8"
    )
    print(
        f"wrote standings for {len(data['models'])} models and {sum(data['probes'].values())} probes, "
        f"and {len(images)} README images"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
