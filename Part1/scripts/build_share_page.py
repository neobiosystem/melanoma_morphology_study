#!/usr/bin/env python3
"""Build the shareable Part 1 findings page.

Reads the analysis outputs directly (calibration JSON, route comparison,
confound check, exemplar manifest) so the published numbers cannot drift
away from the numbers the MATLAB pipeline actually produced. Images are
inlined as data URIs because the artifact host blocks external requests.

Usage:  python3 scripts/build_share_page.py
"""

import base64
import csv
import json
import pathlib
import subprocess

ROOT = pathlib.Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
ASSETS = RESULTS / "share" / "assets"
OUT = RESULTS / "share" / "part1-border-morphology.html"

# Source figures, downscaled on first build. Full-resolution dermoscopy runs to
# several MB per file, which would blow past the artifact size cap once base64
# inflates it by a third.
FIGURES = ["fig1_descriptor_distributions_otsu", "fig2_boxcount_fit_otsu",
           "fig3_roc_otsu", "fig5_agreement_by_grade",
           "fig6_perimeter_scale_stability_otsu"]


def sips(src, dst, longest_edge, extra=()):
    dst.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(["sips", "-Z", str(longest_edge), *extra, str(src), "--out", str(dst)],
                   check=True, capture_output=True)


def build_assets():
    """Regenerate the downscaled copies the page embeds."""
    for name in FIGURES:
        dst = ASSETS / f"{name}.png"
        if not dst.exists():
            sips(RESULTS / "figures" / f"{name}.png", dst, 1400)
    for g in GRADES:
        dst = ASSETS / f"{g}.jpg"
        if not dst.exists():
            src = next((ROOT / "classical_examples").glob(f"{g}_01_*_image.jpg"))
            sips(src, dst, 900, ("-s", "format", "jpeg", "-s", "formatOptions", "68"))

GRADES = ["benign", "in_situ", "invasive"]
GRADE_CN = {"benign": "良性痣", "in_situ": "原位黑色素瘤", "invasive": "浸润性黑色素瘤"}
GRADE_EN = {"benign": "Nevus", "in_situ": "Melanoma in situ", "invasive": "Invasive melanoma"}

# Direction a descriptor is expected to move as severity rises.
# +1 = malignant should measure higher than benign.
EXPECTED = {
    "Df": +1, "compactness": +1, "circularity": -1, "asymmetry": +1,
    "solidity": -1, "convexity": -1, "perimeterDecaySlope": -1,
}
DESC_CN = {
    "Df": ("分形维数", "D<sub>f</sub>"),
    "compactness": ("紧致度", "compactness"),
    "asymmetry": ("不对称性", "asymmetry"),
    "solidity": ("实心度", "solidity"),
    "convexity": ("凸度", "convexity"),
    "circularity": ("圆度", "circularity"),
    "perimeterDecaySlope": ("周长衰减斜率", "perimeter decay"),
}
ORDER = ["compactness", "solidity", "convexity", "asymmetry",
         "perimeterDecaySlope", "circularity", "Df"]

ALPHA = 0.05


def b64(path, mime):
    return f"data:{mime};base64,{base64.b64encode(path.read_bytes()).decode()}"


def img(name, mime="image/png"):
    return b64(ASSETS / name, mime)


def fmt(x, n=4):
    return f"{x:.{n}f}"


def pnum(p):
    """Format a p-value with a real minus sign and no leading zero in the exponent."""
    if p == 0 or p < 1e-15:
        return "&lt; 1e&#8722;15"
    if p < 0.001:
        mant, exp = f"{p:.1e}".split("e")
        return f"{mant}e&#8722;{int(exp[1:])}"
    return f"{p:.3f}"


build_assets()

cal = json.loads((RESULTS / "calibration_targets_otsu.json").read_text())
route = {r["descriptor"]: r for r in csv.DictReader((RESULTS / "route_comparison.csv").open())}

# Pairwise Dunn results. A direction claim is a pairwise claim, so it must be
# backed by the pairwise p, never by the three-group omnibus p.
PH = {(r["route"], r["descriptor"], r["groupA"], r["groupB"]): r
      for r in csv.DictReader((RESULTS / "posthoc_all.csv").open())}


def ph(route_name, key, a, b):
    r = PH[(route_name, key, a, b)]
    return float(r["pAdjusted"]), float(r["rankBiserial"])
conf = {r["descriptor"]: r for r in csv.DictReader((RESULTS / "route_confound_check.csv").open())}
exemplars = list(csv.DictReader((ROOT / "classical_examples" / "manifest.csv").open()))


# ---------------------------------------------------------------- interval bars
def interval_bars(key):
    """Three IQR bars on a shared axis: overlap is the point, so show it."""
    d = cal["descriptors"][key]
    q1s = [d[g]["iqr"][0] for g in GRADES]
    q3s = [d[g]["iqr"][1] for g in GRADES]
    lo, hi = min(q1s), max(q3s)
    pad = 0.18 * (hi - lo)
    lo, hi = lo - pad, hi + pad
    span = hi - lo

    def pct(x):
        return max(0.0, min(100.0, (x - lo) / span * 100))

    rows = []
    for gi, g in enumerate(GRADES):
        s = d[g]
        left, right = pct(s["iqr"][0]), pct(s["iqr"][1])
        rows.append(f"""
        <div class="bar-row">
          <span class="bar-label">{GRADE_CN[g]}</span>
          <span class="bar-track">
            <span class="bar-fill g{gi+1}" style="left:{left:.2f}%;width:{max(right-left,0.8):.2f}%"></span>
            <span class="bar-tick g{gi+1}" style="left:{pct(s['median']):.2f}%"></span>
          </span>
          <span class="bar-val">{fmt(s['median'], 4 if key != 'compactness' else 3)}</span>
        </div>""")
    cn, en = DESC_CN[key]
    return f"""
      <figure class="ivl">
        <figcaption><b>{cn}</b> <span class="mono dim">{en}</span></figcaption>
        {''.join(rows)}
        <p class="ivl-foot">条形 = 四分位距 IQR，竖线 = 中位数</p>
      </figure>"""


# ---------------------------------------------------------------- interval table
def interval_rows():
    out = []
    for key in ORDER:
        d = cal["descriptors"][key]
        cn, en = DESC_CN[key]
        cells = []
        for g in GRADES:
            s = d[g]
            cells.append(
                f'<td><span class="med">{fmt(s["median"])}</span>'
                f'<span class="iqr">{fmt(s["iqr"][0])} – {fmt(s["iqr"][1])}</span></td>')
        out.append(
            f'<tr><th scope="row">{cn}<span class="mono dim">{en}</span></th>'
            f'{"".join(cells)}<td class="auc">{d["auc"]:.3f}</td></tr>')
    return "\n".join(out)


# ------------------------------------------------------------- reversal table
def reversal_rows():
    out = []
    for key in ORDER:
        r = route[key]
        cn, en = DESC_CN[key]
        exp = EXPECTED[key]
        cells = []
        for pre in ("expert", "otsu"):
            b, i = float(r[f"{pre}_medBen"]), float(r[f"{pre}_medInv"])
            p, rb = ph(pre, key, "benign", "invasive")
            good = (1 if i > b else -1) == exp
            sig = p < ALPHA
            arrow = "&#8593;" if i > b else "&#8595;"
            cls = "bad" if (not good and sig) else ("warn" if not good else "ok")
            if not good:
                mark = ('<span class="flag">方向相反 &middot; 显著</span>' if sig
                        else '<span class="flag warn-flag">方向相反 &middot; 不显著</span>')
            else:
                mark = ('' if sig
                        else '<span class="flag warn-flag">不显著</span>')
            cells.append(
                f'<td class="{cls}"><span class="arrow">{arrow}</span>'
                f'<span class="pair">{b:.3f} &#8594; {i:.3f}</span>'
                f'<span class="pstat">p<sub>adj</sub> {pnum(p)} &middot; rb {rb:+.2f}</span>{mark}</td>')
        eta_e, eta_o = float(r["expert_eta2"]), float(r["otsu_eta2"])
        out.append(
            f'<tr><th scope="row">{cn}<span class="mono dim">{en}</span></th>'
            f'{"".join(cells)}'
            f'<td class="ratio">{eta_o/eta_e:.0f}&#215;</td></tr>')
    return "\n".join(out)


# ------------------------------------------------------------- two-tier table
def tier_rows():
    out = []
    for key in ORDER:
        cn, _ = DESC_CN[key]
        cells = []
        for a, b in (("benign", "in_situ"), ("benign", "invasive"), ("in_situ", "invasive")):
            p, _ = ph("otsu", key, a, b)
            if p < ALPHA:
                cells.append(f'<td class="sig">{pnum(p)}</td>')
            else:
                cells.append(f'<td class="ns">{p:.3f}<span class="flag ns-flag">n.s.</span></td>')
        out.append(f'<tr><th scope="row">{cn}</th>{"".join(cells)}</tr>')
    return "\n".join(out)


# ----------------------------------------------------------------- exemplars
def exemplar_cards():
    cards = []
    for g in GRADES:
        row = next(e for e in exemplars if e["gradeKey"] == g and e["rank"] == "1")
        cards.append(f"""
      <figure class="ex">
        <img src="{img(g + '.jpg', 'image/jpeg')}" alt="{GRADE_EN[g]} 典型示例 {row['isic_id']}">
        <figcaption>
          <span class="ex-grade">{GRADE_CN[g]}</span>
          <span class="mono dim">{row['isic_id']}</span>
          <dl>
            <dt>D<sub>f</sub></dt><dd>{float(row['Df']):.4f}</dd>
            <dt>紧致度</dt><dd>{float(row['compactness']):.3f}</dd>
            <dt>Dice</dt><dd>{float(row['dice']):.3f}</dd>
          </dl>
        </figcaption>
      </figure>""")
    return "\n".join(cards)


# ----------------------------------------------------------------- confounds
def confound_rows():
    out = []
    for key in ORDER:
        c = conf[key]
        cn, _ = DESC_CN[key]
        a, g = float(c["allEta2"]), float(c["goodEta2"])
        out.append(
            f'<tr><th scope="row">{cn}</th><td>{a:.4f}</td><td>{g:.4f}</td>'
            f'<td class="{"ok" if g > 0.6 * a else "bad"}">{g/a*100:.0f}%</td></tr>')
    return "\n".join(out)


HTML = f"""<title>黑色素瘤边界形态图谱</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Newsreader:ital,opsz,wght@0,6..72,400;0,6..72,600;1,6..72,400&amp;family=Archivo:wght@400;500;600&amp;family=IBM+Plex+Mono:wght@400;500&amp;family=Noto+Serif+SC:wght@400;600&amp;family=Noto+Sans+SC:wght@400;500&amp;display=swap">
<style>
:root {{
  --paper:#F1F3F5; --surface:#FCFCFD; --ink:#15181D; --muted:#5B636F;
  --rule:#D8DDE3; --rule-soft:#E7EAEF;
  --g1:#C0854A; --g2:#9E4E2B; --g3:#4E2318;
  --accent:#0E7A85; --accent-soft:#D7EAEC;
  --ok:#2C7A57; --bad:#AF3628; --bad-soft:#F6E2DF;
  --shadow:0 1px 2px rgba(20,24,30,.05), 0 8px 24px -16px rgba(20,24,30,.28);
  --serif:'Newsreader','Noto Serif SC',Georgia,'Songti SC',serif;
  --sans:'Archivo','Noto Sans SC','PingFang SC',system-ui,sans-serif;
  --mono:'IBM Plex Mono','SF Mono',Menlo,monospace;
}}
@media (prefers-color-scheme:dark) {{
  :root:not([data-theme="light"]) {{
    --paper:#0E1116; --surface:#161A21; --ink:#E7EAEE; --muted:#98A2AE;
    --rule:#272E38; --rule-soft:#1E242B;
    --g1:#D6A366; --g2:#C36F41; --g3:#8C5140;
    --accent:#46BCC7; --accent-soft:#12333A;
    --ok:#4FB286; --bad:#E07A69; --bad-soft:#3A211D;
    --shadow:0 1px 2px rgba(0,0,0,.4), 0 8px 24px -16px rgba(0,0,0,.8);
  }}
}}
:root[data-theme="dark"] {{
  --paper:#0E1116; --surface:#161A21; --ink:#E7EAEE; --muted:#98A2AE;
  --rule:#272E38; --rule-soft:#1E242B;
  --g1:#D6A366; --g2:#C36F41; --g3:#8C5140;
  --accent:#46BCC7; --accent-soft:#12333A;
  --ok:#4FB286; --bad:#E07A69; --bad-soft:#3A211D;
  --shadow:0 1px 2px rgba(0,0,0,.4), 0 8px 24px -16px rgba(0,0,0,.8);
}}

*, *::before, *::after {{ box-sizing:border-box; }}
body {{
  margin:0; background:var(--paper); color:var(--ink);
  font-family:var(--sans); font-size:16.5px; line-height:1.72;
  -webkit-font-smoothing:antialiased;
}}
.wrap {{ max-width:1080px; margin:0 auto; padding:0 28px 96px; }}
.col {{ max-width:62ch; }}
p {{ margin:0 0 1.05em; }}
a {{ color:var(--accent); text-underline-offset:3px; }}
:focus-visible {{ outline:2px solid var(--accent); outline-offset:3px; border-radius:3px; }}
.mono {{ font-family:var(--mono); font-variant-numeric:tabular-nums; }}
.dim {{ color:var(--muted); }}
sub {{ font-size:.72em; }}

/* ---- masthead ---- */
header.top {{ padding:76px 0 40px; border-bottom:1px solid var(--rule); }}
.eyebrow {{
  font-family:var(--mono); font-size:11.5px; letter-spacing:.16em;
  text-transform:uppercase; color:var(--accent); margin:0 0 20px;
}}
h1 {{
  font-family:var(--serif); font-weight:600; font-size:clamp(38px,6.2vw,62px);
  line-height:1.08; letter-spacing:-.015em; margin:0 0 22px; text-wrap:balance;
}}
.lede {{ font-family:var(--serif); font-size:20.5px; line-height:1.62; color:var(--muted); max-width:56ch; margin:0; }}

.stats {{ display:flex; flex-wrap:wrap; gap:0; margin:44px 0 0; border-top:1px solid var(--rule-soft); }}
.stat {{ flex:1 1 150px; padding:20px 22px 18px; border-right:1px solid var(--rule-soft); }}
.stat:last-child {{ border-right:0; }}
.stat b {{ display:block; font-family:var(--mono); font-size:27px; font-weight:500; letter-spacing:-.02em; }}
.stat span {{ display:block; font-size:12.5px; color:var(--muted); margin-top:3px; }}

/* ---- sections ---- */
section {{ padding:60px 0 0; }}
h2 {{
  font-family:var(--serif); font-weight:600; font-size:31px; line-height:1.22;
  letter-spacing:-.012em; margin:0 0 8px; text-wrap:balance;
}}
h2 .num {{ font-family:var(--mono); font-size:13px; color:var(--accent); margin-right:12px; font-weight:400; vertical-align:2px; }}
.sub {{ color:var(--muted); font-size:15px; margin:0 0 28px; max-width:60ch; }}
h3 {{ font-family:var(--sans); font-weight:600; font-size:16.5px; margin:34px 0 10px; }}

/* ---- verdict cards ---- */
.verdicts {{ display:grid; gap:14px; grid-template-columns:repeat(auto-fit,minmax(268px,1fr)); margin-top:8px; }}
.v {{
  background:var(--surface); border:1px solid var(--rule); border-radius:4px;
  padding:22px 22px 20px; box-shadow:var(--shadow);
}}
.v h3 {{ margin:0 0 8px; font-size:16px; }}
.v p {{ margin:0; font-size:14.5px; line-height:1.62; color:var(--muted); }}
.v .tag {{
  font-family:var(--mono); font-size:10.5px; letter-spacing:.13em; text-transform:uppercase;
  color:var(--accent); display:block; margin-bottom:11px;
}}

/* ---- tables ---- */
.scroll {{ overflow-x:auto; margin:26px 0 0; -webkit-overflow-scrolling:touch; }}
table {{ border-collapse:collapse; width:100%; min-width:600px; font-size:14.5px; }}
caption {{ text-align:left; font-size:13px; color:var(--muted); padding-bottom:12px; }}
th, td {{ text-align:left; padding:11px 14px; border-bottom:1px solid var(--rule-soft); vertical-align:top; }}
thead th {{
  font-family:var(--mono); font-size:11px; letter-spacing:.09em; text-transform:uppercase;
  color:var(--muted); font-weight:500; border-bottom:1px solid var(--rule);
  white-space:nowrap;
}}
tbody th {{ font-weight:500; white-space:nowrap; }}
tbody th .mono {{ display:block; font-size:11px; font-weight:400; }}
td {{ font-family:var(--mono); font-variant-numeric:tabular-nums; }}
.med {{ display:block; }}
.iqr {{ display:block; font-size:11.5px; color:var(--muted); }}
.auc {{ font-weight:500; }}
tbody tr:hover {{ background:var(--rule-soft); }}

/* ---- reversal ---- */
.arrow {{ font-size:17px; margin-right:8px; }}
td.ok .arrow {{ color:var(--ok); }}
td.bad .arrow, td.warn .arrow {{ color:var(--bad); }}
td.bad {{ background:var(--bad-soft); }}
.pair {{ font-size:13px; }}
.pstat {{ display:block; font-size:11px; color:var(--muted); margin-top:2px; }}
.flag {{
  display:block; font-family:var(--sans); font-size:10.5px; letter-spacing:.06em;
  color:var(--bad); margin-top:3px; font-weight:600;
}}
.warn-flag {{ color:var(--muted); }}
.ratio {{ color:var(--accent); font-weight:500; }}
td.sig {{ color:var(--ok); }}
td.ns {{ color:var(--bad); }}
.ns-flag {{ color:var(--muted); font-weight:500; }}

/* ---- interval bars ---- */
.ivls {{ display:grid; gap:22px; grid-template-columns:repeat(auto-fit,minmax(330px,1fr)); margin-top:28px; }}
.ivl {{ margin:0; background:var(--surface); border:1px solid var(--rule); border-radius:4px; padding:20px 22px 16px; }}
.ivl figcaption {{ font-size:15px; margin-bottom:16px; }}
.ivl figcaption .mono {{ font-size:11.5px; margin-left:7px; }}
.bar-row {{ display:flex; align-items:center; gap:12px; margin-bottom:9px; }}
.bar-label {{ flex:0 0 92px; font-size:12px; color:var(--muted); }}
.bar-track {{ position:relative; flex:1 1 auto; height:15px; background:var(--rule-soft); border-radius:2px; }}
.bar-fill {{ position:absolute; top:3px; height:9px; border-radius:2px; opacity:.55; }}
.bar-tick {{ position:absolute; top:0; width:2px; height:15px; border-radius:1px; }}
.bar-fill.g1, .bar-tick.g1 {{ background:var(--g1); }}
.bar-fill.g2, .bar-tick.g2 {{ background:var(--g2); }}
.bar-fill.g3, .bar-tick.g3 {{ background:var(--g3); }}
.bar-val {{ flex:0 0 58px; text-align:right; font-family:var(--mono); font-size:12.5px; font-variant-numeric:tabular-nums; }}
.ivl-foot {{ font-size:11.5px; color:var(--muted); margin:12px 0 0; }}

/* ---- exemplars ---- */
.exs {{ display:grid; gap:18px; grid-template-columns:repeat(auto-fit,minmax(252px,1fr)); margin-top:28px; }}
.ex {{ margin:0; background:var(--surface); border:1px solid var(--rule); border-radius:4px; overflow:hidden; box-shadow:var(--shadow); }}
.ex img {{ display:block; width:100%; height:212px; object-fit:cover; }}
.ex figcaption {{ padding:16px 18px 18px; }}
.ex-grade {{ display:block; font-weight:600; font-size:15px; }}
.ex figcaption .mono {{ font-size:11px; }}
.ex dl {{ display:grid; grid-template-columns:auto auto; gap:2px 14px; margin:12px 0 0; font-size:13px; }}
.ex dt {{ color:var(--muted); }}
.ex dd {{ margin:0; font-family:var(--mono); font-variant-numeric:tabular-nums; text-align:right; }}

/* ---- figures ---- */
figure.plot {{ margin:30px 0 0; }}
figure.plot img {{ width:100%; display:block; border:1px solid var(--rule); border-radius:4px; background:#fff; }}
figure.plot figcaption {{ font-size:13px; color:var(--muted); margin-top:11px; max-width:62ch; }}

/* ---- callout ---- */
.note {{
  border-left:2px solid var(--accent); background:var(--accent-soft);
  padding:18px 22px; margin:28px 0 0; border-radius:0 4px 4px 0; font-size:14.5px;
}}
.note p:last-child {{ margin-bottom:0; }}
.note b {{ font-weight:600; }}

footer {{ margin-top:70px; padding-top:26px; border-top:1px solid var(--rule); font-size:13px; color:var(--muted); }}
footer code {{ font-family:var(--mono); font-size:12px; }}

@media (max-width:640px) {{
  .wrap {{ padding:0 18px 64px; }}
  header.top {{ padding:48px 0 30px; }}
  .stat {{ flex-basis:50%; border-bottom:1px solid var(--rule-soft); }}
}}
@media (prefers-reduced-motion:reduce) {{ * {{ transition:none !important; animation:none !important; }} }}
</style>

<div class="wrap">

<header class="top">
  <p class="eyebrow">Part 1 &middot; 数据源阶段 &middot; ISIC 2018</p>
  <h1>黑色素瘤边界的形态学指标</h1>
  <p class="lede">从 2223 例皮肤镜图像中提取分形维数与紧致度，为第三阶段的细胞动力学仿真提供标定靶值。本页记录了这些数值，以及一个会让结论完全反号的陷阱。</p>
  <div class="stats">
    <div class="stat"><b>2223</b><span>入组皮损</span></div>
    <div class="stat"><b>1581</b><span>主分析可用</span></div>
    <div class="stat"><b>0.819</b><span>双路 Dice 中位数</span></div>
    <div class="stat"><b>34/34</b><span>验证测试通过</span></div>
  </div>
</header>

<section>
  <h2>三个核心结论</h2>
  <div class="verdicts">
    <div class="v">
      <span class="tag">测量陷阱</span>
      <h3>专家掩膜会把结论测反</h3>
      <p>ISIC 官方标注掩膜大量为多边形手绘。在其上测量，7 项指标中 <b>5 项方向与生物学预期相反，且 5 项全部统计显著</b>——浸润性黑色素瘤反而比良性痣更规整。改用真实像素分割后，7 项全部方向正确。</p>
    </div>
    <div class="v">
      <span class="tag">分层结构</span>
      <h3>形态学只分两级，不是三级</h3>
      <p>良性 vs 黑色素瘤的差异极显著，但<b>原位与浸润之间 7 项指标全部不显著</b>（p<sub>adj</sub> 0.057–0.40）。边界不规则性反映径向生长，而浸润是垂直深度事件。</p>
    </div>
    <div class="v">
      <span class="tag">指标权重</span>
      <h3>D<sub>f</sub> 是最弱的描述子</h3>
      <p>分形维数效应量 &eta;&sup2;=0.024，而紧致度达 0.087。真正承载信号的是周长与凸包类指标。原计划以 D<sub>f</sub> 为核心，数据不支持这一预设。</p>
    </div>
  </div>
</section>

<section>
  <h2><span class="num">01</span>形态学指标区间</h2>
  <p class="sub">主分析路线（自动分割，n = 1581）。上行为中位数，下行为四分位距。AUC 为良性 vs 黑色素瘤的判别力。</p>
  <div class="scroll">
    <table>
      <caption>各分级形态学指标：中位数与 IQR</caption>
      <thead><tr><th>描述子</th><th>良性痣<br>n = 1348</th><th>原位<br>n = 124</th><th>浸润<br>n = 109</th><th>AUC</th></tr></thead>
      <tbody>{interval_rows()}</tbody>
    </table>
  </div>
  <div class="ivls">
    {interval_bars("compactness")}
    {interval_bars("Df")}
  </div>
  <div class="note">
    <p><b>关于 D<sub>f</sub> 的一处更正。</b>早期试点得到的 1.10–1.14 中包含 Chan-Vese 主动轮廓的伪影——它沿局部强度抖动边界，会虚增盒维数。关闭后的诚实值为 <span class="mono">1.0400 / 1.0638 / 1.0479</span>。三组均接近 1，说明这些边界在盒计数尺度上接近可求长曲线，而非强分形。</p>
  </div>
</section>

<section>
  <h2><span class="num">02</span>两级而非三级</h2>
  <p class="sub">用户最初的设想是良性 / 轻度恶性 / 重度恶性三级区间。数据给出的是两级。</p>
  <div class="scroll">
    <table>
      <caption>Dunn 事后检验（Holm 校正），主分析路线</caption>
      <thead><tr><th>描述子</th><th>良性 vs 原位</th><th>良性 vs 浸润</th><th>原位 vs 浸润</th></tr></thead>
      <tbody>{tier_rows()}</tbody>
    </table>
  </div>
  <p style="margin-top:24px;max-width:62ch">这不是样本量不足的问题——原位组 n=124、浸润组 n=109，足以检出中等效应。两组中位数本身就几乎重合（紧致度 1.800 vs 1.810）。D<sub>f</sub> 甚至<b>非单调</b>：原位组 1.0638 高于浸润组 1.0479。</p>
  <div class="note">
    <p><b>对第三阶段的直接含义。</b>基于二维边界形态标定的细胞动力学仿真，可以区分痣与黑色素瘤，但<b>不应被期望编码浸润深度</b>。若 Part 3 试图用边界形状复现原位→浸润的转变，它将在标定阶段就失去可分辨的靶值。</p>
  </div>
</section>

<section>
  <h2><span class="num">03</span>专家掩膜的反号效应</h2>
  <p class="sub">两条分割路线走完全相同的测量代码，唯一差别是掩膜来源。箭头表示良性→浸润的变化方向，p<sub>adj</sub> 为该配对的 Dunn 检验（Holm 校正），rb 为 rank-biserial 效应量。</p>
  <div class="scroll">
    <table>
      <caption>良性 vs 浸润：方向、显著性与效应量</caption>
      <thead><tr><th>描述子</th><th>专家掩膜（多边形手绘）</th><th>自动分割（真实像素）</th><th>&eta;&sup2; 倍数</th></tr></thead>
      <tbody>{reversal_rows()}</tbody>
    </table>
  </div>
  <div class="note">
    <p><b>反号是显著的，但很小。</b>5 项反号的配对检验全部通过 Holm 校正（p<sub>adj</sub> 8.6e&#8722;5 至 1.4e&#8722;2），所以这不是噪声中的随机方向。但效应量确实很小（rank-biserial 0.13–0.21，&eta;&sup2; 0.003–0.009）——在 n = 2088 下，小效应很容易达到显著。准确的说法是<b>系统性、可复现</b>，而不是"差异很大"。反过来，方向正确的两项里只有不对称性显著，圆度并不显著（p<sub>adj</sub> = 0.22）。</p>
  </div>
  <p style="margin-top:24px;max-width:62ch">掩膜的共线边界顶点比中位数为 0.587，24.7% 的掩膜超过 0.85——多边形边在构造上局部笔直，其盒维数被钉在 1 附近。但<b>反号无法用"细节更少"解释</b>：细节缺失只会把效应压向零，不会把它推过零。最可能的机制是标注粗糙度随病灶尺寸放大：更大、更不规则的病灶得到成比例更粗的多边形，反而比小的良性痣更平滑，从而翻转依赖边界细节的指标。</p>
  <p style="max-width:62ch">这一机制有一个可检验的推论，数据也支持：<b>没有反号的两项恰好是不依赖边界细节的两项</b>。不对称性由主轴折叠后的掩膜重合面积计算，圆度由质心到边界的径向距离离散度计算——两者都是"质量分布在哪里"的度量，对局部边界抖动不敏感。而反号的 5 项全部依赖边界细节：周长长度（紧致度、凸度、周长衰减）、凹陷结构（实心度）、或精细边界本身（D<sub>f</sub>）。</p>
  <p style="max-width:62ch">因此主分析路线改为自动分割，专家掩膜降为敏感性臂。这与最初计划中"专家路线为主"的设定相反——是数据推翻了计划假设。</p>
  <figure class="plot">
    <img src="{img('fig3_roc_otsu.png')}" alt="七个描述子的 ROC 曲线，紧致度、实心度、凸度 AUC 约 0.74">
    <figcaption>良性 vs 黑色素瘤的 ROC。紧致度、实心度、凸度三者 AUC 均约 0.74，D<sub>f</sub> 仅 0.625。</figcaption>
  </figure>
</section>

<section>
  <h2><span class="num">04</span>典型示例</h2>
  <p class="sub">每级选取 3 例，此处展示各级第 1 例。典型性定义为：在全部 7 个描述子上同时接近本级中位数（稳健 z 的中位数最小），且掩膜质量达标。完整 9 例见 <span class="mono">classical_examples/</span>。</p>
  <div class="exs">{exemplar_cards()}</div>
  <div class="note">
    <p><b>选取规则中的一处偏倚修正。</b>候选被限制在分割质量达标（Dice &ge; 0.75）的病灶内，但参照中位数取自<b>全部</b>入组病灶。两者不可互换：Dice 门槛会优先保留更平滑的病灶，因为不规则边界恰是两条路线分歧最大的地方。若以过滤后子集为参照，浸润组示例的紧致度会落在 1.35 附近，而公布的中位数是 1.81。</p>
  </div>
</section>

<section>
  <h2><span class="num">05</span>两个混淆，都已检验</h2>
  <p class="sub">自动路线的优势可能来自两种非生物学解释。两者都做了检验，都没有推翻结论。</p>

  <h3>分割质量随恶性度下降</h3>
  <p style="max-width:62ch">Jonckheere-Terpstra 单调趋势检验 z = &minus;3.38，p = 0.0007——一致性确实随严重程度下降。但限定在分割良好的病灶（Dice &ge; 0.75，n = 1232）后，7 项指标全部保持显著，效应量仅温和衰减。若差异源于分割伪影，效应量应当坍塌。</p>
  <div class="scroll">
    <table>
      <caption>效应量 &eta;&sup2;：全体 vs 仅分割良好者</caption>
      <thead><tr><th>描述子</th><th>全体 n=1581</th><th>Dice &ge; 0.75 n=1232</th><th>保留比例</th></tr></thead>
      <tbody>{confound_rows()}</tbody>
    </table>
  </div>

  <h3>残留毛发</h3>
  <p style="max-width:62ch;margin-top:22px">形态学开运算因侵蚀病灶周缘而被关闭，毛发抑制只剩 DullRazor。毛发细、暗、附着于边界——正是会虚增盒维数的结构，且 D<sub>f</sub> 确实随毛发比例上升（Spearman &rho; = 0.168，p = 1.6e&#8722;11），毛发负荷本身也随分级不同。但限定在毛发最少的 40% 后，D<sub>f</sub> 的效应量<b>近乎翻倍</b>（&eta;&sup2; 0.024 &rarr; 0.046）。毛发是在掩盖信号，而非制造信号。</p>
</section>

<section>
  <h2><span class="num">06</span>栅格周长偏差</h2>
  <p class="sub">本研究最大的假阳性风险：栅格周长被系统性高估，且越不规则的边界高估越多。</p>
  <p style="max-width:62ch">若不加处理，黑色素瘤"更不紧致"可能纯粹是数字化伪影——这恰好就是本研究的核心论点以伪影形式出现。缓解方式是在固定物理尺度（窗宽 D<sub>eq</sub>/20）上平滑轮廓后再量周长，并跨多个尺度验证稳定性。</p>
  <figure class="plot">
    <img src="{img('fig6_perimeter_scale_stability_otsu.png')}" alt="平滑前后各分级紧致度中位数对比柱状图">
    <figcaption>良性→浸润的紧致度差：原始 +0.8716，固定尺度平滑后 +0.4583。差距未坍塌，说明约一半是几何真实的；另一半确为栅格伪影——这正是必须做这项检验的原因。</figcaption>
  </figure>
</section>

<section>
  <h2><span class="num">07</span>方法要点</h2>
  <p class="sub">三个设计选择承担了大部分效度负担。</p>
  <div class="verdicts">
    <div class="v">
      <span class="tag">盒计数</span>
      <h3>取下确界，不取均值</h3>
      <p>盒维数定义为覆盖的下确界，因此 N(&epsilon;) 在 3&times;3 网格原点上取<b>最小值</b>。取均值估计的是另一个量。</p>
    </div>
    <div class="v">
      <span class="tag">边界处理</span>
      <h3>D<sub>f</sub> 只在原始边界上计算</h3>
      <p>平滑后再盒计数会直接抹掉要测的分形信号。只有周长派生指标才使用平滑轮廓。</p>
    </div>
    <div class="v">
      <span class="tag">分辨率</span>
      <h3>&epsilon; 无量纲化，不重采样</h3>
      <p>把 4500 px 病灶降到 512 px 会优先摧毁高分辨率图像上的精细边界，等于把研究对象本身均质化。改用 &epsilon;/D<sub>eq</sub>。</p>
    </div>
  </div>
  <div class="note">
    <p><b>最关键的一条验证。</b>喂入实心圆盘时，正确实现应返回 D<sub>f</sub> &asymp; 2.0（它覆盖的是二维集合）。若实现悄悄测的是填充区域的边界，会返回 &asymp; 1.0，而后续所有数字都将失去意义。该测试已通过（<span class="mono">testFilledDiskIsDimensionTwo</span>）。</p>
  </div>
</section>

<section>
  <h2><span class="num">08</span>交付给第三阶段</h2>
  <p class="sub">机器可读的标定靶值，避免在多阶段研究之间手抄数字。</p>
  <div class="scroll">
    <table>
      <thead><tr><th>文件</th><th>内容</th></tr></thead>
      <tbody>
        <tr><th scope="row" class="mono">calibration_targets_otsu.json</th><td>主标定靶值：各级各描述子的中位数、IQR、2.5–97.5 参考区间</td></tr>
        <tr><th scope="row" class="mono">calibration_targets_expert.json</th><td>敏感性分析臂，<b>不可</b>用作主靶值</td></tr>
        <tr><th scope="row" class="mono">descriptors.csv</th><td>2223 例 &times; 7 描述子 &times; 2 路线的逐病灶原始测量</td></tr>
        <tr><th scope="row" class="mono">classical_examples/</th><td>9 例典型示例，含图像、掩膜、四联面板与许可信息</td></tr>
      </tbody>
    </table>
  </div>
  <p style="margin-top:26px;max-width:62ch">ISIC 图像为<b>逐图授权</b>（CC-0 / CC-BY / CC-BY-NC 混合）。本页展示的 9 例典型示例均为 CC-0 且全部经组织病理确诊，可无限制复用；其余图像的许可与署名随结果一并留存于 <span class="mono">manifest.csv</span>。</p>
</section>

<footer>
  <p>ISIC 2018 Task1-2 Training（collection 63），MATLAB R2023a + Image Processing Toolbox。分级由 ISIC Archive API 按 <code>diagnosis_3</code> 回连，旧字段 <code>benign_malignant</code> / <code>mel_class</code> 已于 2024 年废弃。</p>
  <p>全部数值由 <code>results/</code> 下的分析产物直接生成，未经手工转录。</p>
</footer>

</div>
"""

OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(HTML, encoding="utf-8")
print(f"wrote {OUT}  ({len(HTML.encode()) / 1024 / 1024:.2f} MB)")
