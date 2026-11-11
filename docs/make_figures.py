"""Draw the README graphics from the model (run with `make figures`).

docs/rs25-schematic.svg  thrust chamber and bell nozzle to scale, with the givens
docs/composition.svg     equilibrium chamber composition of both engines

Plain SVG text, no plotting library. Every number comes from the model, so the
figures stay in step with the code.
"""
import contextlib
import io
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..'))
from exchanger import (FEE, OCC, ChamberTemperature, FEE_F, FCC,  # noqa: E402
                       FluorineChamberTemperature, OxNozzlePerf, MOLAR_MASS)

# ------------------------------------------------------------------ givens
P_C = 204*101325
FO = 0.166
T_H2, T_OX = 850, 530
R_T, R_E = 5.15, 45.35          # in
THETA_P = 32.0                  # deg, initial bell angle
THETA_E = 8.0                   # deg, exit angle used only for the drawing
PHI = (FO*MOLAR_MASS['O2']/MOLAR_MASS['H2'])/2


def run_model():
    with contextlib.redirect_stdout(io.StringIO()):
        T_fee, _, _, _, h_r = FEE(T_H2, P_C, T_OX, PHI)
        T_c = ChamberTemperature(T_fee, P_C, FO, h_r)
        g, R, _, _, df = OCC(T_c, P_C, FO)
        isp = OxNozzlePerf(R_T, R_E, g, T_c, R, P_C)[5]
        T_fee_f, _, h_rf = FEE_F(T_H2, P_C, T_OX, PHI)
        T_cf = FluorineChamberTemperature(T_fee_f, P_C, PHI, h_rf)
        gf, Rf, _, _, dff = FCC(T_cf, P_C, PHI)
        isp_f = OxNozzlePerf(R_T, R_E, gf, T_cf, Rf, P_C)[5]
    ox = dict(T_c=T_c, isp=isp, X=dict(zip(['H₂O', 'H₂', 'H', 'OH'], (float(x) for x in df['X_(i)'][0:4]))))
    fl = dict(T_c=T_cf, isp=isp_f, X=dict(zip(['HF', 'H₂', 'H', 'F'], (float(x) for x in dff['X_(i)'][0:4]))))
    return ox, fl


# ------------------------------------------------------------------ schematic
def schematic(ox, fl):
    W, H, cy, s = 960, 440, 220, 3.6          # px per inch
    x_inj, x_cyl, x_t = 120, 170, 212         # injector face, end of cylinder, throat
    r_c, r_t, r_e = 8.85*s, R_T*s, R_E*s      # chamber radius ~ 8.85 in (17.7 in diameter)
    L_bell = 0.8*(R_E - R_T)/math.tan(math.radians(15))*s   # 80% bell
    x_e = x_t + L_bell
    # bell: quadratic Bezier with tangents of theta_p at the throat and theta_e at the exit
    t1, t2 = math.tan(math.radians(THETA_P)), math.tan(math.radians(THETA_E))
    xq = (r_e - r_t + t1*x_t - t2*x_e)/(t1 - t2)
    rq = r_t + t1*(xq - x_t)

    def wall(sign):
        y = lambda r: cy - sign*r
        return (f'M{x_inj},{y(r_c):.1f} L{x_cyl},{y(r_c):.1f} '
                f'C{x_cyl+22},{y(r_c):.1f} {x_t-20},{y(r_t):.1f} {x_t},{y(r_t):.1f} '
                f'Q{xq:.1f},{y(rq):.1f} {x_e:.1f},{y(r_e):.1f}')

    # interior: top wall forward, bottom wall backward
    top = wall(1)
    bot_pts = wall(-1)
    interior = (f'{top} L{x_e:.1f},{cy + r_e:.1f} '
                f'Q{xq:.1f},{cy + rq:.1f} {x_t},{cy + r_t:.1f} '
                f'C{x_t-20},{cy + r_t:.1f} {x_cyl+22},{cy + r_c:.1f} {x_cyl},{cy + r_c:.1f} '
                f'L{x_inj},{cy + r_c:.1f} Z')
    plume = (f'M{x_e:.1f},{cy - r_e*0.96:.1f} C{x_e+120:.1f},{cy - r_e*0.9:.1f} {x_e+230:.1f},{cy - 70} {W-30},{cy - 46} '
             f'L{W-30},{cy + 46} C{x_e+230:.1f},{cy + 70} {x_e+120:.1f},{cy + r_e*0.9:.1f} {x_e:.1f},{cy + r_e*0.96:.1f} Z')

    def station(x, label, sub):
        return (f'<line class="station" x1="{x:.1f}" y1="34" x2="{x:.1f}" y2="{H-40}"/>'
                f'<circle class="stbadge" cx="{x:.1f}" cy="26" r="11"/>'
                f'<text class="sttxt" x="{x:.1f}" y="30">{label}</text>'
                f'<text class="muted small" x="{x:.1f}" y="{H-24}" text-anchor="middle">{sub}</text>')

    ang_x, ang_y = x_t + 34, cy - r_t - 34*t1
    return f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" font-family="-apple-system, 'Segoe UI', Helvetica, Arial, sans-serif">
  <title>RS-25 thrust chamber and nozzle</title>
  <desc>To-scale half-section of the RS-25 thrust chamber and bell nozzle used by the model: hydrogen at {T_H2} K and oxidizer at {T_OX} K enter at the injector, the chamber runs at 204 atm, the throat radius is {R_T} in and the exit radius {R_E} in (area ratio {(R_E/R_T)**2:.1f}).</desc>
  <style>
    .bg      {{ fill: #ffffff; }}
    .ink     {{ fill: #1f2933; }}
    .muted   {{ fill: #52606d; }}
    .small   {{ font-size: 12px; }}
    .wall    {{ fill: none; stroke: #3e4c59; stroke-width: 5; stroke-linejoin: round; stroke-linecap: round; }}
    .inj     {{ fill: #cbd2d9; stroke: #3e4c59; stroke-width: 1.5; }}
    .gas     {{ fill: url(#hot); }}
    .plume   {{ fill: url(#plume); }}
    .axis    {{ stroke: #9aa5b1; stroke-width: 1; stroke-dasharray: 10 4 2 4; }}
    .station {{ stroke: #9aa5b1; stroke-width: 1; stroke-dasharray: 3 4; }}
    .stbadge {{ fill: #1f2933; }}
    .sttxt   {{ fill: #ffffff; font-size: 12px; font-weight: 700; text-anchor: middle; }}
    .h2      {{ fill: none; stroke: #2680c2; stroke-width: 2.4; }}
    .ox      {{ fill: none; stroke: #c65d21; stroke-width: 2.4; }}
    .h2txt   {{ fill: #1f5f94; font-size: 13px; font-weight: 600; }}
    .oxtxt   {{ fill: #9a4210; font-size: 13px; font-weight: 600; }}
    .dim     {{ stroke: #52606d; stroke-width: 1; fill: none; }}
    .lbl     {{ fill: #1f2933; font-size: 13px; }}
    .lblb    {{ fill: #1f2933; font-size: 13px; font-weight: 600; }}
  </style>
  <defs>
    <linearGradient id="hot" x1="0" x2="1" y1="0" y2="0">
      <stop offset="0" stop-color="#fbe3c4"/><stop offset="0.12" stop-color="#f6a46b"/>
      <stop offset="0.2" stop-color="#f08a4b"/><stop offset="0.45" stop-color="#f9c79a"/>
      <stop offset="1" stop-color="#fdf1e4"/>
    </linearGradient>
    <linearGradient id="plume" x1="0" x2="1" y1="0" y2="0">
      <stop offset="0" stop-color="#fde6cf" stop-opacity="0.9"/><stop offset="1" stop-color="#fde6cf" stop-opacity="0"/>
    </linearGradient>
    <marker id="ah2" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto"><path d="M0,0 L10,5 L0,10 z" fill="#2680c2"/></marker>
    <marker id="aox" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto"><path d="M0,0 L10,5 L0,10 z" fill="#c65d21"/></marker>
  </defs>
  <rect class="bg" width="{W}" height="{H}"/>
  <path class="plume" d="{plume}"/>
  <path class="gas" d="{interior}"/>
  <line class="axis" x1="60" y1="{cy}" x2="{W-20}" y2="{cy}"/>
  <path class="wall" d="{top}"/>
  <path class="wall" d="{bot_pts}"/>
  <rect class="inj" x="{x_inj-12}" y="{cy - r_c - 6:.1f}" width="12" height="{2*r_c + 12:.1f}" rx="2"/>
  {station((x_inj + x_cyl)/2, 'c', 'chamber')}
  {station(x_t, 't', 'throat')}
  {station(x_e, 'e', 'exit')}
  <path class="h2" d="M24,{cy-52} L60,{cy-52} Q86,{cy-52} {x_inj-16},{cy-14}" marker-end="url(#ah2)"/>
  <text class="h2txt" x="20" y="{cy-62}">H₂ · {T_H2} K</text>
  <path class="ox" d="M24,{cy+52} L60,{cy+52} Q86,{cy+52} {x_inj-16},{cy+14}" marker-end="url(#aox)"/>
  <text class="oxtxt" x="20" y="{cy+74}">O₂ or F₂ · {T_OX} K</text>
  <text class="lblb" x="{x_inj-6}" y="{cy - r_c - 58:.1f}">Chamber · 204 atm</text>
  <text class="lbl" x="{x_inj-6}" y="{cy - r_c - 40:.1f}">T<tspan baseline-shift="sub" font-size="10">c</tspan> {ox['T_c']:.0f} K (O₂) · {fl['T_c']:.0f} K (F₂)</text>
  <path class="dim" d="M{ang_x:.1f},{cy - r_t:.1f} A34,34 0 0,0 {x_t + 34*math.cos(math.radians(THETA_P)):.1f},{cy - r_t - 34*math.sin(math.radians(THETA_P)):.1f}"/>
  <text class="lbl" x="{ang_x + 8:.1f}" y="{cy - r_t - 4:.1f}">θ<tspan baseline-shift="sub" font-size="10">p</tspan> = {THETA_P:.0f}°</text>
  <line class="dim" x1="{x_t}" y1="{cy + r_t + 4:.1f}" x2="{x_t - 26}" y2="{cy + r_c + 40:.1f}"/>
  <text class="lbl" x="{x_t - 30}" y="{cy + r_c + 54:.1f}" text-anchor="middle">R<tspan baseline-shift="sub" font-size="10">t</tspan> = {R_T} in</text>
  <line class="dim" x1="{x_e + 14:.1f}" y1="{cy}" x2="{x_e + 14:.1f}" y2="{cy - r_e:.1f}"/>
  <text class="lbl" x="{x_e + 20:.1f}" y="{cy - r_e/2:.1f}">R<tspan baseline-shift="sub" font-size="10">e</tspan> = {R_E} in</text>
  <text class="lbl" x="{x_e + 20:.1f}" y="{cy - r_e/2 + 18:.1f}">ε = {(R_E/R_T)**2:.1f}</text>
  <text class="lblb" x="{W-40}" y="{cy + 92}" text-anchor="end">Vacuum I<tspan baseline-shift="sub" font-size="10">sp</tspan></text>
  <text class="lbl" x="{W-40}" y="{cy + 110}" text-anchor="end">{ox['isp']:.1f} s (O₂) · {fl['isp']:.1f} s (F₂)</text>
</svg>
'''


# ------------------------------------------------------------------ composition chart
def composition(ox, fl):
    W, H = 960, 330
    panel_w, gap, left, top = 380, 90, 80, 96
    bar_h, bar_gap = 26, 16
    x_max = 0.8
    scale = panel_w/x_max

    def panel(x0, engine, data, cls):
        out = [f'<text class="ptitle" x="{x0}" y="{top-30}">{engine}</text>',
               f'<text class="psub" x="{x0}" y="{top-12}">T<tspan baseline-shift="sub" font-size="10">c</tspan> = {data["T_c"]:.0f} K</text>']
        for v in (0, 0.2, 0.4, 0.6, 0.8):
            gx = x0 + v*scale
            out.append(f'<line class="grid" x1="{gx:.1f}" y1="{top}" x2="{gx:.1f}" y2="{top + 4*(bar_h+bar_gap) - bar_gap}"/>')
            out.append(f'<text class="tick" x="{gx:.1f}" y="{top + 4*(bar_h+bar_gap) + 8}" text-anchor="middle">{v:.1f}</text>')
        for i, (sp, x) in enumerate(data['X'].items()):
            y = top + i*(bar_h + bar_gap)
            w = max(x*scale, 2)
            r = min(4, w/2)
            out.append(f'<text class="cat" x="{x0 - 10}" y="{y + bar_h/2 + 5}" text-anchor="end">{sp}</text>')
            out.append(f'<path class="{cls}" d="M{x0},{y} H{x0 + w - r:.1f} Q{x0 + w:.1f},{y} {x0 + w:.1f},{y + r} '
                       f'V{y + bar_h - r} Q{x0 + w:.1f},{y + bar_h} {x0 + w - r:.1f},{y + bar_h} H{x0} Z">'
                       f'<title>{engine}: {sp} mole fraction {x:.3f}</title></path>')
            out.append(f'<text class="val" x="{x0 + w + 8:.1f}" y="{y + bar_h/2 + 5}">{x:.3f}</text>')
        return '\n  '.join(out)

    p1 = panel(left, 'H₂/O₂', ox, 's1')
    p2 = panel(left + panel_w + gap, 'H₂/F₂', fl, 's2')
    return f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" font-family="-apple-system, 'Segoe UI', Helvetica, Arial, sans-serif">
  <title>Equilibrium chamber composition</title>
  <desc>Mole fractions at the equilibrium chamber temperature. H2/O2 at {ox["T_c"]:.0f} K: {", ".join(f"{k} {v:.3f}" for k, v in ox["X"].items())}. H2/F2 at {fl["T_c"]:.0f} K: {", ".join(f"{k} {v:.3f}" for k, v in fl["X"].items())}.</desc>
  <style>
    .surface {{ fill: #fcfcfb; }}
    .title   {{ fill: #0b0b0b; font-size: 16px; font-weight: 600; }}
    .ptitle  {{ fill: #0b0b0b; font-size: 14px; font-weight: 600; }}
    .psub, .cat {{ fill: #52514e; font-size: 13px; }}
    .tick    {{ fill: #52514e; font-size: 11px; }}
    .val     {{ fill: #0b0b0b; font-size: 12px; }}
    .grid    {{ stroke: #e6e5e0; stroke-width: 1; }}
    .s1      {{ fill: #2a78d6; }}
    .s2      {{ fill: #eb6834; }}
    @media (prefers-color-scheme: dark) {{
      .surface {{ fill: #1a1a19; }}
      .title, .ptitle, .val {{ fill: #f0efec; }}
      .psub, .cat, .tick {{ fill: #c3c2b7; }}
      .grid {{ stroke: #3a3a37; }}
      .s1 {{ fill: #3987e5; }}
      .s2 {{ fill: #d95926; }}
    }}
  </style>
  <rect class="surface" width="{W}" height="{H}" rx="8"/>
  <text class="title" x="{left - 50}" y="32">Equilibrium chamber composition (mole fraction)</text>
  {p1}
  {p2}
  <text class="tick" x="{W/2}" y="{H - 14}" text-anchor="middle">Mole fraction X</text>
</svg>
'''


if __name__ == '__main__':
    ox, fl = run_model()
    for name, svg in (('rs25-schematic.svg', schematic(ox, fl)), ('composition.svg', composition(ox, fl))):
        with open(os.path.join(HERE, name), 'w') as f:
            f.write(svg)
        print('wrote docs/' + name)
