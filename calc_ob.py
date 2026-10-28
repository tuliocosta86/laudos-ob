#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CENTRAL DE CÁLCULOS — LAUDOS DE OBSTETRÍCIA (versão 2, auditada contra as fontes em 21/10/2026)
===========================================================================================
Uso (linha de comando):
    python3 calc_ob.py '{"data_exame":"30/09/2026","dum":"12/02/2026","dbp":78,"cc":290,"ca":265,"cf":58,"pfe":null,
                         "ip_au":1.05,"ip_acm":1.62,"pvs_acm":48,"ip_aut_d":0.95,"ip_aut_e":1.10,"ip_dv":0.55,
                         "mbv":5.2,"ila":13,"cl":33}'

    python3 calc_ob.py '{"data_exame":"30/09/2026","ccn":23.4,"dum":"10/07/2026","fc":158}'
    python3 calc_ob.py '{"data_exame":"30/09/2026","dmsg":18}'
    python3 calc_ob.py '{"peso_a":1850,"peso_b":1420}'
    python3 calc_ob.py '{"tn":2.9,"ccn":62}'
    python3 calc_ob.py '{"pbf":{"mrf":2,"mc":2,"tonus":2,"mbv":1.8,"ctg":2}}'

Todas as chaves são opcionais. Medidas em mm (DBP, CC, CA, CF, CCN, DMSG, TN, CL, VV), cm (MBV, ILA),
g (PFE, peso_a, peso_b), cm/s (PVS ACM). Datas dd/mm/aaaa. "ig_ref" (ex.: "32+2" ou "32s2d") = idade
gestacional de referência para percentis; se ausente, usa a IG pela DUM; se não houver DUM, usa a IG
ultrassonográfica calculada. Saída em português, pronta para o laudo. Referências ao final.

Chaves aceitas:
  data_exame, dum, ig_ref, ig_us_previa ("12+3"), data_us_previa
  dmsg (mm), vv (mm), ccn (mm), fc (bpm), tn (mm), bcf_ausente (true quando ditado "sem batimentos")
  ccn_a, ccn_b (mm, gemelar 1º trimestre)
  dbp, cc, ca, cf (mm), pfe (g, ditado)  → IG por Hadlock, PFE por Hadlock, percentil por Hadlock 1991;
      relação CC/CA com percentil (Campbell & Thoms 1977), relação CF/CA com percentil (Hadlock 1983)
      e altura fetal estimada pelo fêmur (Vintzileos 1984). Percentil das relações só com IG de referência
      do início do exame (ig_ref, US prévia ou DUM) — nunca pela biometria atual.
  ip_au, ip_acm, pvs_acm, ip_aut_d, ip_aut_e, ip_dv, incisura ("bilateral"/"unilateral"/"ausente")
  mbv (cm), ila (cm), cl (mm)
  peso_a, peso_b (g)
  pbf: {mrf, mc, tonus, mbv (cm) ou la (0/2), ctg (0/2, opcional)}
  metodo_ig: "robinson" (padrão) ou "hadlock" para CCN; "media" (padrão) ou "regressao" para biometria

Fontes conferidas (21/10/2026):
  CCN: Robinson & Fleming 1975 com o fator 1,037 (BMUS 2009); Hadlock 1992.
  IG pela biometria: Hadlock 1984 (coeficientes conferidos no original). Datação/redatação: CCN; depois CC
      (ISUOG 2019; INTERGROWTH-21st 2016), DBP se faltar a CC, média dos parâmetros só se faltarem as duas;
      limites da ACOG CO 700; não redata pela biometria se o Doppler estiver alterado pela DUM.
  PFE: Hadlock 1985 (4 fórmulas conferidas no original). Percentil: Hadlock 1991 — mediana pela equação do
      artigo, distribuição normal com DP de 13,28% da mediana, que reproduz exatamente a tabela publicada.
  IP AU: Acharya 2005. IP ACM e RCP: FMF, Ciobanu 2019. IP uterinas: Gómez 2008. PVS-ACM: Mari 2000
      (mediana = e^(2,31 + 0,04643 × IG)). IPV DV: Kessler 2006 (medianas conferidas; p5/p95 não conferidos).
"""

import json, math, sys
from datetime import date, datetime, timedelta

def br(x, casas=1):
    if x is None:
        return "—"
    if casas == 0:
        s = f"{int(round(x)):,}".replace(",", ".")
        return s
    s = f"{x:,.{casas}f}"
    return s.replace(",", "X").replace(".", ",").replace("X", ".")

def semanas_dias(dias):
    dias = int(round(dias))
    s, d = divmod(dias, 7)
    return f"{s} semana{'s' if s != 1 else ''} e {d} dia{'s' if d != 1 else ''}"

def sd_curto(dias):
    dias = int(round(dias))
    s, d = divmod(dias, 7)
    return f"{s}s{d}d"

def parse_data(s):
    if not s:
        return None
    return datetime.strptime(s.strip(), "%d/%m/%Y").date()

def parse_ig(s):
    if s is None:
        return None
    if isinstance(s, (int, float)):
        return int(round(float(s) * 7))
    t = str(s).lower().replace("semanas", "s").replace("semana", "s").replace("dias", "d").replace("dia", "d")
    t = t.replace("+", " ").replace("s", " ").replace("d", " ").replace("e", " ")
    partes = [p for p in t.split() if p.strip()]
    if not partes:
        return None
    sem = float(partes[0].replace(",", "."))
    dia = float(partes[1].replace(",", ".")) if len(partes) > 1 else 0
    return int(round(sem * 7 + dia))

def fmt_data(d):
    return d.strftime("%d/%m/%Y")

def phi(z):
    return 0.5 * (1 + math.erf(z / math.sqrt(2)))

def interp(tabela, x):
    ks = sorted(tabela)
    if x < ks[0] or x > ks[-1]:
        return None
    for i in range(len(ks) - 1):
        a, b = ks[i], ks[i + 1]
        if a <= x <= b:
            fa, fb = tabela[a], tabela[b]
            if b == a:
                return fa
            return fa + (fb - fa) * (x - a) / (b - a)
    return tabela[ks[-1]]

def percentil_de(valor, p5, p50, p95):
    """Log-normal com dispersão própria abaixo e acima da mediana: cai exatamente em p5, p50 e p95."""
    if valor is None or p50 is None or p5 is None or p95 is None or valor <= 0:
        return None
    sigma = math.log(p95 / p50) / 1.645 if valor >= p50 else math.log(p50 / p5) / 1.645
    z = (math.log(valor) - math.log(p50)) / sigma
    return max(0.1, min(99.9, 100 * phi(z)))

def fmt_pct(p, v=None, p5=None, p95=None):
    if p is None:
        return "—"
    if v is not None and p95 is not None and v > p95:
        return "> 95"
    if v is not None and p5 is not None and v < p5:
        return "< 5"
    if p < 1:
        return "< 1"
    if p > 99:
        return "> 99"
    return str(int(round(p)))

def ig_ccn_robinson(ccn_mm):
    return 8.052 * math.sqrt(ccn_mm * 1.037) + 23.73

def ig_ccn_hadlock(ccn_mm):
    c = ccn_mm / 10.0
    w = math.exp(1.684969 + 0.315646 * c - 0.049306 * c ** 2 + 0.004057 * c ** 3 - 0.000120456 * c ** 4)
    return w * 7

def ig_dmsg(dmsg_mm):
    return dmsg_mm + 30

def tn_p95(ccn_mm):
    """FMF (Snijders 1998): p95 de 2,1 mm com CCN 45 mm e 2,7 mm com 84 mm; entre eles, interpolação linear."""
    if ccn_mm < 45 or ccn_mm > 84:
        return None
    return 2.1 + (2.7 - 2.1) * (ccn_mm - 45) / (84 - 45)

def ig_hadlock_param(nome, mm):
    x = mm / 10.0
    if nome == "dbp":
        return 9.54 + 1.482 * x + 0.1676 * x ** 2
    if nome == "cc":
        return 8.96 + 0.540 * x + 0.0003 * x ** 3
    if nome == "ca":
        return 8.14 + 0.753 * x + 0.0036 * x ** 2
    if nome == "cf":
        return 10.35 + 2.460 * x + 0.170 * x ** 2
    return None

def ig_hadlock_regressao(dbp, cc, ca, cf):
    b, h, a, f = dbp / 10, cc / 10, ca / 10, cf / 10
    return 10.85 + 0.060 * h * f + 0.6700 * b + 0.1680 * a

def pfe_hadlock(dbp=None, cc=None, ca=None, cf=None):
    b = dbp / 10 if dbp else None
    h = cc / 10 if cc else None
    a = ca / 10 if ca else None
    f = cf / 10 if cf else None
    if a is None or f is None:
        return None, None
    if b and h:
        log10w = 1.3596 - 0.00386 * a * f + 0.0064 * h + 0.00061 * b * a + 0.0424 * a + 0.174 * f
        formula = "Hadlock (DBP, CC, CA, CF)"
    elif h:
        log10w = 1.326 - 0.00326 * a * f + 0.0107 * h + 0.0438 * a + 0.158 * f
        formula = "Hadlock (CC, CA, CF)"
    elif b:
        log10w = 1.335 - 0.0034 * a * f + 0.0316 * b + 0.0457 * a + 0.1623 * f
        formula = "Hadlock (DBP, CA, CF)"
    else:
        log10w = 1.304 + 0.05281 * a + 0.1938 * f - 0.004 * a * f
        formula = "Hadlock (CA, CF)"
    return 10 ** log10w, formula

# Tabela publicada de Hadlock 1991 (p3, p10, p50, p90, p97), mantida só para conferência.
# Semana 29: p97 = 1.724 (corrigido; a versão anterior tinha 1.754). Semana 30: o original imprime 1.649
# (erro de impressão, abaixo do p90); o valor coerente com a tabela é 1.949.
PFE_HADLOCK_1991 = {
    10: (26, 29, 35, 41, 44), 11: (34, 37, 45, 53, 56), 12: (43, 48, 58, 68, 73), 13: (55, 61, 73, 85, 91),
    14: (70, 77, 93, 109, 116), 15: (88, 97, 117, 137, 146), 16: (110, 121, 146, 171, 183),
    17: (136, 150, 181, 212, 226), 18: (167, 185, 223, 261, 279), 19: (205, 227, 273, 319, 341),
    20: (248, 275, 331, 387, 414), 21: (299, 331, 399, 467, 499), 22: (359, 398, 478, 559, 598),
    23: (426, 471, 568, 665, 710), 24: (503, 556, 670, 784, 838), 25: (589, 652, 785, 918, 981),
    26: (685, 758, 913, 1068, 1141), 27: (791, 876, 1055, 1234, 1319), 28: (908, 1004, 1210, 1416, 1513),
    29: (1034, 1145, 1379, 1613, 1724), 30: (1169, 1294, 1559, 1824, 1949), 31: (1313, 1453, 1751, 2049, 2189),
    32: (1465, 1621, 1953, 2285, 2441), 33: (1622, 1794, 2162, 2530, 2703), 34: (1783, 1973, 2377, 2781, 2971),
    35: (1946, 2154, 2595, 3036, 3244), 36: (2110, 2335, 2813, 3291, 3516), 37: (2271, 2513, 3028, 3543, 3785),
    38: (2427, 2686, 3236, 3786, 4045), 39: (2576, 2851, 3435, 4019, 4294), 40: (2714, 3004, 3619, 4234, 4524),
}
PFE_DP_REL = 0.1328

def pfe_mediana(semanas):
    """Hadlock 1991: ln(peso) = 0,578 + 0,332·IG − 0,00354·IG² (IG em semanas)."""
    return math.exp(0.578 + 0.332 * semanas - 0.00354 * semanas ** 2)

def percentil_pfe(peso_g, ig_dias):
    """Percentil exato pela curva de Hadlock 1991: normal com DP = 13,28% da mediana (reproduz a tabela publicada:
    p3 = 0,75; p10 = 0,83; p90 = 1,17; p97 = 1,25 × p50)."""
    w = ig_dias / 7.0
    if w < 10 or w > 40:
        return None, None
    p50 = pfe_mediana(w)
    z = (peso_g / p50 - 1) / PFE_DP_REL
    return max(0.1, min(99.9, 100 * phi(z))), (0.75 * p50, 0.83 * p50, p50, 1.17 * p50, 1.25 * p50)

# IP da artéria umbilical — Acharya G et al., Am J Obstet Gynecol 2005;192:937–944 (p5, p50, p95)
IP_AU = {
    19: (1.02, 1.30, 1.66), 20: (0.99, 1.27, 1.62), 22: (0.92, 1.19, 1.54), 24: (0.86, 1.12, 1.47),
    26: (0.80, 1.06, 1.41), 28: (0.75, 1.00, 1.35), 30: (0.70, 0.95, 1.29), 32: (0.66, 0.90, 1.25),
    34: (0.62, 0.86, 1.20), 36: (0.58, 0.82, 1.16), 38: (0.55, 0.78, 1.12), 40: (0.51, 0.75, 1.09),
}
# IP da artéria cerebral média — FMF, Ciobanu A et al., Ultrasound Obstet Gynecol 2019;53:465–472
IP_ACM = {
    20: (1.16, 1.49, 1.90), 21: (1.21, 1.54, 1.96), 22: (1.26, 1.59, 2.01), 23: (1.31, 1.65, 2.07),
    24: (1.36, 1.70, 2.14), 25: (1.40, 1.76, 2.20), 26: (1.44, 1.80, 2.25), 27: (1.48, 1.85, 2.31),
    28: (1.50, 1.88, 2.36), 29: (1.52, 1.91, 2.39), 30: (1.52, 1.92, 2.42), 31: (1.52, 1.92, 2.44),
    32: (1.50, 1.91, 2.44), 33: (1.47, 1.89, 2.42), 34: (1.42, 1.85, 2.39), 35: (1.37, 1.79, 2.34),
    36: (1.30, 1.72, 2.27), 37: (1.22, 1.63, 2.18), 38: (1.13, 1.53, 2.07), 39: (1.03, 1.42, 1.95),
    40: (0.93, 1.30, 1.82), 41: (0.83, 1.17, 1.67),
}
# Relação cerebroplacentária — FMF, Ciobanu 2019
RCP = {
    20: (0.87, 1.21, 1.69), 21: (0.93, 1.29, 1.78), 22: (1.00, 1.37, 1.88), 23: (1.06, 1.45, 1.98),
    24: (1.12, 1.53, 2.08), 25: (1.18, 1.60, 2.18), 26: (1.24, 1.68, 2.28), 27: (1.29, 1.75, 2.38),
    28: (1.34, 1.82, 2.47), 29: (1.38, 1.87, 2.56), 30: (1.41, 1.92, 2.63), 31: (1.43, 1.96, 2.70),
    32: (1.44, 1.99, 2.75), 33: (1.43, 2.00, 2.79), 34: (1.42, 2.00, 2.81), 35: (1.39, 1.98, 2.81),
    36: (1.35, 1.94, 2.79), 37: (1.30, 1.89, 2.75), 38: (1.24, 1.83, 2.69), 39: (1.17, 1.75, 2.61),
    40: (1.09, 1.65, 2.51), 41: (1.00, 1.55, 2.39),
}
# IP médio das artérias uterinas — Gómez O et al., Ultrasound Obstet Gynecol 2008;32:128–132 (p5, p50, p95).
# p5 publicado de 11 a 31 semanas; de 32 a 41, estimado pela razão p5/p50 da semana 31 (só afeta o percentil
# impresso abaixo da mediana; a classificação usa o p95).
IP_AUT = {
    11: (1.18, 1.79, 2.70), 12: (1.11, 1.68, 2.53), 13: (1.05, 1.58, 2.38), 14: (0.99, 1.49, 2.24), 15: (0.94, 1.41, 2.11),
    16: (0.89, 1.33, 1.99), 17: (0.85, 1.27, 1.88), 18: (0.81, 1.20, 1.79), 19: (0.78, 1.15, 1.70), 20: (0.74, 1.10, 1.61),
    21: (0.71, 1.05, 1.54), 22: (0.69, 1.00, 1.47), 23: (0.66, 0.96, 1.41), 24: (0.64, 0.93, 1.35), 25: (0.62, 0.89, 1.30),
    26: (0.60, 0.86, 1.25), 27: (0.58, 0.84, 1.21), 28: (0.56, 0.81, 1.17), 29: (0.55, 0.79, 1.13), 30: (0.54, 0.77, 1.10),
    31: (0.52, 0.75, 1.06), 32: (0.51, 0.73, 1.04), 33: (0.49, 0.71, 1.01), 34: (0.49, 0.70, 0.99), 35: (0.48, 0.69, 0.97),
    36: (0.47, 0.68, 0.95), 37: (0.46, 0.67, 0.94), 38: (0.46, 0.66, 0.92), 39: (0.45, 0.65, 0.91), 40: (0.45, 0.65, 0.90),
    41: (0.45, 0.65, 0.89),
}
# IPV do ducto venoso — Kessler 2006 (medianas conferidas: 0,57 em 21 e 0,44 em 40 semanas; p5/p95 sem fonte acessível)
IP_DV = {
    20: (0.32, 0.56, 0.86), 24: (0.30, 0.54, 0.83), 28: (0.28, 0.52, 0.80), 32: (0.26, 0.50, 0.77),
    36: (0.24, 0.47, 0.73), 40: (0.22, 0.44, 0.70), 41: (0.22, 0.43, 0.69),
}

def pvs_acm_mediana(ig_dias):
    """Mari 2000: mediana da PVS-ACM = e^(2,31 + 0,04643 × IG em semanas)."""
    return math.exp(2.31 + 0.04643 * ig_dias / 7.0)

def ref3(tabela, ig_dias):
    w = ig_dias / 7.0
    ks = sorted(tabela)
    if w < ks[0] or w > ks[-1]:
        return None
    return tuple(interp(dict((k, tabela[k][i]) for k in ks), w) for i in range(len(tabela[ks[0]])))

# ---------------------------------------------------------------------------
# ÍNDICES BIOMÉTRICOS E ALTURA FETAL
# ---------------------------------------------------------------------------
# Relação CC/CA — Campbell & Thoms 1977 (BJOG 84:165): média, p5, p95 por faixa de 2 semanas
# (conferida com a reprodução da GLOWM, tabela 7).
_HC_AC_FAIXAS = [
    ((13, 14), (1.23, 1.14, 1.31)), ((15, 16), (1.22, 1.05, 1.39)), ((17, 18), (1.18, 1.07, 1.29)),
    ((19, 20), (1.18, 1.09, 1.26)), ((21, 22), (1.15, 1.06, 1.25)), ((23, 24), (1.13, 1.05, 1.21)),
    ((25, 26), (1.13, 1.04, 1.22)), ((27, 28), (1.13, 1.05, 1.22)), ((29, 30), (1.10, 0.99, 1.21)),
    ((31, 32), (1.07, 0.96, 1.17)), ((33, 34), (1.04, 0.96, 1.11)), ((35, 36), (1.02, 0.93, 1.11)),
    ((37, 38), (0.98, 0.92, 1.05)), ((39, 40), (0.97, 0.87, 1.06)), ((41, 42), (0.96, 0.93, 1.00)),
]
HC_AC_CAMPBELL = {s: v for (a, b), v in _HC_AC_FAIXAS for s in (a, b)}

# Relação CF/CA — Hadlock 1983 (AJR 141:979): a partir de 21 semanas, 22 ± 2% (média ± 2 DP → DP = 1%);
# limite superior da normalidade 23,5% (p90).
CF_CA_MEDIA, CF_CA_DP, CF_CA_IG_MIN = 22.0, 1.0, 21 * 7

def pct_assimetrico(v, media, p5, p95):
    """Percentil supondo normal com DP distinto abaixo e acima da média, ajustado para cair exatamente em p5/p95."""
    dp = (p95 - media) / 1.645 if v >= media else (media - p5) / 1.645
    return max(0.1, min(99.9, 100 * phi((v - media) / dp)))

def pct_cc_ca(rel, ig_dias):
    ref = HC_AC_CAMPBELL.get(int(ig_dias // 7))
    if not ref:
        return None
    return pct_assimetrico(rel, *ref), ref

def pct_cf_ca(rel_pct, ig_dias):
    if ig_dias < CF_CA_IG_MIN:
        return None
    return max(0.1, min(99.9, 100 * phi((rel_pct - CF_CA_MEDIA) / CF_CA_DP)))

def altura_fetal_cm(cf_mm):
    """Vintzileos 1984 (Obstet Gynecol 64:779): comprimento fetal (cm) = 6,18 + 0,59 × CF (mm)."""
    return 6.18 + 0.59 * cf_mm

def limite_acog(ig_dias, por_ccn):
    """ACOG CO 700: diferença máxima entre a DUM e o ultrassom antes de redatar."""
    if por_ccn:
        return 5 if ig_dias <= 62 else 7
    w = ig_dias / 7
    return 7 if w < 16 else 10 if w < 22 else 14 if w < 28 else 21

def main(d):
    out = []
    refs = set()
    alertas = []
    hoje = parse_data(d.get("data_exame")) or date.today()
    metodo_ccn = d.get("metodo_ig", "robinson")
    metodo_bio = d.get("metodo_biometria", "media")

    ig_dum = None
    dum = parse_data(d.get("dum"))
    if dum:
        ig_dum = (hoje - dum).days
        out.append(f"IG pela DUM: {semanas_dias(ig_dum)} ({sd_curto(ig_dum)}). DPP pela DUM: {fmt_data(dum + timedelta(days=280))}.")

    ig_us_previa = None
    if d.get("ig_us_previa") and d.get("data_us_previa"):
        dprev = parse_data(d["data_us_previa"])
        ig_us_previa = parse_ig(d["ig_us_previa"]) + (hoje - dprev).days
        out.append(f"IG por US prévia (datação anterior projetada para hoje): {semanas_dias(ig_us_previa)} ({sd_curto(ig_us_previa)}).")

    ig_us = None
    datacao = None  # (dias, método) usado para comparar com a DUM

    if d.get("dmsg"):
        ig = ig_dmsg(float(d["dmsg"]))
        out.append(f"DMSG {br(float(d['dmsg']),1)} mm → IG aproximada {semanas_dias(ig)} ({sd_curto(ig)}) [DMSG + 30 dias; usar só sem embrião].")
        if ig_us is None:
            ig_us = ig
        refs.add("Saco gestacional: IG ≈ DMSG (mm) + 30 dias")
        if float(d["dmsg"]) >= 25 and not d.get("ccn"):
            alertas.append("DMSG ≥ 25 mm sem embrião → critério definitivo de inviabilidade (SRU 2013).")
        elif 16 <= float(d["dmsg"]) < 25 and not d.get("ccn"):
            alertas.append("DMSG 16–24 mm sem embrião → suspeito, mas não definitivo de inviabilidade (SRU 2013); controle em 7–14 dias.")

    if d.get("vv") and float(d["vv"]) > 7:
        alertas.append(f"Vesícula vitelínica de {br(float(d['vv']),1)} mm (> 7 mm) → achado suspeito de mau prognóstico (SRU 2013).")

    if d.get("ccn"):
        ccn = float(d["ccn"])
        ig = ig_ccn_hadlock(ccn) if metodo_ccn == "hadlock" else ig_ccn_robinson(ccn)
        ig_us = ig
        datacao = (ig, "ccn")
        fonte = "Hadlock 1992" if metodo_ccn == "hadlock" else "Robinson & Fleming 1975"
        out.append(f"CCN {br(ccn,1)} mm → IG {semanas_dias(ig)} ({sd_curto(ig)}) [{fonte}]. DPP pelo CCN: {fmt_data(hoje + timedelta(days=280 - round(ig)))}.")
        refs.add(f"CCN: {fonte}")
        if d.get("fc") is None and d.get("bcf_ausente"):
            if ccn >= 7:
                alertas.append("CCN ≥ 7 mm sem batimentos cardíacos → critério definitivo de inviabilidade (morte embrionária) (SRU 2013).")
            else:
                alertas.append("CCN < 7 mm sem batimentos → suspeito, não definitivo; controle em 7–14 dias (SRU 2013).")
        if d.get("fc") is not None:
            fc = float(d["fc"])
            if ig < 6 * 7 + 2 and fc < 100:
                alertas.append(f"FC {int(fc)} bpm com IG < 6s2d → bradicardia embrionária (< 100 bpm).")
            elif 6 * 7 + 2 <= ig <= 7 * 7 and fc < 120:
                alertas.append(f"FC {int(fc)} bpm entre 6s2d e 7s0d → bradicardia embrionária (< 120 bpm).")
            elif ig > 7 * 7 and fc < 110:
                alertas.append(f"FC {int(fc)} bpm após 7 semanas → bradicardia embrionária (< 110 bpm).")
        if d.get("dmsg"):
            dif = float(d["dmsg"]) - ccn
            if dif < 5:
                alertas.append(f"DMSG − CCN = {br(dif,1)} mm (< 5 mm) → saco gestacional pequeno (achado suspeito, SRU 2013).")

    if d.get("tn") is not None:
        tn = float(d["tn"])
        p95 = tn_p95(float(d["ccn"])) if d.get("ccn") else None
        linha = f"TN {br(tn,1)} mm"
        if p95:
            linha += f" (p95 estimado para o CCN: {br(p95,1)} mm)"
        if tn >= 3.5:
            linha += " — TN AUMENTADA (≥ 3,5 mm)"
        elif p95 and tn > p95:
            linha += " — acima do p95 para o CCN (relatar como aumentada)"
        elif p95:
            linha += " — abaixo do p95 para o CCN"
        out.append(linha + ".")
        refs.add("TN: corte ≥ 3,5 mm e p95 aproximado pelo CCN")

    bio = {k: float(d[k]) for k in ("dbp", "cc", "ca", "cf") if d.get(k) is not None}
    rel_ccca = rel_cfca = None
    if bio:
        igs = {k: ig_hadlock_param(k, v) for k, v in bio.items()}
        nomes = {"dbp": "DBP", "cc": "CC", "ca": "CA", "cf": "CF"}
        partes = [f"{nomes[k]} {br(bio[k],1)} mm → {sd_curto(igs[k]*7)}" for k in ("dbp", "cc", "ca", "cf") if k in igs]
        out.append("Biometria (Hadlock 1984): " + "; ".join(partes) + ".")
        if len(bio) == 4 and metodo_bio == "regressao":
            ig = ig_hadlock_regressao(bio["dbp"], bio["cc"], bio["ca"], bio["cf"]) * 7
            out.append(f"IG pela biometria (regressão Hadlock 4 parâmetros): {semanas_dias(ig)} ({sd_curto(ig)}).")
        else:
            ig = sum(igs.values()) / len(igs) * 7
            out.append(f"IG pela biometria (média dos parâmetros): {semanas_dias(ig)} ({sd_curto(ig)}).")
        ig_us = ig
        if datacao is None:
            if "cc" in igs:
                datacao = (igs["cc"] * 7, "cc")
            elif "dbp" in igs:
                datacao = (igs["dbp"] * 7, "dbp")
            else:
                datacao = (ig, "media")
            nome_dat = {"cc": "circunferência cefálica", "dbp": "diâmetro biparietal", "media": "média dos parâmetros"}[datacao[1]]
            out.append(f"IG de datação (ISUOG 2019: CC; DBP se faltar a CC; média só sem os dois): {semanas_dias(datacao[0])} "
                       f"({sd_curto(datacao[0])}), pela {nome_dat}. DPP pelo ultrassom: {fmt_data(hoje + timedelta(days=280 - round(datacao[0])))}.")
        refs.add("Biometria e IG: Hadlock 1984")
        pfe, formula = pfe_hadlock(bio.get("dbp"), bio.get("cc"), bio.get("ca"), bio.get("cf"))
        if d.get("pfe"):
            pfe = float(d["pfe"]); formula = "ditado"
        if pfe:
            pfe_r = int(round(pfe / 5.0) * 5)
            out.append(f"PFE: {br(pfe_r,0)} g [{formula}].")
            refs.add("Peso fetal: Hadlock 1985")
        if "cf" in bio:
            out.append(f"Altura fetal estimada: {br(altura_fetal_cm(bio['cf']),1)} cm [Vintzileos 1984: 6,18 + 0,59 × CF (mm)].")
            refs.add("Altura fetal: Vintzileos 1984")
        # relações arredondadas como aparecem no laudo (o percentil usa o valor arredondado)
        if "cc" in bio and "ca" in bio:
            rel_ccca = round(bio["cc"] / bio["ca"], 2)
        if "cf" in bio and "ca" in bio:
            rel_cfca = round(100 * bio["cf"] / bio["ca"], 1)
    elif datacao is None and ig_us is not None:
        datacao = (ig_us, "dmsg")

    ig_ref = parse_ig(d.get("ig_ref")) if d.get("ig_ref") else None
    origem_ref = "ditada"
    if ig_ref is None and ig_us_previa is not None:
        ig_ref, origem_ref = ig_us_previa, "US prévia"
    if ig_ref is None and ig_dum is not None:
        ig_ref, origem_ref = ig_dum, "DUM"
    if ig_ref is None and datacao is not None:
        ig_ref, origem_ref = datacao[0], "ultrassom atual (sem DUM/US prévia)"
    if ig_ref is not None:
        out.append(f"IG de referência para percentis: {semanas_dias(ig_ref)} ({sd_curto(ig_ref)}) — origem: {origem_ref}.")

    # Índices biométricos — percentil só com IG de referência do início do exame (nunca pela biometria atual)
    ref_inicio = ig_ref is not None and not origem_ref.startswith("ultrassom atual")
    sem_ref = " → sem IG de referência no início do exame: percentil não calculado (escrever só o valor)"
    if rel_ccca is not None:
        linha = f"Relação CC/CA: {br(rel_ccca,2)}"
        if ref_inicio:
            r = pct_cc_ca(rel_ccca, ig_ref)
            if r:
                p, (m, p5, p95) = r
                linha += (f" → percentil {fmt_pct(p)} para {sd_curto(ig_ref)} "
                          f"[Campbell & Thoms: média {br(m,2)}; p5 {br(p5,2)}; p95 {br(p95,2)}]")
                refs.add("Relação CC/CA: Campbell & Thoms 1977")
            else:
                linha += " → percentil não calculado (tabela de 13 a 42 semanas)"
        else:
            linha += sem_ref
        out.append(linha + ".")
    if rel_cfca is not None:
        linha = f"Relação CF/CA: {br(rel_cfca,1)}%"
        if ref_inicio:
            p = pct_cf_ca(rel_cfca, ig_ref)
            if p is not None:
                linha += (f" → percentil {fmt_pct(p)} [Hadlock 1983: {br(CF_CA_MEDIA,0)} ± {br(2*CF_CA_DP,0)}% (média ± 2 DP), "
                          f"a partir de 21 semanas; limite superior 23,5% (p90)]")
                refs.add("Relação CF/CA: Hadlock 1983")
            else:
                linha += " → percentil não calculado (referência válida a partir de 21 semanas)"
        else:
            linha += sem_ref
        out.append(linha + ".")

    if bio and ig_ref is not None:
        pfe_val = float(d["pfe"]) if d.get("pfe") else pfe_hadlock(bio.get("dbp"), bio.get("cc"), bio.get("ca"), bio.get("cf"))[0]
        if pfe_val:
            p, curva = percentil_pfe(pfe_val, ig_ref)
            if p is not None:
                p3, p10, p50, p90, p97 = curva
                classe = ("abaixo do percentil 3" if p < 3 else "abaixo do percentil 10" if p < 10 else
                          "acima do percentil 97" if p > 97 else "acima do percentil 90" if p > 90 else "entre os percentis 10 e 90")
                out.append(f"Percentil do PFE para {sd_curto(ig_ref)}: percentil {fmt_pct(p)} ({classe}) [Hadlock 1991: p10 {br(p10,0)} g; p50 {br(p50,0)} g; p90 {br(p90,0)} g].")
                refs.add("Percentil de peso: Hadlock 1991")
            else:
                out.append("Percentil do PFE: fora da faixa da curva (10–40 semanas).")

    doppler_alterado = False
    if ig_ref is not None:
        if d.get("ip_au") is not None:
            v = float(d["ip_au"]); r = ref3(IP_AU, ig_ref)
            if r:
                p = percentil_de(v, *r)
                doppler_alterado |= v > r[2]
                out.append(f"Artéria umbilical: IP {br(v,2)} → percentil {fmt_pct(p, v, r[0], r[2])} (p95 = {br(r[2],2)}){' — AUMENTO DA RESISTÊNCIA (> p95)' if v > r[2] else ''}.")
                refs.add("IP AU: Acharya 2005")
            else:
                out.append("Artéria umbilical: IG fora da faixa da referência (19–40 semanas): percentil não calculado.")
        if d.get("ip_acm") is not None:
            v = float(d["ip_acm"]); r = ref3(IP_ACM, ig_ref)
            if r:
                p = percentil_de(v, *r)
                doppler_alterado |= v < r[0]
                out.append(f"Artéria cerebral média: IP {br(v,2)} → percentil {fmt_pct(p, v, r[0], r[2])} (p5 = {br(r[0],2)}){' — REDUÇÃO DA RESISTÊNCIA (< p5)' if v < r[0] else ''}.")
                refs.add("IP ACM: FMF, Ciobanu 2019")
        if d.get("ip_au") is not None and d.get("ip_acm") is not None and float(d["ip_au"]) > 0:
            rcp = float(d["ip_acm"]) / float(d["ip_au"]); r = ref3(RCP, ig_ref)
            linha = f"Relação cerebroplacentária (IP ACM / IP AU): {br(rcp,2)}"
            if r:
                p = percentil_de(rcp, *r)
                linha += f" → percentil {fmt_pct(p, rcp, r[0], r[2])} (p5 = {br(r[0],2)})"
                if rcp < r[0]:
                    linha += " — ABAIXO DO p5 (centralização hemodinâmica)"
                    doppler_alterado = True
            if rcp < 1.0:
                linha += " — RCP < 1,0"
                doppler_alterado = True
            out.append(linha + ".")
            refs.add("RCP: FMF, Ciobanu 2019")
        if d.get("pvs_acm") is not None:
            v = float(d["pvs_acm"]); med = pvs_acm_mediana(ig_ref); mom = v / med
            out.append(f"PVS-ACM: {br(v,1)} cm/s → {br(mom,2)} MoM (mediana {br(med,1)} cm/s; 1,5 MoM = {br(1.5*med,1)} cm/s)"
                       + (" — ELEVADA (> 1,5 MoM: sugere anemia fetal)" if mom > 1.5 else "") + ".")
            refs.add("PVS-ACM: Mari 2000")
        if d.get("ip_aut_d") is not None or d.get("ip_aut_e") is not None:
            vals = [float(d[k]) for k in ("ip_aut_d", "ip_aut_e") if d.get(k) is not None]
            med = sum(vals) / len(vals); r = ref3(IP_AUT, ig_ref)
            linha = f"Artérias uterinas: IP médio {br(med,2)}"
            if r:
                p5, p50, p95 = r
                p = percentil_de(med, p5, p50, p95)
                linha += f" → percentil {fmt_pct(p, med, None, p95)} (p95 = {br(p95,2)})"
                if med > p95:
                    linha += " — AUMENTO DA RESISTÊNCIA (> p95)"
                    doppler_alterado = True
            inc = d.get("incisura")
            if inc and inc != "ausente":
                linha += f"; incisura protodiastólica {inc}"
            out.append(linha + ".")
            refs.add("IP artérias uterinas: Gómez 2008")
        if d.get("ip_dv") is not None:
            v = float(d["ip_dv"]); r = ref3(IP_DV, ig_ref)
            if r:
                p = percentil_de(v, *r)
                doppler_alterado |= v > r[2]
                out.append(f"Ducto venoso: IPV {br(v,2)} → percentil {fmt_pct(p, v, r[0], r[2])} (p95 = {br(r[2],2)}){' — AUMENTADO (> p95)' if v > r[2] else ''}.")
                refs.add("IPV ducto venoso: Kessler 2006")
    elif any(d.get(k) is not None for k in ("ip_au", "ip_acm", "pvs_acm", "ip_aut_d", "ip_aut_e", "ip_dv")):
        out.append("Doppler informado sem IG de referência (dite ig_ref, DUM ou biometria): percentis não calculados.")

    # Redatação (ACOG CO 700) — só quando a referência é a DUM (sem US prévio nem IG ditada)
    if origem_ref == "DUM" and datacao is not None and datacao[1] != "dmsg":
        dias_us, met = datacao
        lim = limite_acog(dias_us, met == "ccn")
        dif = abs(ig_dum - dias_us)
        nome = {"ccn": "CCN", "cc": "CC", "dbp": "DBP", "media": "biometria (média)"}[met]
        linha = f"Discrepância DUM × {nome}: {int(round(dif))} dia(s); limite para redatação (ACOG): {lim} dias → "
        if dif <= lim:
            linha += "manter DUM."
        elif met != "ccn" and doppler_alterado:
            linha += ("manter DUM: Doppler alterado pela DUM (possível restrição de crescimento) — não redatar pela "
                      "biometria; percentis pela DUM.")
        else:
            linha += "REDATAR pelo ultrassom" + (" (só se não houver datação prévia por CCN)." if met != "ccn" else ".")
            if met != "ccn":
                linha += " Percentis e Doppler devem ser recalculados com ig_ref = IG de datação."
        out.append(linha)
        refs.add("Redatação: ACOG CO 700")

    if d.get("mbv") is not None:
        v = float(d["mbv"])
        c = "oligoidrâmnio" if v < 2 else "polidrâmnio" if v >= 8 else "normal"
        out.append(f"Maior bolsão vertical: {br(v,1)} cm → {c}.")
    if d.get("ila") is not None:
        v = float(d["ila"])
        c = "oligoidrâmnio" if v <= 5 else "polidrâmnio" if v >= 24 else "normal"
        if 5 < v <= 8:
            c += " (limítrofe inferior, 5–8 cm)"
        out.append(f"ILA: {br(v,1)} cm → {c}.")

    if d.get("cl") is not None:
        v = float(d["cl"])
        c = "normal (≥ 25 mm)" if v >= 25 else "COLO CURTO (< 25 mm)" if v >= 15 else "COLO MUITO CURTO (< 15 mm)"
        out.append(f"Comprimento cervical: {br(v,1)} mm → {c}.")

    if d.get("peso_a") is not None and d.get("peso_b") is not None:
        a, b = float(d["peso_a"]), float(d["peso_b"])
        maior, menor = max(a, b), min(a, b)
        disc = (maior - menor) / maior * 100
        c = "≥ 25% — discordância significativa (ISUOG)" if disc >= 25 else "≥ 20% — discordância a relatar" if disc >= 20 else "< 20%"
        out.append(f"Discordância de peso entre os fetos: {br(disc,1)}% ({c}).")
    if d.get("ccn_a") is not None and d.get("ccn_b") is not None:
        a, b = float(d["ccn_a"]), float(d["ccn_b"])
        disc = abs(a - b) / max(a, b) * 100
        out.append(f"Discordância de CCN: {br(disc,1)}%" + (" (≥ 10% — relatar)" if disc >= 10 else "") + ".")

    if d.get("pbf"):
        p = d["pbf"]
        pts = {}
        for k in ("mrf", "mc", "tonus"):
            pts[k] = 2 if p.get(k, 2) in (2, True, "presente") else 0
        if p.get("mbv") is not None:
            pts["la"] = 2 if float(p["mbv"]) >= 2 else 0
        else:
            pts["la"] = 2 if p.get("la", 2) in (2, True, "normal") else 0
        total = sum(pts.values()); maximo = 8
        if p.get("ctg") is not None:
            pts["ctg"] = 2 if p["ctg"] in (2, True, "reativa") else 0
            total += pts["ctg"]; maximo = 10
        if total in (8, 10) or (total == 8 and maximo == 10):
            classe = "normal" if pts["la"] == 2 else "duvidoso (líquido alterado)"
        elif total == 6:
            classe = "duvidoso"
        else:
            classe = "alterado"
        out.append(f"Perfil biofísico fetal: {total}/{maximo} — {classe}. (MRF {pts['mrf']}, MC {pts['mc']}, tônus {pts['tonus']}, LA {pts['la']}"
                   + (f", CTG {pts['ctg']}" if 'ctg' in pts else "") + ")")

    print("=== RESULTADOS ===")
    for l in out:
        print("- " + l)
    if alertas:
        print("=== ALERTAS ===")
        for a in alertas:
            print("! " + a)
    if refs:
        print("=== REFERÊNCIAS USADAS ===")
        print("; ".join(sorted(refs)))

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__); sys.exit(0)
    main(json.loads(sys.argv[1]))
