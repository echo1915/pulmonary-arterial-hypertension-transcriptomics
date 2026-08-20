#!/usr/bin/env python3
"""Generate manuscript-ready PAH bulk/scRNA figures from preserved result tables."""

from __future__ import annotations

import math
import os
from pathlib import Path

import numpy as np
import pandas as pd
from reportlab.graphics import renderPDF, renderSVG
from reportlab.graphics.shapes import Circle, Drawing, Line, Path as RPath, Polygon, Rect, String
from reportlab.lib import colors


PROJECT = "project-001-pulmonary-arterial-hypertension-transcriptomics"
WORKDATA = Path(os.environ.get("WORKDATA_ROOT", r"G:\workdata"))
PROJECT_DATA = Path(os.environ.get("PAH_DATA_ROOT", WORKDATA / "projects" / PROJECT))
RESULTS = PROJECT_DATA / "outputs" / "current-results" / "pah_audit"
OUT = PROJECT_DATA / "outputs" / "manuscript-figures"
OUT.mkdir(parents=True, exist_ok=True)

W, H = 1600, 1100
NAVY = colors.HexColor("#17324D")
BLUE = colors.HexColor("#2F6B9A")
TEAL = colors.HexColor("#2A9D8F")
ORANGE = colors.HexColor("#E6953B")
RED = colors.HexColor("#C94C4C")
PURPLE = colors.HexColor("#7A5AA6")
GREY = colors.HexColor("#75808A")
LIGHT = colors.HexColor("#EEF3F6")
GRID = colors.HexColor("#D8E0E5")
INK = colors.HexColor("#263238")
WHITE = colors.white


def add_text(d, x, y, text, size=24, color=INK, anchor="start", bold=False):
    d.add(String(x, y, str(text), fontName="Helvetica-Bold" if bold else "Helvetica",
                 fontSize=size, fillColor=color, textAnchor=anchor))


def title(d, text, subtitle=None):
    add_text(d, 55, H - 62, text, 30, NAVY, bold=True)
    if subtitle:
        add_text(d, 55, H - 96, subtitle, 16, GREY)


def panel_label(d, x, y, label):
    add_text(d, x, y, label, 25, NAVY, bold=True)


def panel_box(d, x, y, w, h):
    d.add(Rect(x, y, w, h, rx=10, ry=10, fillColor=WHITE, strokeColor=GRID, strokeWidth=1.2))


def scale(v, lo, hi, a, b):
    if hi == lo:
        return (a + b) / 2
    return a + (float(v) - lo) / (hi - lo) * (b - a)


def color_interp(value, lo=-1, hi=1):
    value = max(lo, min(hi, float(value)))
    if value < 0:
        t = (value - lo) / (0 - lo)
        return colors.Color(0.20 + 0.80*t, 0.42 + 0.58*t, 0.72 + 0.28*t)
    t = value / hi if hi else 0
    return colors.Color(1.0, 1.0 - 0.70*t, 1.0 - 0.70*t)


def axes(d, x, y, w, h, xmin, xmax, ymin, ymax, xticks, yticks, xlabel="", ylabel=""):
    d.add(Line(x, y, x+w, y, strokeColor=INK, strokeWidth=1.2))
    d.add(Line(x, y, x, y+h, strokeColor=INK, strokeWidth=1.2))
    for v in xticks:
        px = scale(v, xmin, xmax, x, x+w)
        d.add(Line(px, y, px, y-5, strokeColor=INK))
        add_text(d, px, y-24, f"{v:g}", 12, GREY, "middle")
    for v in yticks:
        py = scale(v, ymin, ymax, y, y+h)
        d.add(Line(x-5, py, x, py, strokeColor=INK))
        d.add(Line(x, py, x+w, py, strokeColor=GRID, strokeWidth=0.5))
        add_text(d, x-10, py-4, f"{v:g}", 12, GREY, "end")
    if xlabel:
        add_text(d, x+w/2, y-48, xlabel, 14, INK, "middle")
    if ylabel:
        add_text(d, x, y+h+14, ylabel, 13, INK)


def save(d, stem):
    renderSVG.drawToFile(d, str(OUT / f"{stem}.svg"))
    renderPDF.drawToFile(d, str(OUT / f"{stem}.pdf"))


def load():
    meta = pd.read_csv(RESULTS / "lung_mechanism_random_effects_meta.csv")
    pri = pd.read_csv(RESULTS / "lung_mechanism_candidate_prioritization.csv")
    pb = pd.read_csv(RESULTS / "scrna_gse210248_candidate_pseudobulk.csv")
    ce = pd.read_csv(RESULTS / "scrna_gse210248_candidate_celltype_effects.csv")
    rep = pd.read_csv(RESULTS / "scrna_gse293580_candidate_effects.csv")
    spec = pd.read_csv(RESULTS / "scrna_pah_vs_phpf_candidate_specificity.csv")
    return meta, pri, pb, ce, rep, spec


def fig1(meta, pri):
    d = Drawing(W, H)
    d.add(Rect(0, 0, W, H, fillColor=colors.HexColor("#F8FAFB"), strokeColor=None))
    title(d, "Figure 1. Study design and stepwise identification of PDE8B",
          "Independent lung cohorts, disease-specific filtering, and single-cell replication")
    panel_label(d, 45, 945, "A")
    panel_box(d, 80, 575, 1440, 390)
    cohorts = [
        ("GSE15197", "18 PAH / 13 donor", BLUE),
        ("GSE113439", "14 PAH / 11 donor", TEAL),
        ("GSE254617", "82 Group 1 PAH / 52 donor", ORANGE),
        ("GSE208592", "15 PAH / 18 donor", PURPLE),
    ]
    for i, (name, n, c) in enumerate(cohorts):
        x = 125 + i*350
        d.add(Rect(x, 820, 260, 90, rx=12, ry=12, fillColor=colors.Color(c.red,c.green,c.blue,alpha=0.13), strokeColor=c, strokeWidth=2))
        add_text(d, x+130, 870, name, 21, c, "middle", True)
        add_text(d, x+130, 840, n, 14, INK, "middle")
        d.add(Line(x+130, 820, x+130, 770, strokeColor=GREY, strokeWidth=1.5))
    d.add(Rect(165, 690, 1270, 80, rx=12, ry=12, fillColor=LIGHT, strokeColor=NAVY, strokeWidth=2))
    add_text(d, 800, 735, "Random-effects meta-analysis across 15,165 common genes", 22, NAVY, "middle", True)
    add_text(d, 800, 705, "Cohort is the unit of replication; overlapping PHBI subjects excluded", 14, GREY, "middle")
    checks = [
        ("Disease specificity", "PAH vs IPF-PH (GSE15197)"),
        ("Secondary bulk", "GSE53408 filtered lung matrix"),
        ("Pulmonary-artery scRNA", "GSE210248: 3 PAH / 3 donor"),
        ("Independent whole-lung scRNA", "GSE293580: IPAH + SSc-PAH"),
    ]
    for i, (a,b) in enumerate(checks):
        x = 120 + i*355
        d.add(Rect(x, 600, 300, 65, rx=8, ry=8, fillColor=WHITE, strokeColor=GRID))
        add_text(d, x+150, 638, a, 15, NAVY, "middle", True)
        add_text(d, x+150, 615, b, 12, GREY, "middle")

    panel_label(d, 45, 515, "B")
    panel_box(d, 80, 80, 1440, 450)
    stages = [
        (15165, "Common genes", "4 cohorts"),
        (251, "Meta-significant", "FDR < 0.05"),
        (188, "Robust signature", "direction + I2 + |g|"),
        (77, "Priority genes", "specificity filters"),
        (3, "External bulk", "GSE53408"),
        (1, "Lead", "PDE8B"),
    ]
    widths = [220, 190, 170, 150, 130, 115]
    x = 110
    for i, ((n, lab, rule), ww) in enumerate(zip(stages, widths)):
        c = RED if i == len(stages)-1 else [NAVY,BLUE,TEAL,ORANGE,PURPLE][min(i,4)]
        yy = 335 - i*38
        d.add(Polygon([x,yy+95,x+ww,yy+95,x+ww-25,yy,x+25,yy], fillColor=colors.Color(c.red,c.green,c.blue,alpha=0.16), strokeColor=c, strokeWidth=2))
        add_text(d, x+ww/2, yy+60, f"{n:,}", 25, c, "middle", True)
        add_text(d, x+ww/2, yy+34, lab, 15, INK, "middle", True)
        add_text(d, x+ww/2, yy+14, rule, 10.5, GREY, "middle")
        x += ww - 10
    add_text(d, 1060, 455, "Why PDE8B advanced", 19, NAVY, bold=True)
    bullets = [
        "4/4 bulk cohorts: concordant increase",
        "Meta g = 0.744; I2 = 0%",
        "PAH vs IPF-PH: g = 0.901",
        "GSE53408 replication: g = 1.077",
        "Localized to SMC 1/2 in pulmonary artery",
        "Replicated in IPAH and SSc-PAH whole lung",
    ]
    for i,b in enumerate(bullets):
        d.add(Circle(1070, 415-i*48, 5, fillColor=RED, strokeColor=None))
        add_text(d, 1085, 409-i*48, b, 14, INK)
    save(d, "figure1_study_design_pde8b_selection")


def fig2(meta, pri):
    d = Drawing(W, H)
    d.add(Rect(0, 0, W, H, fillColor=colors.HexColor("#F8FAFB"), strokeColor=None))
    title(d, "Figure 2. Cross-cohort bulk lung transcriptomic evidence",
          "Random-effects meta-analysis and prioritization of robust PAH-associated genes")

    # A: volcano
    panel_label(d, 35, 960, "A")
    panel_box(d, 70, 565, 720, 410)
    x0,y0,w,h = 145,635,590,290
    yy = -np.log10(meta.fdr.clip(lower=1e-300))
    xmax = max(2.2, float(np.nanpercentile(np.abs(meta.meta_g), 99.5)))
    ymax = min(18, max(8, float(np.nanpercentile(yy, 99.7))))
    axes(d,x0,y0,w,h,-xmax,xmax,0,ymax,[-2,-1,0,1,2],[0,3,6,9,12,15],"Meta Hedges' g","-log10(FDR)")
    robust = meta.same_direction & (meta.fdr<0.05) & (meta.I2<50) & (meta.meta_g.abs()>=0.5)
    for idx,row in meta.iterrows():
        yv=min(ymax,float(yy.iloc[idx])); xv=max(-xmax,min(xmax,float(row.meta_g)))
        c = TEAL if robust.iloc[idx] else colors.HexColor("#B9C2C9")
        rr = 2.2 if robust.iloc[idx] else 1.1
        if row.symbol == "PDE8B": c,rr = RED,6
        d.add(Circle(scale(xv,-xmax,xmax,x0,x0+w),scale(yv,0,ymax,y0,y0+h),rr,fillColor=c,strokeColor=None,fillOpacity=0.65))
    pde=meta[meta.symbol.eq("PDE8B")].iloc[0]
    px=scale(pde.meta_g,-xmax,xmax,x0,x0+w); py=scale(min(ymax,-math.log10(max(pde.fdr,1e-300))),0,ymax,y0,y0+h)
    d.add(Line(px+5,py+5,px+65,py+38,strokeColor=RED,strokeWidth=1.5))
    add_text(d,px+70,py+34,"PDE8B",15,RED,bold=True)
    add_text(d,180,585,f"251 FDR-significant | 188 robust | PDE8B FDR={pde.fdr:.1e}",13,GREY)

    # B: heatmap
    panel_label(d, 805, 960, "B")
    panel_box(d, 840, 565, 700, 410)
    focus=["PDE8B","PIEZO2","SLC16A12","DCBLD1","NCAM1","SOX9","AQP5","CCL21","SELP","ABCG2","PDE3A","LPL"]
    hp=pri.set_index("symbol").reindex(focus).dropna(how="all")
    cols=["g_GSE15197","g_GSE113439","g_GSE254617","g_GSE208592","g_PAH_vs_IPF_PH_GSE15197","g_GSE53408_secondary"]
    labs=["15197","113439","254617","208592","PAH vs\nIPF-PH","53408"]
    hx,hy,cw,ch=1005,625,76,25
    for j,lab in enumerate(labs):
        for k,line in enumerate(lab.split("\n")):
            add_text(d,hx+j*cw+cw/2,930-k*15,line,11,NAVY,"middle",True)
    for i,(gene,row) in enumerate(hp.iterrows()):
        yy0=895-i*ch
        add_text(d,hx-15,yy0+5,gene,12,RED if gene=="PDE8B" else INK,"end",gene=="PDE8B")
        for j,cname in enumerate(cols):
            v=row[cname]
            fill=colors.HexColor("#ECEFF1") if pd.isna(v) else color_interp(v,-1.5,1.5)
            d.add(Rect(hx+j*cw,yy0-8,cw-3,ch-2,fillColor=fill,strokeColor=WHITE,strokeWidth=0.5))
            if not pd.isna(v): add_text(d,hx+j*cw+(cw-3)/2,yy0,f"{v:.2f}",9,INK,"middle")
    add_text(d,1175,590,"Blue: decreased   Red: increased",12,GREY,"middle")

    # C: PDE8B effects
    panel_label(d, 35, 505, "C")
    panel_box(d, 70, 75, 720, 440)
    prow=pri[pri.symbol.eq("PDE8B")].iloc[0]
    effects=[("GSE15197",prow.g_GSE15197),("GSE113439",prow.g_GSE113439),("GSE254617",prow.g_GSE254617),("GSE208592",prow.g_GSE208592),("PAH vs IPF-PH",prow.g_PAH_vs_IPF_PH_GSE15197),("GSE53408",prow.g_GSE53408_secondary)]
    fx,fy,fw,fh=300,140,430,310
    axes(d,fx,fy,fw,fh,-0.2,1.5,0,len(effects)+1,[-0.0,0.5,1.0,1.5],[],"Hedges' g","")
    d.add(Line(scale(0,-.2,1.5,fx,fx+fw),fy,scale(0,-.2,1.5,fx,fx+fw),fy+fh,strokeColor=GREY,strokeDashArray=[4,4]))
    for i,(lab,val) in enumerate(effects):
        yy0=fy+fh-(i+1)*42
        add_text(d,fx-15,yy0-4,lab,13,INK,"end")
        d.add(Circle(scale(val,-.2,1.5,fx,fx+fw),yy0,6,fillColor=BLUE if i<4 else ORANGE,strokeColor=WHITE))
        add_text(d,scale(val,-.2,1.5,fx,fx+fw)+12,yy0-4,f"{val:.2f}",11,INK)
    my=fy+20
    d.add(Line(scale(prow.ci_low,-.2,1.5,fx,fx+fw),my,scale(prow.ci_high,-.2,1.5,fx,fx+fw),my,strokeColor=RED,strokeWidth=3))
    d.add(Rect(scale(prow.meta_g,-.2,1.5,fx,fx+fw)-6,my-6,12,12,fillColor=RED,strokeColor=None))
    add_text(d,fx-15,my-4,"Random-effects meta",13,RED,"end",True)

    # D: prioritization matrix
    panel_label(d, 805, 505, "D")
    panel_box(d, 840, 75, 700, 440)
    attrs=[("4/4 direction",True),("FDR < 0.05",prow.fdr<.05),("I2 < 50%",prow.I2<50),("|g| >= 0.5",abs(prow.meta_g)>=.5),("PAH vs IPF-PH",prow.specificity_same_direction),("GSE53408",prow.GSE53408_direction_concordant),("SMC localization",True),("IPAH replication",True),("SSc-PAH replication",True)]
    add_text(d,890,465,"PDE8B evidence checklist",20,NAVY,bold=True)
    for i,(lab,ok) in enumerate(attrs):
        yy0=425-i*36
        d.add(Circle(910,yy0,9,fillColor=TEAL if ok else RED,strokeColor=None))
        add_text(d,910,yy0-5,"OK" if ok else "X",8,WHITE,"middle",True)
        add_text(d,935,yy0-5,lab,14,INK)
    add_text(d,1210,425,"Interpretation",18,NAVY,bold=True)
    lines=["Stable association", "PAH-biased signal", "SMC-associated", "Independent replication", "Causality not established"]
    for i,t in enumerate(lines):
        c=RED if "not" in t else TEAL
        add_text(d,1210,390-i*48,t,15,c,bold=True)
    save(d,"figure2_bulk_meta_and_pde8b")


def fig3(pb, ce):
    d=Drawing(W,H); d.add(Rect(0,0,W,H,fillColor=colors.HexColor("#F8FAFB"),strokeColor=None))
    title(d,"Figure 3. Pulmonary-artery single-cell transcriptomic overview",
          "GSE210248: subject-aware summary of 22,704 cells from 3 PAH and 3 donor lungs")
    base=pb[["sample","disease","cell_type","n_cells"]].drop_duplicates()
    counts=base.groupby(["sample","disease"],as_index=False).n_cells.sum().sort_values(["disease","sample"])
    samples=counts["sample"].tolist()

    # A counts
    panel_label(d,35,960,"A"); panel_box(d,70,600,710,375)
    bx,by,bw,bh=145,670,580,245; vmax=counts.n_cells.max()*1.08
    axes(d,bx,by,bw,bh,0,len(samples),0,vmax,[],[0,2000,4000,6000,8000],"","Cells")
    for i,row in counts.reset_index(drop=True).iterrows():
        c=RED if row.disease=="PAH" else BLUE
        x=bx+(i+.18)*bw/len(samples); ww=.64*bw/len(samples); hh=row.n_cells/vmax*bh
        d.add(Rect(x,by,ww,hh,fillColor=c,strokeColor=None))
        add_text(d,x+ww/2,by-23,row["sample"].replace("Donor_","D").replace("PAH_","P"),12,INK,"middle")
        add_text(d,x+ww/2,by+hh+8,f"{int(row.n_cells):,}",10,c,"middle",True)
    add_text(d,435,620,"D: donor   P: PAH",12,GREY,"middle")

    # B composition
    panel_label(d,795,960,"B"); panel_box(d,830,600,710,375)
    piv=base.pivot_table(index="sample",columns="cell_type",values="n_cells",aggfunc="sum",fill_value=0).reindex(samples)
    prop=piv.div(piv.sum(axis=1),axis=0)
    ctypes=prop.sum().sort_values(ascending=False).index.tolist()
    pal=[NAVY,BLUE,TEAL,ORANGE,PURPLE,RED,colors.HexColor("#70A288"),colors.HexColor("#C49A6C"),colors.HexColor("#7C8DA5"),colors.HexColor("#C97B9B"),colors.HexColor("#9B8E5E"),colors.HexColor("#5F9EA0"),colors.HexColor("#A0A0A0"),colors.HexColor("#4C6A92")]
    sx,sy,sw,sh=890,670,380,245
    for i,s in enumerate(samples):
        x=sx+i*sw/len(samples); bottom=sy
        for j,ct in enumerate(ctypes):
            hh=prop.loc[s,ct]*sh
            d.add(Rect(x,bottom,.72*sw/len(samples),hh,fillColor=pal[j%len(pal)],strokeColor=WHITE,strokeWidth=.2))
            bottom+=hh
        add_text(d,x+.36*sw/len(samples),sy-22,s.replace("Donor_","D").replace("PAH_","P"),11,INK,"middle")
    add_text(d,sx-15,sy-4,"0",11,GREY,"end");add_text(d,sx-15,sy+sh-4,"100%",11,GREY,"end")
    for j,ct in enumerate(ctypes[:10]):
        xx=1305+(j//5)*125; yy=900-(j%5)*42
        d.add(Rect(xx,yy,13,13,fillColor=pal[j],strokeColor=None));add_text(d,xx+19,yy,ct,10.5,INK)

    # C dot plot
    panel_label(d,35,545,"C"); panel_box(d,70,70,870,490)
    genes=["PDE8B","PIEZO2","DCBLD1","NCAM1","SELP","ABCG2"]
    celltypes=["SMC 1","SMC 2","Fibro","Endo 1","Endo 2","Mono / Macs","T cells","NK cells"]
    sub=ce[ce.gene.isin(genes)&ce.cell_type.isin(celltypes)].copy()
    gx,gy,cw,ch=240,145,100,42
    for j,g in enumerate(genes): add_text(d,gx+j*cw,515,g,13,RED if g=="PDE8B" else NAVY,"middle",True)
    for i,ct in enumerate(celltypes):
        yy0=475-i*ch; add_text(d,gx-65,yy0-4,ct,13,INK,"end")
        for j,g in enumerate(genes):
            r=sub[(sub.gene==g)&(sub.cell_type==ct)]
            if r.empty: continue
            rr=r.iloc[0]; det=max(rr.pah_detection,rr.donor_detection); rad=3+17*math.sqrt(max(0,det))
            d.add(Circle(gx+j*cw,yy0,rad,fillColor=color_interp(rr.delta_log_cp10k,-.4,.4),strokeColor=WHITE,strokeWidth=.6))
    add_text(d,275,102,"Dot size: detection fraction",12,GREY)
    for k,v in enumerate([.02,.10,.30]):
        rad=3+17*math.sqrt(v); xx=470+k*115;d.add(Circle(xx,107,rad,fillColor=colors.HexColor("#B0BEC5"),strokeColor=WHITE));add_text(d,xx+24,102,f"{v:.0%}",11,GREY)
    add_text(d,740,102,"Color: PAH - donor expression",12,GREY)
    for k,v in enumerate([-.3,0,.3]):d.add(Rect(855+k*30,98,28,16,fillColor=color_interp(v,-.4,.4),strokeColor=None))

    # D PDE8B subject points
    panel_label(d,955,545,"D"); panel_box(d,990,70,550,490)
    p=pb[(pb.gene=="PDE8B")&pb.cell_type.isin(["SMC 1","SMC 2"])].copy()
    add_text(d,1265,515,"PDE8B subject-level pseudobulk",18,NAVY,"middle",True)
    groups=[("SMC 1","Donor"),("SMC 1","PAH"),("SMC 2","Donor"),("SMC 2","PAH")]
    px,py,pw,ph=1060,150,410,300; ymax=max(.55,p.log_cp10k.max()*1.2)
    axes(d,px,py,pw,ph,0,4,0,ymax,[],[0,.1,.2,.3,.4,.5],"","log(CP10K + 1)")
    rng=np.random.default_rng(8)
    for i,(ct,dis) in enumerate(groups):
        vals=p[(p.cell_type==ct)&(p.disease==dis)].log_cp10k.values
        x=px+(i+.5)*pw/4
        if len(vals):
            mean=float(np.mean(vals)); d.add(Line(x-24,scale(mean,0,ymax,py,py+ph),x+24,scale(mean,0,ymax,py,py+ph),strokeColor=RED if dis=="PAH" else BLUE,strokeWidth=3))
        for j,v in enumerate(vals):
            jitter=(-1+2*(j/(max(1,len(vals)-1))))*13
            d.add(Circle(x+jitter,scale(v,0,ymax,py,py+ph),6,fillColor=RED if dis=="PAH" else BLUE,strokeColor=WHITE))
        add_text(d,x,py-25,"D" if dis=="Donor" else "P",12,INK,"middle",True)
        if dis=="Donor": add_text(d,x+pw/8,py-48,ct,12,NAVY,"middle")
    add_text(d,1265,93,"Each dot is one subject; horizontal line is the mean",11,GREY,"middle")
    save(d,"figure3_scrna_overview_and_localization")


def fig4(pri,pb,ce,rep,spec):
    d=Drawing(W,H);d.add(Rect(0,0,W,H,fillColor=colors.HexColor("#F8FAFB"),strokeColor=None))
    title(d,"Figure 4. Integrated evidence supporting an SMC-associated PDE8B axis",
          "Bulk robustness, disease context, single-cell localization, and independent subtype replication")
    prow=pri[pri.symbol.eq("PDE8B")].iloc[0]

    # A evidence effects
    panel_label(d,35,960,"A");panel_box(d,70,555,730,420)
    entries=[("GSE15197",prow.g_GSE15197,BLUE),("GSE113439",prow.g_GSE113439,BLUE),("GSE254617",prow.g_GSE254617,BLUE),("GSE208592",prow.g_GSE208592,BLUE),("Bulk meta",prow.meta_g,RED),("PAH vs IPF-PH",prow.g_PAH_vs_IPF_PH_GSE15197,ORANGE),("GSE53408",prow.g_GSE53408_secondary,PURPLE)]
    reps=rep[rep.gene.eq("PDE8B")]
    for _,r in reps.iterrows():entries.append((r.contrast.replace("_vs_Donor",""),r.hedges_g,TEAL))
    x0,y0,w,h=310,620,430,300; axes(d,x0,y0,w,h,-.2,1.6,0,len(entries)+1,[0,.5,1,1.5],[],"Standardized effect (Hedges' g)","")
    zx=scale(0,-.2,1.6,x0,x0+w);d.add(Line(zx,y0,zx,y0+h,strokeColor=GREY,strokeDashArray=[4,4]))
    for i,(lab,val,c) in enumerate(entries):
        yy=y0+h-(i+1)*29
        add_text(d,x0-15,yy-4,lab,12,INK,"end",lab=="Bulk meta")
        d.add(Circle(scale(val,-.2,1.6,x0,x0+w),yy,6,fillColor=c,strokeColor=WHITE))
        add_text(d,scale(val,-.2,1.6,x0,x0+w)+12,yy-4,f"{val:.2f}",10,c)
    add_text(d,105,580,"All standardized-effect comparisons are positive",12,GREY)

    # B cell localization
    panel_label(d,815,960,"B");panel_box(d,850,555,690,420)
    pce=ce[ce.gene.eq("PDE8B")].sort_values("delta_log_cp10k",ascending=False)
    x0,y0,w,h=1110,630,370,285; xmin=min(-.05,pce.delta_log_cp10k.min()*1.1);xmax=max(.12,pce.delta_log_cp10k.max()*1.15)
    show=pce.head(10)
    axes(d,x0,y0,w,h,xmin,xmax,0,len(show)+1,[-.05,0,.05,.1],[],"PAH - donor log(CP10K + 1)","")
    d.add(Line(scale(0,xmin,xmax,x0,x0+w),y0,scale(0,xmin,xmax,x0,x0+w),y0+h,strokeColor=GREY,strokeDashArray=[3,3]))
    for i,(_,r) in enumerate(show.iterrows()):
        yy=y0+h-(i+1)*25
        c=RED if r.cell_type.startswith("SMC") else GREY
        add_text(d,x0-12,yy-4,r.cell_type,11,c,"end",r.cell_type.startswith("SMC"))
        d.add(Circle(scale(r.delta_log_cp10k,xmin,xmax,x0,x0+w),yy,4+18*math.sqrt(max(r.pah_detection,r.donor_detection)),fillColor=c,fillOpacity=.75,strokeColor=WHITE))
    add_text(d,900,580,"Largest positive localization: SMC 2; concordant signal in SMC 1",12,GREY)

    # C disease-context contrast
    panel_label(d,35,505,"C");panel_box(d,70,75,730,440)
    sp=spec[spec.gene.eq("PDE8B")].iloc[0]
    add_text(d,435,465,"Disease-context direction",19,NAVY,"middle",True)
    items=[("PAH pulmonary artery",sp.PAH_vs_Donor,RED),("PHPF-PH pulmonary artery",sp["PHPF-PH_vs_Donor"],BLUE)]
    x0,y0,w,h=210,250,470,150; mx=max(abs(v) for _,v,_ in items)*1.35
    d.add(Line(scale(0,-mx,mx,x0,x0+w),y0-20,scale(0,-mx,mx,x0,x0+w),y0+h+20,strokeColor=INK))
    for i,(lab,val,c) in enumerate(items):
        yy=y0+h-i*85
        xz=scale(0,-mx,mx,x0,x0+w); xv=scale(val,-mx,mx,x0,x0+w)
        d.add(Rect(min(xz,xv),yy-14,abs(xv-xz),28,fillColor=c,strokeColor=None))
        add_text(d,x0-20,yy-4,lab,13,INK,"end")
        add_text(d,xv+(10 if val>0 else -10),yy-4,f"{val:+.3f}",12,c,"start" if val>0 else "end",True)
    add_text(d,435,170,"Opposite directions suggest a PAH-biased rather than universal PH signal",13,GREY,"middle")
    add_text(d,435,125,"Scale is within-cohort pseudobulk delta, not a cross-study standardized effect",11,GREY,"middle")

    # D model
    panel_label(d,815,505,"D");panel_box(d,850,75,690,440)
    add_text(d,1195,465,"Working model (association, not causality)",19,NAVY,"middle",True)
    boxes=[(900,330,160,80,"PAH vascular\nremodeling",NAVY),(1110,330,160,80,"SMC-associated\nPDE8B increase",RED),(1320,330,160,80,"Local cAMP\nregulation",ORANGE)]
    for x,y,w0,h0,lab,c in boxes:
        d.add(Rect(x,y,w0,h0,rx=10,ry=10,fillColor=colors.Color(c.red,c.green,c.blue,alpha=.15),strokeColor=c,strokeWidth=2))
        for k,line in enumerate(lab.split("\n")):add_text(d,x+w0/2,y+49-k*22,line,14,c,"middle",True)
    for x1,x2 in [(1060,1110),(1270,1320)]:
        d.add(Line(x1,370,x2-8,370,strokeColor=GREY,strokeWidth=2));d.add(Polygon([x2-8,376,x2,370,x2-8,364],fillColor=GREY,strokeColor=None))
    outcomes=["contractile tone","SMC proliferation","phenotypic switching","cell-cell signaling"]
    for i,t in enumerate(outcomes):
        x=925+(i%2)*280;y=220-(i//2)*85
        d.add(Rect(x,y,235,55,rx=8,ry=8,fillColor=WHITE,strokeColor=GRID))
        add_text(d,x+117,y+20,t,13,INK,"middle")
        d.add(Line(1400,330,x+117,y+55,strokeColor=GRID,strokeWidth=1))
    add_text(d,1195,105,"Required next: protein localization + PDE8B perturbation + cAMP/PKA readout",12,RED,"middle",True)
    save(d,"figure4_integrated_pde8b_evidence")


def main():
    meta,pri,pb,ce,rep,spec=load()
    fig1(meta,pri)
    fig2(meta,pri)
    fig3(pb,ce)
    fig4(pri,pb,ce,rep,spec)
    print(f"Generated SVG/PDF figures in {OUT}")


if __name__ == "__main__":
    main()



