#!/usr/bin/env python3
import math
from pathlib import Path

import pandas as pd
from reportlab.graphics import renderPDF, renderSVG
from reportlab.graphics.shapes import Circle, Drawing, Line, Rect, String
from reportlab.lib import colors

ROOT = Path(r"G:\workdata\projects\project-001-pulmonary-arterial-hypertension-transcriptomics\outputs\pde8b-coexpression")
OUT = ROOT / "figures"
W, H = 1600, 1050
NAVY, BLUE, TEAL = colors.HexColor("#17324D"), colors.HexColor("#2F6B9A"), colors.HexColor("#2A9D8F")
ORANGE, RED, GREY = colors.HexColor("#E6953B"), colors.HexColor("#C94C4C"), colors.HexColor("#75808A")
GRID, INK, BG = colors.HexColor("#D8E0E5"), colors.HexColor("#263238"), colors.HexColor("#F8FAFB")

def text(d,x,y,s,size=15,c=INK,anchor="start",bold=False):
    d.add(String(x,y,str(s),fontName="Helvetica-Bold" if bold else "Helvetica",fontSize=size,fillColor=c,textAnchor=anchor))

def box(d,x,y,w,h): d.add(Rect(x,y,w,h,rx=9,ry=9,fillColor=colors.white,strokeColor=GRID))
def sx(v,a,b,x,w): return x+(v-a)/(b-a)*w
def sy(v,a,b,y,h): return y+(v-a)/(b-a)*h

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    m=pd.read_csv(ROOT/"pde8b_cross_cohort_coexpression_meta.csv")
    e=pd.read_csv(ROOT/"enrichment"/"pde8b_core_pathway_enrichment.csv")
    d=Drawing(W,H); d.add(Rect(0,0,W,H,fillColor=BG,strokeColor=None))
    text(d,45,1000,"PDE8B defines a reproducible vascular contractile co-expression program",28,NAVY,bold=True)
    text(d,45,968,"Four lung cohorts; group-adjusted and PAH-only concordance; custom 15,015-gene enrichment background",15,GREY)
    # A volcano
    text(d,40,915,"A",23,NAVY,bold=True); box(d,75,515,720,410)
    x,y,w,h=135,580,610,285
    d.add(Line(x,y,x+w,y,strokeColor=INK)); d.add(Line(x,y,x,y+h,strokeColor=INK))
    m["mlog"]=-m.adjusted_fdr.clip(lower=1e-300).map(math.log10)
    bg=m[~m.mechanism_module.astype(bool)].sample(min(1800,(~m.mechanism_module.astype(bool)).sum()),random_state=7)
    for _,r in bg.iterrows(): d.add(Circle(sx(r.adjusted_meta_r,-.65,.65,x,w),sy(min(r.mlog,12),0,12,y,h),1.3,fillColor=GRID,strokeColor=None))
    for _,r in m[m.mechanism_module.astype(bool)].iterrows():
        c=RED if bool(r.core_module) else ORANGE
        d.add(Circle(sx(r.adjusted_meta_r,-.65,.65,x,w),sy(min(r.mlog,12),0,12,y,h),2.4,fillColor=c,strokeColor=None))
    for g in ["TRPC3","PRKG1","PDE5A","PRKAA2","PRKAB2","PLCB4"]:
        r=m[m.symbol.eq(g)].iloc[0]; px=sx(r.adjusted_meta_r,-.65,.65,x,w); py=sy(min(r.mlog,12),0,12,y,h)
        text(d,px+5,py+4,g,10,NAVY,bold=True)
    for v in [-.6,-.3,0,.3,.6]: text(d,sx(v,-.65,.65,x,w),y-22,f"{v:.1f}",11,GREY,"middle")
    for v in [0,3,6,9,12]: text(d,x-10,sy(v,0,12,y,h)-4,str(v),11,GREY,"end")
    text(d,x+w/2,y-48,"Group-adjusted meta-correlation with PDE8B",13,INK,"middle")
    text(d,x,y+h+14,"-log10(FDR)",12,INK)
    text(d,115,535,"543 stable genes; 201 high-confidence core (|r| >= 0.40)",13,NAVY,bold=True)
    # B enrichment
    text(d,820,915,"B",23,NAVY,bold=True); box(d,855,515,700,410)
    wanted=["cilium organization","muscle system process","actin filament-based process","muscle contraction","cGMP-PKG signaling pathway","smooth muscle contraction","cellular response to calcium ion"]
    q=e[e.term_name.isin(wanted)].drop_duplicates("term_name").sort_values("fdr",ascending=False)
    bx,by,bw,bh=1120,570,370,290
    vmax=max(4.0,float((-q.fdr.map(math.log10)).max())+.3)
    for i,(_,r) in enumerate(q.iterrows()):
        py=by+i*38; score=-math.log10(r.fdr)
        text(d,bx-15,py-4,r.term_name,12,INK,"end")
        d.add(Line(bx,py,bx+bw,py,strokeColor=GRID,strokeWidth=.6))
        rad=4+min(10,float(r.intersection_size)/2)
        d.add(Circle(sx(score,0,vmax,bx,bw),py,rad,fillColor=TEAL,strokeColor=colors.white))
    for v in range(0,math.ceil(vmax)+1): text(d,sx(v,0,vmax,bx,bw),by-35,str(v),11,GREY,"middle")
    text(d,bx+bw/2,by-62,"-log10(FDR); dot size = contributing genes",12,INK,"middle")
    # C mechanism chain
    text(d,40,455,"C",23,NAVY,bold=True); box(d,75,80,1480,385)
    nodes=[("PDE8B","cAMP hydrolysis",RED),("TRPC3 / Ca2+","Ca2+ entry",ORANGE),("PLCB4","Ca2+ mobilization",ORANGE),("PRKG1 / PDE5A","cGMP-PKG",TEAL),("PPP1R12B / CALD1","contractile apparatus",BLUE)]
    for i,(a,b,c) in enumerate(nodes):
        nx=120+i*285
        d.add(Rect(nx,260,225,95,rx=13,ry=13,fillColor=colors.Color(c.red,c.green,c.blue,alpha=.13),strokeColor=c,strokeWidth=2))
        text(d,nx+112,315,a,15,c,"middle",True); text(d,nx+112,285,b,12,GREY,"middle")
        if i<4:
            d.add(Line(nx+225,307,nx+275,307,strokeColor=NAVY,strokeWidth=2));
            d.add(Line(nx+265,313,nx+275,307,strokeColor=NAVY,strokeWidth=2)); d.add(Line(nx+265,301,nx+275,307,strokeColor=NAVY,strokeWidth=2))
    text(d,120,195,"Interpretation",14,NAVY,bold=True)
    text(d,120,165,"PDE8B is a reproducible marker of a smooth-muscle-enriched cyclic-nucleotide / calcium / contractile state.",16,INK)
    text(d,120,130,"Co-expression does not prove that PDE8B causally controls this program; cell-composition and cilium signals remain explicit alternatives.",14,GREY)
    stem=OUT/"Figure_5_PDE8B_coexpression_and_pathways"
    renderSVG.drawToFile(d,str(stem.with_suffix('.svg'))); renderPDF.drawToFile(d,str(stem.with_suffix('.pdf')))
    print(stem)

if __name__=="__main__": main()
