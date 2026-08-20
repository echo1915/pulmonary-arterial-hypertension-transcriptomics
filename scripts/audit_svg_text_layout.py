import re, xml.etree.ElementTree as ET
from pathlib import Path

ROOT=Path(r"G:\workdata\projects\project-001-pulmonary-arterial-hypertension-transcriptomics\outputs\submission-draft-figures")
NS={'s':'http://www.w3.org/2000/svg'}
def main():
 rows=[]
 for f in sorted(ROOT.glob('*.svg')):
  root=ET.parse(f).getroot(); vb=list(map(float,root.attrib['viewBox'].split())); W,H=vb[2:]
  boxes=[]; out=[]
  for t in root.findall('.//s:text',NS):
   txt=''.join(t.itertext()).strip(); st=t.attrib.get('style',''); tr=t.attrib.get('transform','')
   m=re.search(r'font-size:\s*([0-9.]+)px',st); z=re.search(r'(?:translate|matrix)\(([^)]+)\)',tr)
   if not txt or not m: continue
   fs=float(m.group(1)); nums=[float(x) for x in re.split(r'[ ,]+',z.group(1).strip())]
   x,y=(nums[-2],nums[-1]) if len(nums)>=2 else (0,0); width=max(fs*.45*len(txt),fs); height=fs*1.2
   anchor='middle' if 'text-anchor: middle' in st else ('end' if 'text-anchor: end' in st else 'start')
   l=x-width/2 if anchor=='middle' else (x-width if anchor=='end' else x); b=y-height; box=(l,b,l+width,y,txt,fs)
   boxes.append(box)
   if l<-.5 or b<-.5 or l+width>W+.5 or y>H+.5: out.append(txt)
  overlaps=[]
  for i,a in enumerate(boxes):
   for b in boxes[i+1:]:
    inter=max(0,min(a[2],b[2])-max(a[0],b[0]))*max(0,min(a[3],b[3])-max(a[1],b[1]))
    if inter>0.30*min((a[2]-a[0])*(a[3]-a[1]),(b[2]-b[0])*(b[3]-b[1])) and abs(a[5]-b[5])<8:
     if a[4]!=b[4]: overlaps.append((a[4],b[4]))
  rows.append((f.name,len(boxes),len(out),len(overlaps),out[:5],overlaps[:5]))
 for r in rows: print(r)
if __name__=='__main__': main()
