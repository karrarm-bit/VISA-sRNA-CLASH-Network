#!/usr/bin/env python3
"""Re-extract JKD6009 Term-seq windows after sequence-verified coordinate mapping.
Run: python scripts/correct_termseq.py --project-root /path/to/Bac-sRNA --output corrected_results
Python standard library only. Original raw tracks are not modified.
"""
import argparse,csv,gzip,hashlib,io,json,re,statistics,sys
from pathlib import Path
from zipfile import ZipFile
P=argparse.ArgumentParser();P.add_argument('--project-root',type=Path,required=True);P.add_argument('--output',type=Path,required=True);A=P.parse_args()
A.output.mkdir(parents=True,exist_ok=True)
package=Path(__file__).resolve().parents[1]
def get(rel,archive):
 for base in [A.project_root,A.project_root/'upload']:
  p=base/rel
  if p.is_file():return p.read_bytes()
  zp=base/archive
  if zp.is_file():
   with ZipFile(zp) as z:return z.read(rel)
 raise FileNotFoundError(rel)
def csvout(name,rows):
 with (A.output/name).open('w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def fasta(path):
 lines=path.read_text().splitlines();return ''.join(lines[1:]).upper()
s=fasta(package/'references/JKD6008_CP002120.1.fa');t=fasta(package/'references/JKD6009_NZ_LR027876.1.fa')
assert len(s)==2924344 and len(t)==2930439
specs=[('3UTR-RS04120',811819,811847,'+',812356,812384),('3UTR-RS05585',1084557,1084613,'-',1085084,1085140)]
mapping=[]
for name,a,b,strand,c,d in specs:
 q=s[a-51:b+50];rc=q.translate(str.maketrans('ACGTN','TGCAN'))[::-1];p=t.find(q)
 assert p==c-51 and t.find(q,p+1)==-1 and t.find(rc)==-1,(name,'Nonunique or changed reference mapping')
 assert s[a-1:b]==t[c-1:d]
 mapping.append(dict(candidate=name,source_reference='CP002120.1 (JKD6008; equivalent NC_017341.1 sequence)',source_start=a,source_end=b,target_reference='NZ_LR027876.1 (JKD6009)',target_start=c,target_end=d,strand=strand,coordinate_system='1-based inclusive',source_3prime=b if strand=='+' else a,target_3prime=d if strand=='+' else c,offset=c-a,flank_each_side_nt=50,exact_alignment_length_nt=len(q),unique_forward_match=True,reverse_complement_match=False))
csvout('Candidate_Coordinate_Mapping.csv',mapping)
csvout('Candidate_Coordinates_JKD6009.csv',[dict(candidate=m['candidate'],start=m['target_start'],end=m['target_end'],strand=m['strand'],expected_3prime=m['target_3prime'],reference='NZ_LR027876.1') for m in mapping])
windows=[];oldwindows=[];provenance=[]
for rep,gsm in [('rep1','GSM4811618'),('rep2','GSM4811619')]:
 for strand,word in [('+','plus'),('-','minus')]:
  name=f'{gsm}_Staph_termseq_cont{1 if rep=="rep1" else 2}_{word}_strand_single_3end_base_only.gr.gz'
  path='07_processed_data/downloaded_files/GSE158830/'+name
  blob=get(path,'07_processed_data.zip');provenance.append(dict(path=path,sha256=hashlib.sha256(blob).hexdigest(),sample=gsm,replicate=rep,strand=strand))
  regions=[]
  for m in mapping:
   if m['strand']==strand:
    regions.extend([(m,m['target_3prime'],windows),(m,m['source_3prime'],oldwindows)])
  with gzip.GzipFile(fileobj=io.BytesIO(blob)) as f:
   for line in f:
    chrom,pos,signal=line.split()[:3];pos=int(pos)
    for m,bound,dest in regions:
     if abs(pos-bound)<=100:
      delta=pos-bound
      dest.append(dict(chromosome=chrom.decode(),position=pos,signal=float(signal),replicate=rep,strand=strand,source_file=name,candidate=m['candidate'],expected_3prime=bound,relative_position=delta,region='Expected boundary ±5 nt' if abs(delta)<=5 else 'Local background' if abs(delta)>=20 else 'Intermediate region'))
assert len(windows)==len(oldwindows)==804
assert all(r['chromosome']=='Saa_6009_C1' for r in windows)
def metrics(rows):
 result=[]
 for m in mapping:
  for rep in ['rep1','rep2']:
   rr=sorted([r for r in rows if r['candidate']==m['candidate'] and r['replicate']==rep],key=lambda r:r['position']);assert len(rr)==201
   near=[r['signal'] for r in rr if abs(r['relative_position'])<=5];bg=[r['signal'] for r in rr if 20<=abs(r['relative_position'])<=100];peak=max(rr,key=lambda r:r['signal'])
   ma=max(near);med=statistics.median(bg);ratio=(ma+1)/(med+1);pn=abs(peak['relative_position'])<=5
   support='Strong local 3′-end signal' if ma>0 and ratio>=5 and pn else 'Local 3′-end enrichment' if ma>0 and ratio>=2 else '3′-end signal detected' if ma>0 else 'No 3′-end signal detected'
   result.append(dict(candidate=m['candidate'],replicate=rep,strand=m['strand'],expected_3prime=rr[0]['expected_3prime'],boundary_max_signal=ma,boundary_total_signal=sum(near),background_median_signal=med,background_mean_signal=statistics.mean(bg),strongest_position=peak['position'],strongest_relative_position=peak['relative_position'],strongest_signal=peak['signal'],boundary_to_background_ratio=ratio,peak_near_expected_boundary=pn,replicate_boundary_support=support))
 return result
met=metrics(windows);oldmet=metrics(oldwindows)
# Validate this implementation against the supplied original per-replicate results.
try:
 original=list(csv.DictReader(io.StringIO(get('17A_TermSeq_3UTR_boundary_audit/tables/TermSeq_Per_Replicate_Boundary_Metrics.csv','17A_TermSeq_3UTR_boundary_audit.zip').decode('utf-8-sig'))))
 for r in oldmet:
  ref=next(x for x in original if x['candidate']==r['candidate'] and x['replicate']==r['replicate'])
  for k in ['boundary_max_signal','boundary_total_signal','background_median_signal','background_mean_signal','strongest_position','strongest_relative_position','strongest_signal','boundary_to_background_ratio']:
   assert abs(float(r[k])-float(ref[k]))<1e-9,(r['candidate'],r['replicate'],k,r[k],ref[k])
 validation='PASS: recalculation at old coordinates reproduces every original numeric boundary metric for all four candidate-replicate combinations'
except FileNotFoundError:validation='Original audit not supplied; comparison not run'
summary=[]
for m in mapping:
 rr=[r for r in met if r['candidate']==m['candidate']];signals=sum(r['boundary_max_signal']>0 for r in rr);enrich=sum(r['boundary_max_signal']>0 and r['boundary_to_background_ratio']>=2 for r in rr);strong=sum(r['boundary_max_signal']>0 and r['boundary_to_background_ratio']>=5 and r['peak_near_expected_boundary'] for r in rr)
 category='Reproducible strong 3′-boundary support' if strong>=2 else 'Reproducible local 3′-boundary enrichment' if enrich>=2 else 'Reproducible 3′-end signal without strong enrichment' if signals>=2 else 'Single-replicate 3′-end signal' if signals==1 else 'No Term-seq boundary signal detected'
 interp='Reproducible local boundary support in parental JKD6009; independent processing in JKD6008 remains unestablished.' if enrich>=2 else 'No reproducible discrete candidate boundary established in parental JKD6009; independent processing in JKD6008 remains unestablished.'
 summary.append(dict(candidate=m['candidate'],strand=m['strand'],expected_3prime=m['target_3prime'],replicates_examined=2,replicates_with_signal=signals,replicates_with_enrichment=enrich,replicates_with_strong_local_peak=strong,median_boundary_signal=statistics.median(r['boundary_max_signal'] for r in rr),median_enrichment_ratio=statistics.median(r['boundary_to_background_ratio'] for r in rr),TermSeq_boundary_class=category,manuscript_interpretation=interp))
csvout('TermSeq_Local_Windows.csv',windows);csvout('TermSeq_Per_Replicate_Boundary_Metrics.csv',met);csvout('TermSeq_Candidate_Boundary_Summary.csv',summary);csvout('Original_Coordinate_Metrics_Reproduced.csv',oldmet);csvout('Source_Track_Checksums.csv',provenance)
report=dict(method='Same formulas and thresholds as archived Script17A; coordinates changed after unique exact sequence mapping, not optimized against signal',python_version=sys.version,validation=validation,coordinate_mapping=mapping,corrected_summary=summary,scope='Candidate-boundary audit only. RNA-seq, CLASH, genomic context and IntaRNA were not rerun. Parental-strain untreated Term-seq is contextual evidence, not strain-matched or vancomycin-dependent validation.')
(A.output/'Correction_Audit.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report,indent=2))
