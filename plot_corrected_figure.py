#!/usr/bin/env python3
"""Regenerate corrected Figure 3 from package tables. Requires matplotlib."""
from pathlib import Path
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parents[1];F=P/'figures';F.mkdir(exist_ok=True)
def load(p):
 with p.open(encoding='utf-8-sig') as f:return list(csv.DictReader(f))
win=load(P/'tables/TermSeq_Local_Windows.csv')
metrics=load(P/'tables/TermSeq_Per_Replicate_Boundary_Metrics.csv')
summ=load(P/'tables/TermSeq_Candidate_Boundary_Summary.csv')
mp=load(P/'tables/Candidate_Coordinate_Mapping.csv')
neighborhood=load(P/'inputs/Genomic_Neighborhood_JKD6008.csv');neighbor=neighborhood
evidence=load(P/'tables/FINAL_3UTR_Candidate_Evidence_Integration.csv')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
fig=plt.figure(figsize=(11.5,10.6));gs=fig.add_gridspec(3,2,height_ratios=[1.05,1.25,1.05],hspace=.50,wspace=.27)
for j,(short,name) in enumerate([('RS04120','3UTR-RS04120'),('RS05585','3UTR-RS05585')]):
 e=next(e for e in evidence if e['candidate_short']==name);ax=fig.add_subplot(gs[0,j]);rr=[r for r in neighbor if r['candidate']==name]
 center=(int(e['candidate_start'])+int(e['candidate_end']))/2
 for r in rr:
  a=int(r['genomic_start']);b=int(r['genomic_end']);strand=r['strand'];y=.64 if strand=='+' else .32
  if b<center-3900 or a>center+4600: continue
  start=a if strand=='+' else b;dx=b-a if strand=='+' else a-b
  ax.arrow(start,y,dx,0,width=.07,head_width=.12,head_length=min(110,abs(dx)*.25),length_includes_head=True,color='#555555')
  gene=r['gene_symbol']
  if gene not in ['','NA']:ax.text((a+b)/2,y+.12 if strand=='+' else y-.14,gene,ha='center',va='center',fontsize=7,fontstyle='italic',rotation=35 if j==1 and strand=='+' else 0)
 ax.axvspan(int(e['candidate_start']),int(e['candidate_end']),color='#b22222',alpha=.7)
 ax.annotate(short,xy=(center,.75),xytext=(center,.99),ha='center',fontsize=8,fontweight='bold',arrowprops=dict(arrowstyle='-',lw=.8))
 ax.set_xlim(center-3900,center+4600);ax.set_ylim(.02,1.15);ax.set_yticks([]);ax.spines['left'].set_visible(False);ax.ticklabel_format(axis='x',style='plain',useOffset=False);ax.tick_params(axis='x',labelsize=7)
 ax.set_xlabel('JKD6008 genomic coordinate (bp)',fontsize=8)
 ax.set_title(f'{"ab"[j]}  Genomic context of {short}',loc='left',fontweight='bold',fontsize=10,pad=16)
 ax.text(0,1.03,f"{e['distance_from_parent_3prime_nt']} nt from the annotated {e['putative_parent_gene']} 3′ end",transform=ax.transAxes,fontsize=8)
 ax=fig.add_subplot(gs[1,j]);
 for rep,ls,col in [('rep1','-','#222222'),('rep2','--','#687a89')]:
  rr=sorted([r for r in win if r['candidate']==name and r['replicate']==rep],key=lambda r:int(r['position']))
  ax.plot([int(r['relative_position']) for r in rr],[float(r['signal']) for r in rr],ls=ls,color=col,lw=1.1,label=rep)
 ax.axvspan(-5,5,color='#e2e6eb');ax.axvline(0,color='#777777',lw=.7,ls=':');ax.set_xlim(-100,100);ax.set_ylim(bottom=0)
 ax.set_xlabel('Position relative to mapped 3′ boundary (nt)',fontsize=8);ax.set_ylabel('JKD6009 Term-seq 3′-end signal',fontsize=8)
 ax.set_title(f'{"cd"[j]}  Term-seq profile of {short}',loc='left',fontweight='bold',fontsize=10,pad=24)
 m=next(x for x in mp if x['candidate']==name)
 ax.text(0,1.035,f"Mapped JKD6009 boundary: {int(m['target_3prime']):,} ({m['strand']})",transform=ax.transAxes,fontsize=8)
 ax.legend(frameon=False,fontsize=7,loc='upper right')
 if j==0:ax.annotate('Replicated peak\n+2 nt',xy=(2,1034),xytext=(30,820),fontsize=8,arrowprops=dict(arrowstyle='->',lw=.8))
 else:ax.annotate('Dominant local peak\n−86 nt',xy=(-86,411),xytext=(-48,335),fontsize=8,arrowprops=dict(arrowstyle='->',lw=.8))
ax=fig.add_subplot(gs[2,:]);ax.axis('off');ax.set_title('e  Independent evidence layers',loc='left',fontsize=10,fontweight='bold',pad=12)
labels=['Putative parental transcript (JKD6008)','Distance from parental 3′ end','CLASH hybrids / experiments','nrdF target log₂FC / adjusted P','Mapped JKD6009 candidate 3′ end','Boundary-window maximum: rep1 / rep2','Strong replicated peak near boundary','Independent processing in JKD6008']
vals=[['nrdF','folD'],['53 nt','177 nt'],['15 / 2','7 / 1'],['−1.118 / 3.04 × 10⁻¹¹','−1.118 / 3.04 × 10⁻¹¹'],['812384 (+)','1085084 (−)'],['849 / 1034','3 / 3'],['Yes: +2 nt in both replicates','No: dominant peak at −86 nt'],['Not established','Not established']]
table=ax.table(cellText=[[labels[i],*vals[i]] for i in range(len(labels))],colLabels=['Evidence','RS04120','RS05585'],colWidths=[.43,.285,.285],cellLoc='center',bbox=[0,-.02,1,.96]);table.auto_set_font_size(False);table.set_fontsize(8)
for (r,c),cell in table.get_celld().items():
 cell.set_edgecolor('#dddddd');cell.set_linewidth(.35)
 if r==0:cell.set_facecolor('#e2e6eb');cell.set_text_props(weight='bold')
 elif r%2==0:cell.set_facecolor('#f6f7f8')
 if c==0:cell.set_text_props(ha='left')
fig.subplots_adjust(left=.08,right=.98,bottom=.055,top=.95)
fig.savefig(F/'Fig3_TermSeq_corrected.png',dpi=600);fig.savefig(F/'Fig3_TermSeq_corrected.pdf');plt.close(fig)
