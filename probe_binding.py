#!/usr/bin/env python3
"""Compare probe alignments to selected contigs and their top BLAST references."""
from __future__ import annotations
import argparse, html, re, subprocess, tempfile
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

IUPAC_DNA=set('ACGTRYSWKMBDHVN')
IUPAC_COMPAT = {
    'A': set('A'), 'C': set('C'), 'G': set('G'), 'T': set('T'),
    'R': set('AG'), 'Y': set('CT'), 'S': set('GC'), 'W': set('AT'),
    'K': set('GT'), 'M': set('AC'), 'B': set('CGT'), 'D': set('AGT'),
    'H': set('ACT'), 'V': set('ACG'), 'N': set('ACGT'),
}
def iupac_match(a, b):
    return bool(IUPAC_COMPAT.get(a.upper(), set()) & IUPAC_COMPAT.get(b.upper(), set()))
def revcomp(seq):
    return seq.upper().translate(str.maketrans(
        'ACGTRYSWKMBDHVN', 'TGCAYRSWMKVHDBN'))[::-1]

def read_fasta(path):
    d={}; name=None; seq=[]
    with open(path,encoding='utf-8') as f:
        for line in f:
            s=line.strip()
            if not s: continue
            if s.startswith('>'):
                if name is not None: d[name]=''.join(seq).upper()
                name=s[1:].split()[0]; seq=[]
            else:
                if name is None: raise ValueError(f'Sequence before FASTA header: {path}')
                seq.append(s)
    if name is not None: d[name]=''.join(seq).upper()
    return d

def norm_acc(s):
    s=str(s or '').strip(); parts=[p for p in s.split('|') if p]
    for x in parts+[s]:
        if re.match(r'^(?:[A-Z]{1,6}_)?[A-Z0-9]+\.\d+$',x): return x
    return parts[-1] if parts else s

def safe(s): return re.sub(r'[^A-Za-z0-9_.-]','_',str(s)) or 'unnamed'
def write_fasta(d,p):
    with open(p,'w') as f:
        for n,s in d.items():
            f.write(f'>{n}\n')
            for i in range(0,len(s),80): f.write(s[i:i+80]+'\n')
def cigar_stats(c):
    ops=[(int(n),op) for n,op in re.findall(r'(\d+)([MIDNSHP=X])',c)]
    q=sum(n for n,o in ops if o in 'MI=X'); r=sum(n for n,o in ops if o in 'MDN=X')
    cols=sum(n for n,o in ops if o in 'MID=X'); ins=sum(n for n,o in ops if o=='I'); dels=sum(n for n,o in ops if o=='D')
    return q,r,cols,ins,dels


def binding_class(longest):
    if longest > 70: return 'strong', 'green'
    if longest >= 60: return 'intermediate', 'orange'
    return 'weak', 'red'

def map_sam(probes_fa,target_fa,kind,sample,threads):
    cmd=['minimap2','-ax','sr','--MD','-N','50','-t',str(threads),str(target_fa),str(probes_fa)]
    proc=subprocess.run(cmd,check=True,text=True,capture_output=True)
    qseqs=read_fasta(probes_fa); tseqs=read_fasta(target_fa); out=[]
    for line in proc.stdout.splitlines():
        if not line or line.startswith('@'): continue
        f=line.split('\t')
        if len(f)<11: continue
        q,flag_s,target,pos_s,mapq_s,cigar=f[:6]; flag=int(flag_s)
        if flag&4 or target=='*' or cigar=='*': continue
        tags={}
        for t in f[11:]:
            z=t.split(':',2)
            if len(z)==3: tags[z[0]]=z[2]
        nm=int(tags.get('NM','0'))
        qaln,raln,cols,ins,dels=cigar_stats(cigar)
        if not cols: continue
        qseq=qseqs.get(q,''); tseq=tseqs.get(target,'')
        if not tseq: continue
        start=int(pos_s) -1; reverse=bool(flag&16)
        aln=reconstruct_alignment(qseq,tseq,start,cigar,reverse)
        longest=aln['longest_match_run']
        classification,colour=binding_class(longest)
        out.append(dict(
            sample=sample, probe_id=q, target_type=kind, target_id=target,
            strand='-' if reverse else '+', start=start+1, end=start+raln,
            probe_length=len(qseq), aligned_query_bases=qaln,
            probe_fraction_aligned=round(qaln/len(qseq),4) if qseq else 0,
            probe_coverage_pct=round(100*qaln/len(qseq),3) if qseq else 0,
            alignment_columns=cols, match_bases=aln['match_bases'],
            mismatch_bases=aln['mismatch_bases'],
            identity_pct=round(100*aln['match_bases']/cols,3) if cols else 0,
            longest_match_run=longest, binding_class=classification,
            plot_colour=colour, edit_distance=nm, inserted_bases=ins,
            deleted_bases=dels, mapq=int(mapq_s), cigar=cigar,
            md=tags.get('MD',''), alignment_probe=aln['alignment_probe'],
            alignment_match=aln['alignment_match'],
            alignment_target=aln['alignment_target'],
            match_runs=aln['match_runs'], is_secondary=bool(flag&256)))
    return out

def depth(path):
    xs=[]; ys=[]
    if not Path(path).is_file(): return xs,ys
    for line in open(path):
        f=line.rstrip().split('\t')
        if len(f)>=3:
            try: xs.append(int(f[1])); ys.append(int(f[2]))
            except ValueError: pass
    return xs,ys

def reconstruct_alignment(probe_seq, target_seq, start0, cigar, reverse):
    """
    Reconstruct a minimap2 SAM alignment.

    probe_seq is the original probe sequence.
    target_seq is the full forward-oriented target sequence.
    start0 is the 0-based SAM POS.
    reverse indicates FLAG 16.

    Matching is IUPAC-compatible rather than exact character equality.
    """

    query = revcomp(probe_seq) if reverse else probe_seq

    ops = [
        (int(n), op)
        for n, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    ]

    qpos = 0
    rpos = start0

    aligned_probe = []
    aligned_target = []
    match_line = []

    match_runs = []
    current_run_start = None
    current_run_length = 0

    match_bases = 0
    mismatch_bases = 0

    def finish_run():
        nonlocal current_run_start, current_run_length

        if current_run_length:
            match_runs.append(
                (
                    current_run_start,
                    current_run_start + current_run_length,
                    current_run_length,
                )
            )

        current_run_start = None
        current_run_length = 0

    for n, op in ops:

        # Query-only clipping.
        if op == 'S':
            qpos += n

            if qpos > len(query):
                raise ValueError(
                    f"CIGAR {cigar}: soft clipping moves query position "
                    f"to {qpos}, but query length is {len(query)}"
                )

            continue

        # Hard clipping / padding consume neither sequence.
        if op in ('H', 'P'):
            continue

        # Insertion relative to reference.
        if op == 'I':
            if qpos + n > len(query):
                raise ValueError(
                    f"CIGAR {cigar}: insertion consumes query through "
                    f"{qpos+n}, but query length is {len(query)}"
                )

            finish_run()

            for i in range(n):
                aligned_probe.append(query[qpos])
                aligned_target.append('-')
                match_line.append(' ')
                qpos += 1

            continue

        # Deletion / skipped reference sequence.
        if op in ('D', 'N'):
            finish_run()

            for i in range(n):
                if rpos >= len(target_seq):
                    raise ValueError(
                        f"CIGAR {cigar}: {op} consumes reference beyond "
                        f"target length: rpos={rpos}, "
                        f"target_length={len(target_seq)}, "
                        f"start0={start0}"
                    )

                aligned_probe.append('-')
                aligned_target.append(target_seq[rpos])
                match_line.append(' ')
                rpos += 1

            continue

        # Aligned query/reference columns.
        if op in ('M', '=', 'X'):

            if qpos + n > len(query):
                raise ValueError(
                    f"CIGAR {cigar}: {op} consumes query through "
                    f"{qpos+n}, but query length is {len(query)} "
                    f"(probe length={len(probe_seq)}, reverse={reverse})"
                )

            if rpos + n > len(target_seq):
                raise ValueError(
                    f"CIGAR {cigar}: {op} consumes reference through "
                    f"{rpos+n}, but target length is {len(target_seq)} "
                    f"(start0={start0}, target_length={len(target_seq)}, "
                    f"reverse={reverse})"
                )

            for i in range(n):
                q = query[qpos]
                t = target_seq[rpos]

                compatible = iupac_match(q, t)

                aligned_probe.append(q)
                aligned_target.append(t)

                if compatible:
                    match_line.append('.')
                    match_bases += 1

                    if current_run_length == 0:
                        current_run_start = rpos

                    current_run_length += 1

                else:
                    match_line.append('.')
                    mismatch_bases += 1
                    finish_run()

                qpos += 1
                rpos += 1

            continue

        raise ValueError(
            f"Unsupported CIGAR operation {op} in {cigar}"
        )

    finish_run()

    # The CIGAR should account for the complete query sequence,
    # apart from hard clipping.
    if qpos != len(query):
        raise ValueError(
            f"CIGAR/query length mismatch: CIGAR consumed {qpos} "
            f"query bases, but query length is {len(query)}; "
            f"CIGAR={cigar}, reverse={reverse}"
        )

    return {
        'alignment_probe': ''.join(aligned_probe),
        'alignment_match': ''.join(match_line),
        'alignment_target': ''.join(aligned_target),
        'match_bases': match_bases,
        'mismatch_bases': mismatch_bases,
        'longest_match_run': (
            max((x[2] for x in match_runs), default=0)
        ),
        'match_runs': match_runs,
    }
def main():
    p=argparse.ArgumentParser()
    for name in ['sample','probes','reference-fasta','ranked-hits','selected-dir','reference-coverage','analysis-dir','outdir']:
        p.add_argument('--'+name,required=True)
    p.add_argument('--threads',type=int,default=2); a=p.parse_args()
    probes=read_fasta(a.probes)
    if not probes: raise ValueError('No probes found')

    IUPAC_DNA = set("ACGTRYSWKMBDHVN")

    bad = []

    for name, seq in probes.items():
        invalid = sorted(set(seq) - IUPAC_DNA)
        if not seq:
            bad.append(f"{name}: empty sequence")
        elif invalid:
            bad.append(f"{name}: invalid characters {invalid}")

    if bad:
        raise ValueError(
            "Invalid probe sequence(s): " + "; ".join(bad)
        )

    #bad=[k for k,v in probes.items() if not v or set(v)-DNA]
    #if bad: raise ValueError(f'Invalid probe sequence(s): {bad}')
    ranked=pd.read_csv(a.ranked_hits,sep='\t',dtype=str).fillna('')
    required={'genes','contig','top_accession'}
    if not required.issubset(ranked.columns): raise ValueError(f'Missing ranked-hit columns: {required-set(ranked.columns)}')
    ranked=ranked.drop_duplicates('genes',keep='first'); ranked=ranked[(ranked.genes!='')&(ranked.contig!='')]
    ref_all=read_fasta(a.reference_fasta); ref_index={}
    for h,s in ref_all.items():
        acc=norm_acc(h); ref_index[h]=s; ref_index[acc]=s; ref_index[re.sub(r'\.\d+$','',acc)]=s
    refs={}; ref_meta={}; selected={}; sel_meta={}; missing=[]
    for _,r in ranked.iterrows():
        gene,contig=r['genes'],r['contig']; acc=norm_acc(r['top_accession'])
        seq=ref_index.get(acc) or ref_index.get(re.sub(r'\.\d+$','',acc))
        if seq:
            tid=safe(gene)+'__'+safe(acc); refs[tid]=seq; ref_meta[tid]=(gene,acc)
        else: missing.append(f'{gene}\t{r["top_accession"]}')
        sf=Path(a.selected_dir)/(gene+'.fa')
        if not sf.is_file(): sf=Path(a.selected_dir)/(safe(gene)+'.fa')
        if sf.is_file():
            recs=read_fasta(sf)
            if recs:
                tid=safe(gene)+'__'+safe(contig); selected[tid]=next(iter(recs.values()))
                sel_meta[tid]=(gene,contig)
    if not refs: raise ValueError('No top BLAST reference accessions could be resolved in reference FASTA')
    if not selected: raise ValueError(f'No selected contig FASTAs found in {a.selected_dir}')
    out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='probe_binding_',dir=out) as td:
        td=Path(td); qf=td/'probes.fa'; rf=td/'refs.fa'; sf=td/'selected.fa'
        write_fasta(probes,qf); write_fasta(refs,rf); write_fasta(selected,sf)
        rh=map_sam(qf,rf,'reference',a.sample,a.threads); sh=map_sam(qf,sf,'selected_contig',a.sample,a.threads)
    for h in rh: h['gene'],h['reference_accession']=ref_meta[h['target_id']]; h['contig']=''
    for h in sh: h['gene'],h['contig']=sel_meta[h['target_id']]; h['reference_accession']=''
    hits=rh+sh
    cols=['sample','probe_id','gene','target_type','target_id','reference_accession','contig','strand','start','end','probe_length','aligned_query_bases','probe_fraction_aligned','probe_coverage_pct','alignment_columns','match_bases','mismatch_bases','identity_pct','longest_match_run','binding_class','plot_colour','edit_distance','inserted_bases','deleted_bases','mapq','cigar','md','alignment_probe','alignment_match','alignment_target','match_runs','is_secondary']
    pd.DataFrame(hits,columns=cols).to_csv(out/'probe_hits.tsv',sep='\t',index=False)
    cov=pd.read_csv(a.reference_coverage,sep='\t') if Path(a.reference_coverage).is_file() else pd.DataFrame()
    with PdfPages(out/'probe_binding.pdf') as pdf:
        for _,r in ranked.iterrows():
            gene,contig=r['genes'],r['contig']; rgh=[h for h in rh if h['gene']==gene]; shh=[h for h in sh if h['gene']==gene]
            if not rgh and not shh: continue
            fig,axs=plt.subplots(
                2,1,figsize=(13,6),sharex=False,constrained_layout=True
            )
            from matplotlib.patches import Patch

            def draw_probe_track(ax, target_hits, target_length,
                                 covered_positions=None, mark_uncovered=False):
                # Draw a positional track below the x-axis.
                # Coloured regions are contiguous IUPAC-compatible matching runs.
                # For selected contigs, uncovered positions are black.
                track_y = -0.16
                track_h = 0.10

                if mark_uncovered and target_length:
                    # Black is the baseline: it represents zero/missing read coverage.
                    ax.add_patch(
                        plt.Rectangle(
                            (1, track_y),
                            target_length,
                            track_h,
                            transform=ax.get_xaxis_transform(),
                            facecolor='black',
                            edgecolor='none',
                            alpha=.9,
                            clip_on=False,
                        )
                    )

                def covered_parts(start0, end0):
                    if not mark_uncovered or covered_positions is None:
                        yield start0, end0
                        return
                    # Only colour bases having positive read depth.
                    run_start = None
                    prev = None
                    for pos in range(start0 + 1, end0 + 1):
                        if pos in covered_positions:
                            if run_start is None:
                                run_start = pos
                            prev = pos
                        elif run_start is not None:
                            yield run_start, prev + 1
                            run_start = None
                            prev = None
                    if run_start is not None:
                        yield run_start, prev + 1

                for h in target_hits:
                    for rs, run_end, runlen in h['match_runs']:
                        _, colour = binding_class(runlen)
                        for part_start, part_end in covered_parts(rs, run_end):
                            ax.add_patch(
                                plt.Rectangle(
                                    (part_start, track_y),
                                    part_end - part_start,
                                    track_h,
                                    transform=ax.get_xaxis_transform(),
                                    facecolor=colour,
                                    edgecolor='none',
                                    alpha=.9,
                                    clip_on=False,
                                )
                            )

                ax.text(
                    0.0, track_y - 0.035,
                    'Probe-compatible binding',
                    transform=ax.transAxes,
                    ha='left', va='top', fontsize=8
                )

                if target_length:
                    ax.set_xlim(1, target_length)

            # Top: reference depth plus positional probe-compatible track.
            acc=norm_acc(r['top_accession'])
            ax=axs[0]
            reference_length=0
            if not cov.empty and {'accession','position','depth'}.issubset(cov.columns):
                cv=cov[cov['accession'].astype(str).map(norm_acc)==acc]
                if cv.empty and 'gene' in cov.columns:
                    cv=cv[cv['gene'].astype(str)==gene]
                if not cv.empty:
                    ax.plot(cv['position'],cv['depth'],linewidth=.65)
                    reference_length=int(cv['position'].max())

            if not reference_length:
                reference_length=max((h['end'] for h in rgh),default=0)

            draw_probe_track(ax,rgh,reference_length)
            ax.set_title(f'{gene} — top reference ({acc})')
            ax.set_xlabel('Reference position (bp)')
            ax.set_ylabel('Read depth')
            ax.grid(alpha=.25)

            # Bottom: selected-contig depth plus positional probe-compatible track.
            ax=axs[1]
            xs,ys=depth(Path(a.analysis_dir)/f'{contig}.fa.depth.txt')
            if xs:
                ax.plot(xs,ys,linewidth=.65)

            selected_length=max((h['end'] for h in shh),default=0)
            if xs:
                selected_length=max(selected_length,max(xs))

            # A position is considered covered when the depth file reports >0.
            # Missing positions are treated as uncovered, so they remain black.
            covered_positions={x for x,y in zip(xs,ys) if y > 0}
            draw_probe_track(
                ax, shh, selected_length,
                covered_positions=covered_positions,
                mark_uncovered=True,
            )
            ax.set_title(f'{gene} — selected contig ({contig})')
            ax.set_xlabel('Selected-contig position (bp)')
            ax.set_ylabel('Read depth')
            ax.grid(alpha=.25)

            fig.legend(handles=[
                Patch(facecolor='green',alpha=.9,label='>20 bp contiguous match'),
                Patch(facecolor='orange',alpha=.9,label='10–20 bp contiguous match'),
                Patch(facecolor='red',alpha=.9,label='<10 bp contiguous match'),
                Patch(facecolor='black',alpha=.9,label='No contig read coverage'),
            ],loc='upper center',ncol=3,fontsize=8)
            pdf.savefig(fig); plt.close(fig)
    df=pd.DataFrame(hits,columns=cols)
    if not df.empty: df=df.sort_values(['gene','probe_id','target_type','start'])
    table=df.to_html(index=False,escape=True,table_id='hits')
    blocks=[]
    for _,h in df.iterrows():
        blocks.append(
            f"""<details class="alignment" data-search="{html.escape(str(h['probe_id'])+' '+str(h['gene'])+' '+str(h['target_type'])+' '+str(h['contig']))}">
<summary><b>Probe: {html.escape(str(h['probe_id']))}</b> — {html.escape(str(h['gene']))} — {html.escape(str(h['target_type']))} — longest compatible run: <b>{int(h['longest_match_run'])} bp</b> ({html.escape(str(h['binding_class']))})</summary>
<div class="meta">Coverage: {float(h['probe_coverage_pct']):.1f}% | Identity: {float(h['identity_pct']):.1f}% | Matches: {int(h['match_bases'])} | Mismatches: {int(h['mismatch_bases'])} | Indels: {int(h['inserted_bases'])}I/{int(h['deleted_bases'])}D | Strand: {html.escape(str(h['strand']))} | CIGAR: {html.escape(str(h['cigar']))}</div>
<pre class="seq">PROBE  {html.escape(str(h['alignment_probe']))}
        {html.escape(str(h['alignment_match']))}
TARGET  {html.escape(str(h['alignment_target']))}</pre>
</details>"""
        )
    align_html=''.join(blocks)
    page=f"""<!doctype html><html><head><meta charset="utf-8"><title>Probe binding — {html.escape(str(a.sample))}</title>
<style>
body{{font-family:Arial,sans-serif;margin:2rem}}
table{{border-collapse:collapse;font-size:.8rem;width:100%}}
td,th{{border:1px solid #ddd;padding:.35rem;white-space:nowrap}}
th{{position:sticky;top:0;background:#eee}}
.wrap{{overflow:auto;max-height:60vh}}
.alignment{{border:1px solid #ddd;border-radius:6px;margin:.5rem 0;padding:.5rem}}
.alignment summary{{cursor:pointer}}
.seq{{overflow:auto;background:#f7f7f7;padding:1rem;font-family:monospace;line-height:1.4}}
.meta{{margin:.5rem 0;color:#555}}
.legend span{{display:inline-block;padding:.25rem .5rem;margin-right:.5rem;border-radius:3px}}
.green{{background:#90ee90}} .orange{{background:#ffd580}} .red{{background:#ff9999}} .black{{background:#000;color:#fff}}
</style></head><body>
<h1>Probe binding — {html.escape(str(a.sample))}</h1>
<p>Probe alignments against each gene's top BLAST reference and corresponding selected contig.
IUPAC-compatible bases are counted as matches. Classification is based on the
<b>longest contiguous compatible matching run</b>, not overall identity.</p>
<div class="legend"><span class="green">&gt;20 bp contiguous match</span>
<span class="orange">10–20 bp contiguous match</span>
<span class="red">&lt;10 bp contiguous match</span>
<span style="background:#000;color:#fff">No contig read coverage</span></div>
<p><a href="probe_binding.pdf">Coverage and probe PDF</a> | <a href="probe_hits.tsv">Download alignment TSV</a></p>
<input id="q" placeholder="Filter probe, gene, target..." style="width:30rem;padding:.5rem">
<h2>Alignment table</h2><div class="wrap">{table}</div>
<h2>Actual probe alignments</h2>{align_html}
<script>
document.getElementById('q').addEventListener('input',e=>{{
 let q=e.target.value.toLowerCase();
 document.querySelectorAll('#hits tbody tr').forEach(r=>r.style.display=r.textContent.toLowerCase().includes(q)?'':'none');
 document.querySelectorAll('.alignment').forEach(r=>r.style.display=r.dataset.search.toLowerCase().includes(q)?'':'none');
}})
</script></body></html>"""
    (out/'probe_binding.html').write_text(page,encoding='utf-8')
    if missing: (out/'missing_reference_accessions.tsv').write_text('gene\taccession\n'+'\n'.join(missing)+'\n')
    print(f'Wrote {out}/probe_hits.tsv, probe_binding.pdf, probe_binding.html')
if __name__=='__main__': main()
