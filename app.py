import streamlit as st
import primer3
import pandas as pd
import io
import time
import plotly.graph_objects as go
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import Entrez

# NCBI asks that you identify yourself. Put YOUR real address here --
# NCBI blocks anonymous / fake addresses that hammer the server.
Entrez.email = "spud_researcher@example.com"

# Stable output schema. Every row -- success or error -- carries all of these
# keys, so a failed gene can never delete the Tm columns from the table.
DISPLAY_COLS = [
    'Gene_ID', 'Organism', 'Rank',
    'Forward_Seq', 'F_Len', 'F_Tm', 'F_GC',
    'Reverse_Seq', 'R_Len', 'R_Tm', 'R_GC',
    'Mean_Tm', 'Tm_Diff', 'Tm_Flag',
    'Amp_Len', 'Primer_Fold', 'Penalty',
    'BLAST_Flag', 'Status',
]


# --- Small numeric helpers -------------------------------------------------
def r1(x, nd=1):
    """Round only if numeric; None stays None instead of raising TypeError."""
    return round(x, nd) if isinstance(x, (int, float)) else None


def gc_pct(seq):
    if not seq:
        return None
    s = seq.upper()
    return round(100.0 * (s.count('G') + s.count('C')) / len(s), 1)


def rev_comp(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement.get(base, base) for base in reversed(seq.upper()))


def blank_row(gene_id, status, organism=None):
    """An error/placeholder row that still has the full column schema."""
    row = {c: None for c in DISPLAY_COLS}
    row['Gene_ID'] = gene_id
    row['Organism'] = organism
    row['Status'] = status
    return row


# --- BLAST -----------------------------------------------------------------
def build_entrez_query(species):
    """
    Accepts a scientific name ('Solanum tuberosum') or a taxid
    ('4113' / 'txid4113'). Taxid is safer: a misspelled name is a valid
    Entrez query that matches nothing, which looks identical to a
    genuinely unique primer.
    """
    s = (species or "").strip()
    if not s:
        return ""
    if s.lower().startswith("txid") and s[4:].strip().isdigit():
        return f"txid{s[4:].strip()}[ORGN]"
    if s.isdigit():
        return f"txid{s}[ORGN]"
    return f'"{s}"[organism]'


def blast_count_subjects(primer_seq, species, min_identity=0.95, min_cov=0.9):
    """
    Counts DISTINCT subject accessions carrying a near-full-length,
    high-identity, 3'-anchored alignment to the primer.

    Counting distinct accessions (not raw HSPs) matters: S. tuberosum has
    several assemblies plus RefSeq mRNA and genomic records in nt, so one
    unique primer routinely produces many HSPs of the same locus.

    word_size=7 + permissive expect are required for a ~20-mer to be
    reported at all; this mirrors NCBI Primer-BLAST's own settings.

    Returns {'n': int|None, 'error': str|None, 'accessions': [str]}
    """
    try:
        handle = NCBIWWW.qblast(
            "blastn", "nt", primer_seq,
            entrez_query=build_entrez_query(species),
            word_size=7, expect=1000, hitlist_size=50,
        )
        record = NCBIXML.read(handle)
        primer_len = len(primer_seq)
        subjects = set()
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                if hsp.align_length < primer_len * min_cov:
                    continue
                if hsp.identities / hsp.align_length < min_identity:
                    continue
                # The 3' end drives extension -- an alignment that stops
                # short of it will not prime, regardless of overall identity.
                if hsp.query_end < primer_len:
                    continue
                subjects.add(alignment.accession)
                break
        return {"n": len(subjects), "error": None, "accessions": sorted(subjects)[:5]}
    except Exception as e:
        return {"n": None, "error": type(e).__name__, "accessions": []}


def label_hits(res):
    """Human-readable label + machine flag. 0 hits is NOT a pass."""
    if res["error"]:
        return f"BLAST Error: {res['error']}", "ERROR"
    n = res["n"]
    if n == 0:
        return "0 Hits ⚠️ unverified", "UNVERIFIED"
    if n == 1:
        return "1 Hit ✅", "OK"
    return f"{n} Hits ⚠️", "MULTI"


def blast_primer_pair(forward_seq, reverse_seq, species, cache):
    """BLASTs both primers, with a cache so repeated primers cost one call."""
    def cached(seq):
        key = (seq, species)
        if key not in cache:
            cache[key] = blast_count_subjects(seq, species)
            time.sleep(2)  # polite delay; NCBI rate-limits aggressively
        return cache[key]

    f_res, r_res = cached(forward_seq), cached(reverse_seq)
    f_label, f_flag = label_hits(f_res)
    r_label, r_flag = label_hits(r_res)

    # Worst flag wins.
    order = {"OK": 0, "MULTI": 1, "UNVERIFIED": 2, "ERROR": 3}
    pair_flag = f_flag if order[f_flag] >= order[r_flag] else r_flag
    return f"F: {f_label} | R: {r_label}", pair_flag


# --- Target parsing & primer geometry --------------------------------------
def parse_target(target_str):
    """Parses the Smart Input into Global, Junction, or ROI (all 1-based)."""
    target_str = str(target_str).strip()
    if not target_str or target_str.lower() in ['nan', 'none']:
        return "Global", None
    if '-' in target_str:
        try:
            start, end = map(int, target_str.split('-'))
            return "ROI", [start, end]
        except Exception:
            return "Error", None
    if target_str.isdigit():
        return "Junction", int(target_str)
    return "Error", None


def validate_7bp_anchor(primer_start, primer_length, junction, is_reverse=False):
    """At least 7 bp of the primer on each side of the junction."""
    if is_reverse:
        r_end = primer_start - primer_length + 1
        return (primer_start >= junction + 7) and (r_end <= junction - 6)
    f_end = primer_start + primer_length - 1
    return (primer_start <= junction - 7) and (f_end >= junction + 6)


def extract_candidate(raw_res, i, gene_id):
    f_left = raw_res.get(f'PRIMER_LEFT_{i}')
    f_right = raw_res.get(f'PRIMER_RIGHT_{i}')
    return {
        "Gene_ID": gene_id,
        "f_start": f_left[0], "f_len": f_left[1],
        "r_start": f_right[0], "r_len": f_right[1],
        "Forward_Seq": raw_res.get(f'PRIMER_LEFT_{i}_SEQUENCE'),
        "Reverse_Seq": raw_res.get(f'PRIMER_RIGHT_{i}_SEQUENCE'),
        "F_Tm": raw_res.get(f'PRIMER_LEFT_{i}_TM'),
        "R_Tm": raw_res.get(f'PRIMER_RIGHT_{i}_TM'),
        "Amp_Len": raw_res.get(f'PRIMER_PAIR_{i}_PRODUCT_SIZE'),
        "Penalty": raw_res.get(f'PRIMER_PAIR_{i}_PENALTY'),
    }


def build_gene_map(gene_id, seq_len, candidates, junction_0based=None):
    """Amplicon track map; each bar is labelled with its Tm pair."""
    fig = go.Figure()
    fig.add_shape(type="rect", x0=0, y0=0, x1=seq_len, y1=1,
                  line=dict(color="gray"), fillcolor="lightgray")

    for cand in candidates:
        rank = cand['Rank']
        y = 1.5 + ((rank - 1) * 0.6)
        x0, x1 = cand['f_start'], cand['r_start']
        fig.add_shape(type="rect", x0=x0, y0=y - 0.2, x1=x1, y1=y + 0.2,
                      fillcolor="LimeGreen", opacity=0.8, line=dict(color="green"))
        tm_txt = f"Tm {cand['F_Tm']}/{cand['R_Tm']}°C"
        fig.add_annotation(x=(x0 + x1) / 2, y=y,
                           text=f"Rank {rank} · {cand['Amp_Len']} bp · {tm_txt}",
                           showarrow=False, font=dict(size=11))

    if junction_0based is not None:
        top_y = 1.5 + (len(candidates) * 0.6)
        fig.add_shape(type="line", x0=junction_0based, y0=-0.5,
                      x1=junction_0based, y1=top_y, line=dict(color="red", dash="dash"))
        fig.add_annotation(x=junction_0based, y=top_y + 0.2,
                           text="Junction", showarrow=False, font=dict(color="red"))

    fig.update_layout(
        title=f"Amplicon Map: {gene_id}",
        xaxis_title="Position (bp)", yaxis_visible=False,
        height=300 + (len(candidates) * 30), plot_bgcolor="white"
    )
    return fig


# --- Page Setup ------------------------------------------------------------
st.set_page_config(page_title="SPUD - Batch Mode", page_icon="🧬", layout="wide")

with st.sidebar:
    st.title("📖 SPUD Documentation")
    st.markdown("Welcome to **SPUD** (Specific Primer Universal Designer) — **Batch Edition**.")

    with st.expander("🎯 1. Smart Target Input"):
        st.write("""
        The **Target** field auto-detects your design strategy.
        **All positions are 1-based**:
        * **Single number (452):** Junction Mode — one primer spans the exon-exon
          junction with a strict 7-bp anchor on each side.
        * **Range (300-450):** ROI Mode — the whole amplicon sits inside the range.
        * **Empty:** Global Mode — best thermodynamic pairs anywhere.
        """)

    with st.expander("🌡️ 2. Reading the Tm columns"):
        st.write("""
        * **F_Tm / R_Tm:** primer3 Tm under *your* salt and primer concentrations
          (SantaLucia 1998 thermodynamics, salt-corrected).
        * **Mean_Tm / Tm_Diff:** pair mean, and the F–R gap. A gap above 2 °C is
          the practical problem — one primer out-competes the other at anneal.
        * **Tm_Flag:** `Outlier` if the pair mean is >1.5 °C from your Target Tm
          (won't share a plate-wide program), `PairGap` if F and R differ by >2 °C.
        """)

    with st.expander("🔍 3. Reading the BLAST column"):
        st.write("""
        `BLAST_Flag` reports **distinct subject accessions**, not raw HSPs:
        * **OK** — exactly one accession matched. What you want.
        * **MULTI** — several distinct accessions. Often assembly/RefSeq
          redundancy of the same locus rather than a true off-target: open the
          accessions before rejecting the pair.
        * **UNVERIFIED** — zero hits. This is *not* a clean result. The usual
          cause is a misspelled organism name or the gene being absent from `nt`.
          Verify the organism field before trusting the pair.
        * **ERROR** — NCBI call failed (rate limit, timeout).

        A per-primer hit count cannot predict amplification: an off-target
        product needs *both* primers on the same subject, in opposing
        orientation, within your amplicon window. Treat this as a screen.
        """)

    with st.expander("🧬 4. Organism field"):
        st.write("""
        Prefer the **taxid** (`4113` for *Solanum tuberosum*, `3760` for
        *Prunus persica*, `3702` for *Arabidopsis thaliana*). A typo in a
        scientific name is a valid Entrez query that matches nothing —
        producing a silent `0 Hits` that looks like success.

        Leaving the field blank searches all of `nt`, where any 20-mer hits
        thousands of unrelated genomes. Off-target priming can only occur on
        the template in your tube, so keep the filter on.
        """)

st.title("🧬 SPUD: Specific Primer Universal Designer (Batch Edition)")
st.divider()

with st.expander("⚙️ Lab Conditions & Advanced Thermodynamics (Applied to all)", expanded=True):
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        species_name = st.text_input(
            "Default Organism or taxid (fallback):", "txid4113",
            help="Used only when a gene has no organism of its own. Set it "
                 "per gene beside the Gene ID, or in an Organism column in "
                 "your batch file. taxid is safer than a name — a misspelling "
                 "silently returns 0 hits. S. tuberosum = 4113."
        )
        num_returns = st.slider("Max Candidates per Gene", 1, 10, 3, step=1)
    with col2:
        primer_conc = st.number_input("Primer Conc. (nM)", value=250.0, step=10.0)
        mg_conc = st.number_input("Mg2+ Conc. (mM)", value=2.5, step=0.1)
    with col3:
        target_tm = st.slider("Target Tm (°C)", 50.0, 72.0, 60.0, step=0.5)
        max_pair_gap = st.slider("Max F–R Tm gap (°C)", 0.5, 5.0, 2.0, step=0.5)
    with col4:
        min_amp = st.text_input("Min Amplicon Length", value="80")
        max_amp = st.text_input("Max Amplicon Length", value="150")
        blast_mode = st.radio(
            "🔍 BLAST Specificity Check:",
            ["Off", "Rank 1 Only", "All Candidates"],
            index=1,
            help="Rank 1 = top candidate per gene (fast). All = every candidate "
                 "(slow, risk of NCBI rate-limiting)."
        )

if blast_mode != "Off" and not species_name.strip():
    st.warning(
        "No default organism set. Any gene that doesn't specify its own will be "
        "BLASTed against all of `nt`, where every 20-mer looks non-specific."
    )

mode = st.radio("Select Processing Mode:",
                ["Batch Upload (CSV/Excel)", "Single Gene (Quick Test)"], horizontal=True)

genes_data = []

if mode == "Batch Upload (CSV/Excel)":
    st.info("Upload a file with columns: **Gene_ID**, **Sequence**, **Target**, "
            "and optionally **Organism** (name or taxid; blank rows fall back "
            "to the default above).")
    template_df = pd.DataFrame({
        "Gene_ID": ["StPOT1", "StGA2ox", "PpPDC1"],
        "Sequence": ["ATGC...", "TTGC...", "CCGG..."],
        "Target": ["452", "300-450", ""],
        "Organism": ["txid4113", "Solanum tuberosum", "txid3760"]
    })
    csv_buffer = io.BytesIO()
    template_df.to_csv(csv_buffer, index=False)
    st.download_button(label="📥 Download Template", data=csv_buffer.getvalue(),
                       file_name="SPUD_Template.csv", mime="text/csv")

    uploaded_file = st.file_uploader("Upload Batch File", type=['csv', 'xlsx'])
    if uploaded_file:
        try:
            df_input = (pd.read_csv(uploaded_file) if uploaded_file.name.endswith('.csv')
                        else pd.read_excel(uploaded_file))
            for index, row in df_input.iterrows():
                org = str(row.get("Organism", "") or "").strip()
                if org.lower() in ("nan", "none"):
                    org = ""
                genes_data.append({
                    "Gene_ID": str(row.get("Gene_ID", f"Gene_{index}")),
                    "Sequence": str(row.get("Sequence", "")),
                    "Target": str(row.get("Target", "")),
                    "Organism": org
                })
        except Exception as e:
            st.error(f"Error reading file: {e}")
else:
    col_a, col_b = st.columns([2, 1])
    with col_a:
        s_seq = st.text_area("Target Sequence (5' to 3'):", height=150)
    with col_b:
        s_id = st.text_input("Gene ID:", "MyGene")
        s_target = st.text_input("Target (Smart Input):", placeholder="e.g., 452 OR 300-450")
        s_org = st.text_input(
            "Organism (for BLAST):", "",
            placeholder=f"blank = {species_name or 'no filter'}",
            help="Scientific name or taxid for this gene. Leave blank to use "
                 "the default set in Lab Conditions."
        )
    if s_seq:
        genes_data.append({"Gene_ID": s_id, "Sequence": s_seq,
                           "Target": s_target, "Organism": s_org})


# --- Execution Engine ------------------------------------------------------
if st.button("🚀 Run SPUD Engine", type="primary") and genes_data:
    all_results = []
    gene_map_data = {}

    global_args = {
        'PRIMER_OPT_SIZE': 20, 'PRIMER_MIN_SIZE': 18, 'PRIMER_MAX_SIZE': 25,
        'PRIMER_OPT_TM': target_tm, 'PRIMER_MIN_TM': target_tm - 5.0,
        'PRIMER_MAX_TM': target_tm + 5.0,
        'PRIMER_DNA_CONC': primer_conc, 'PRIMER_SALT_DIVALENT': mg_conc,
        'PRIMER_SALT_MONOVALENT': 50.0, 'PRIMER_DNTP_CONC': 0.8,
        'PRIMER_TM_FORMULA': 1, 'PRIMER_SALT_CORRECTIONS': 1,
        'PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT': 1,
        'PRIMER_NUM_RETURN': max(20, num_returns * 5)
    }
    if min_amp.isdigit() and max_amp.isdigit():
        global_args['PRIMER_PRODUCT_SIZE_RANGE'] = [[int(min_amp), int(max_amp)]]

    progress_bar = st.progress(0)
    status_text = st.empty()

    for idx, gene in enumerate(genes_data):
        status_text.text(f"Processing {gene['Gene_ID']} ({idx+1}/{len(genes_data)})...")
        clean_seq = "".join(gene['Sequence'].split()).upper()

        # Per-gene organism wins; the Lab Conditions field is only a fallback.
        gene_org = (gene.get('Organism') or "").strip() or species_name.strip()

        if len(clean_seq) < 50:
            all_results.append(blank_row(gene['Gene_ID'], "Error: Sequence too short", gene_org))
            progress_bar.progress((idx + 1) / len(genes_data))
            continue

        target_type, target_val = parse_target(gene['Target'])
        seq_args = {'SEQUENCE_ID': gene['Gene_ID'], 'SEQUENCE_TEMPLATE': clean_seq}
        junction_0based = None

        if target_type == "Junction":
            if target_val < 20 or target_val > len(clean_seq) - 20:
                all_results.append(blank_row(gene['Gene_ID'], "Error: Junction too close to edge", gene_org))
                progress_bar.progress((idx + 1) / len(genes_data))
                continue
            junction_0based = target_val - 1  # 1-based UI -> 0-based primer3
            seq_args['SEQUENCE_OVERLAP_JUNCTION_LIST'] = [junction_0based]

        elif target_type == "ROI":
            roi_start = max(0, target_val[0] - 1)
            roi_len = target_val[1] - target_val[0] + 1
            if roi_len > 0 and roi_start + roi_len <= len(clean_seq):
                seq_args['SEQUENCE_INCLUDED_REGION'] = [roi_start, roi_len]
            else:
                all_results.append(blank_row(gene['Gene_ID'], "Error: Invalid ROI range", gene_org))
                progress_bar.progress((idx + 1) / len(genes_data))
                continue

        elif target_type == "Error":
            all_results.append(blank_row(gene['Gene_ID'], "Error: Invalid Target format", gene_org))
            progress_bar.progress((idx + 1) / len(genes_data))
            continue

        try:
            raw_res = primer3.bindings.designPrimers(seq_args, global_args)
            valid_candidates = []

            for i in range(raw_res.get('PRIMER_PAIR_NUM_RETURNED', 0)):
                cand = extract_candidate(raw_res, i, gene['Gene_ID'])

                # Skip malformed records rather than crashing the whole gene --
                # this is what used to wipe the Tm columns out of the table.
                if cand['Forward_Seq'] is None or cand['Reverse_Seq'] is None:
                    continue
                if not isinstance(cand['F_Tm'], (int, float)) or \
                   not isinstance(cand['R_Tm'], (int, float)):
                    continue

                if target_type == "Junction":
                    spans = (validate_7bp_anchor(cand['f_start'], cand['f_len'], junction_0based, False)
                             or validate_7bp_anchor(cand['r_start'], cand['r_len'], junction_0based, True))
                    if not spans:
                        continue

                # Hairpin of the primer itself competing with annealing
                f_hp = primer3.calc_hairpin(cand['Forward_Seq']).tm
                r_hp = primer3.calc_hairpin(cand['Reverse_Seq']).tm
                fold_flag = "!" if (f_hp >= cand['F_Tm'] - 3.0 or r_hp >= cand['R_Tm'] - 3.0) else "v"
                penalty = cand['Penalty'] if isinstance(cand['Penalty'], (int, float)) else 999.0
                if fold_flag == "!":
                    penalty += 50.0

                mean_tm = (cand['F_Tm'] + cand['R_Tm']) / 2.0
                tm_diff = abs(cand['F_Tm'] - cand['R_Tm'])

                row = blank_row(gene['Gene_ID'], "Pending", gene_org)
                row.update({
                    "Forward_Seq": cand['Forward_Seq'],
                    "F_Len": cand['f_len'],
                    "F_Tm": r1(cand['F_Tm']),
                    "F_GC": gc_pct(cand['Forward_Seq']),
                    "Reverse_Seq": cand['Reverse_Seq'],
                    "R_Len": cand['r_len'],
                    "R_Tm": r1(cand['R_Tm']),
                    "R_GC": gc_pct(cand['Reverse_Seq']),
                    "Mean_Tm": r1(mean_tm),
                    "Tm_Diff": r1(tm_diff),
                    "Amp_Len": cand['Amp_Len'],
                    "Primer_Fold": fold_flag,
                    "Penalty": r1(penalty, 2),
                })
                row["f_start"] = cand['f_start']
                row["r_start"] = cand['r_start']
                valid_candidates.append(row)

            valid_candidates.sort(key=lambda x: x['Penalty'])
            top_candidates = valid_candidates[:num_returns]

            for rank, row in enumerate(top_candidates):
                row['Rank'] = rank + 1
                all_results.append(row)

            if top_candidates:
                gene_map_data[gene['Gene_ID']] = {
                    "seq_len": len(clean_seq),
                    "junction": junction_0based,
                    "candidates": [dict(c) for c in top_candidates],
                }
            else:
                all_results.append(blank_row(gene['Gene_ID'],
                                             "No candidates found under constraints",
                                             gene_org))

        except Exception as e:
            all_results.append(blank_row(gene['Gene_ID'], f"Error: {e}", gene_org))

        progress_bar.progress((idx + 1) / len(genes_data))

    # --- Final Audit & BLAST ---
    status_text.text("Finalizing quality checks & BLAST...")
    blast_cache = {}

    for row in all_results:
        if row.get('Status') != "Pending":
            continue

        flags = []
        if abs(row['Mean_Tm'] - target_tm) > 1.5:
            flags.append("Outlier")
        if row['Tm_Diff'] > max_pair_gap:
            flags.append("PairGap")
        row['Tm_Flag'] = "/".join(flags) if flags else "OK"

        should_blast = (blast_mode == "All Candidates"
                        or (blast_mode == "Rank 1 Only" and row['Rank'] == 1))
        if should_blast:
            row_org = row.get('Organism') or ""
            status_text.text(
                f"BLASTing {row['Gene_ID']} (Rank {row['Rank']}) "
                f"against {row_org or 'all of nt'}..."
            )
            label, flag = blast_primer_pair(row['Forward_Seq'], row['Reverse_Seq'],
                                            row_org, blast_cache)
            row['Status'] = label
            row['BLAST_Flag'] = flag
        else:
            row['Status'] = "Skipped BLAST"
            row['BLAST_Flag'] = "NOT_RUN"

    status_text.text("Complete!")

    # Force the full schema so the Tm columns exist even in an all-error run.
    df_results = pd.DataFrame(all_results)
    for c in DISPLAY_COLS:
        if c not in df_results.columns:
            df_results[c] = None
    extra = [c for c in df_results.columns
             if c not in DISPLAY_COLS and c not in ('f_start', 'r_start')]
    df_results = df_results[DISPLAY_COLS + extra]

    st.session_state['spud_results'] = df_results
    st.session_state['spud_map_data'] = gene_map_data


# --- Results Display -------------------------------------------------------
if 'spud_results' in st.session_state:
    df_results = st.session_state['spud_results']
    gene_map_data = st.session_state.get('spud_map_data', {})

    st.subheader("📋 Results Table")
    st.dataframe(df_results, use_container_width=True)

    if 'BLAST_Flag' in df_results.columns:
        unver = df_results[df_results['BLAST_Flag'] == "UNVERIFIED"]
        if len(unver):
            bad_orgs = sorted({str(o) for o in unver['Organism'].dropna().unique()})
            st.warning(
                f"{len(unver)} pair(s) returned **0 BLAST hits** — an unverified "
                f"result, not a clean one. Check that these organism entries "
                f"resolve in NCBI taxonomy: {', '.join(f'`{o}`' for o in bad_orgs) or '(none set)'}"
            )

    st.download_button(
        label="📥 Download Results (CSV)",
        data=df_results.to_csv(index=False).encode('utf-8'),
        file_name="SPUD_Results.csv", mime="text/csv", type="primary"
    )

    if gene_map_data:
        st.divider()
        st.subheader("🗺️ Amplicon Maps")
        for gene_id, data in gene_map_data.items():
            st.plotly_chart(
                build_gene_map(gene_id, data["seq_len"], data["candidates"], data["junction"]),
                use_container_width=True
            )
