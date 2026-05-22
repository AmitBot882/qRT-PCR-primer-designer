import streamlit as st
import primer3
import pandas as pd
import io
import time
import plotly.graph_objects as go
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import Entrez

# Configuration for NCBI
Entrez.email = "spud_researcher@example.com"

# --- Helper Functions & Thermodynamics ---
def rev_comp(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement.get(base, base) for base in reversed(seq.upper()))

def check_and_get_blast(primer_seq, species):
    """Runs blastn for a single primer and summarizes significant hits."""
    try:
        entrez_query = f"{species}[organism]" if species else ""
        handle = NCBIWWW.qblast(
            "blastn", "nt", primer_seq,
            entrez_query=entrez_query, word_size=11, hitlist_size=3
        )
        blast_record = NCBIXML.read(handle)
        significant_hits = 0
        primer_len = len(primer_seq)
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.align_length >= primer_len * 0.9:
                    significant_hits += 1
        if significant_hits > 1:
            return f"{significant_hits} Hits"
        elif significant_hits == 1:
            return "1 Hit ✅"
        else:
            return "0 Hits ✅"
    except Exception as e:
        # Return an informative reason instead of a generic "Error"
        return f"BLAST Error: {type(e).__name__}"

def blast_primer_pair(forward_seq, reverse_seq, species):
    """BLASTs both Forward and Reverse primers and returns a combined status."""
    f_status = check_and_get_blast(forward_seq, species)
    time.sleep(2)  # polite delay to avoid NCBI rate-limiting / IP blocking
    r_status = check_and_get_blast(reverse_seq, species)
    time.sleep(2)
    return f"F: {f_status} | R: {r_status}"

def parse_target(target_str):
    """Parses the Smart Input into Global, Junction, or ROI (all 1-based)."""
    target_str = str(target_str).strip()
    if not target_str or target_str.lower() in ['nan', 'none']:
        return "Global", None

    if '-' in target_str:
        try:
            start, end = map(int, target_str.split('-'))
            return "ROI", [start, end]
        except:
            return "Error", None

    if target_str.isdigit():
        return "Junction", int(target_str)

    return "Error", None

def validate_7bp_anchor(primer_start, primer_length, junction, is_reverse=False):
    """Validates that a primer spanning a junction has at least 7bp on each side."""
    if is_reverse:
        r_end = primer_start - primer_length + 1
        return (primer_start >= junction + 7) and (r_end <= junction - 6)
    else:
        f_end = primer_start + primer_length - 1
        return (primer_start <= junction - 7) and (f_end >= junction + 6)

def extract_candidate(raw_res, i, gene_id):
    """Extracts a single primer-pair candidate (index i) from primer3 output."""
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
    """Builds a plotly amplicon track map for a single gene's candidates."""
    fig = go.Figure()

    # Backbone of the gene sequence
    fig.add_shape(type="rect", x0=0, y0=0, x1=seq_len, y1=1,
                  line=dict(color="gray"), fillcolor="lightgray")

    # One track per candidate amplicon
    for cand in candidates:
        rank = cand['Rank']
        y = 1.5 + ((rank - 1) * 0.6)
        x0 = cand['f_start']
        x1 = cand['r_start']
        fig.add_shape(type="rect", x0=x0, y0=y - 0.2, x1=x1, y1=y + 0.2,
                      fillcolor="LimeGreen", opacity=0.8, line=dict(color="green"))
        fig.add_annotation(x=(x0 + x1) / 2, y=y,
                           text=f"Rank {rank} ({cand['Amp_Len']} bp)", showarrow=False)

    # Red dashed line for the exon-exon junction, if any
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

# --- Page Setup ---
st.set_page_config(page_title="SPUD - Batch Mode", page_icon="🧬", layout="wide")

# --- Documentation Sidebar ---
with st.sidebar:
    st.title("📖 SPUD Documentation")
    st.markdown("Welcome to **SPUD** (Specific Primer Universal Designer) - **Batch Edition**.")

    with st.expander("🎯 1. Smart Target Input"):
        st.write("""
        The **Target** field automatically detects your design strategy.
        **All positions are 1-based** (the first base of the sequence is position 1):
        * **Single Number (e.g., 452):** Triggers **Junction Mode**. SPUD ensures the primer spans an exon-exon junction with a strict **7-bp anchor** on each side.
        * **Range (e.g., 300-450):** Triggers **ROI Mode**. Forces the entire amplicon to be designed **inside** this specific sequence range.
        * **Empty:** Triggers **Global Mode**. Searches the entire sequence for the best thermodynamic pairs.
        """)

    with st.expander("📦 2. Batch Processing"):
        st.write("""
        1. **Download Template:** Use the provided CSV to format your data.
        2. **Upload:** SPUD will process every row independently without crashing on errors.
        3. **BLAST Strategy:** Choose Off, Rank 1 only, or All candidates. Fewer BLAST calls = faster and less risk of NCBI blocking.
        4. **Quality Flagging:** If a primer pair's mean Tm deviates by >1.5°C from your **Target Tm**, it will be flagged as an **Outlier**.
        """)

    with st.expander("📊 3. Output Parameters"):
        st.write("""
        * **Rank:** Priority based on penalty score.
        * **Tm_Flag:** Marks primers that deviate from your Target Tm (may not fit a global plate PCR program).
        * **Primer_Fold:** (v / !) Assesses whether the primer itself forms a hairpin that competes with annealing.
        * **Status:** Shows BLAST results (Forward + Reverse) or error messages for failed designs.
        """)

    with st.expander("🗺️ 4. Amplicon Map"):
        st.write("""
        After a run, a visual **track map** is shown per gene: a gray backbone of
        the full sequence, green amplicon bars per candidate (labeled by Rank and length),
        and a red dashed line marking the exon-exon junction (in Junction mode).
        In batch mode, one map is rendered per gene, stacked vertically.
        """)

# --- Main Header ---
st.title("🧬 SPUD: Specific Primer Universal Designer (Batch Edition)")
st.divider()

# --- Top Settings (Global Conditions) ---
with st.expander("⚙️ Lab Conditions & Advanced Thermodynamics (Applied to all)", expanded=True):
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        species_name = st.text_input("Global Species (for BLAST):", "Solanum tuberosum")
        num_returns = st.slider("Max Candidates per Gene", 1, 10, 3, step=1)
    with col2:
        primer_conc = st.number_input("Primer Conc. (nM)", value=250.0, step=10.0)
        mg_conc = st.number_input("Mg2+ Conc. (mM)", value=2.5, step=0.1)
    with col3:
        target_tm = st.slider("Target Tm (°C)", 50.0, 72.0, 60.0, step=0.5)
    with col4:
        min_amp = st.text_input("Min Amplicon Length", value="80")
        max_amp = st.text_input("Max Amplicon Length", value="150")
        blast_mode = st.radio(
            "🔍 BLAST Specificity Check:",
            ["Off", "Rank 1 Only", "All Candidates"],
            index=1,
            help="Off = no BLAST. Rank 1 = check only the top candidate per gene "
                 "(fast, avoids NCBI blocking). All = check every candidate "
                 "(slow, risk of NCBI rate-limiting)."
        )

# --- Mode Selection ---
mode = st.radio("Select Processing Mode:", ["Batch Upload (CSV/Excel)", "Single Gene (Quick Test)"], horizontal=True)

genes_data = []

if mode == "Batch Upload (CSV/Excel)":
    st.info("Upload a file with columns: **Gene_ID**, **Sequence**, **Target**.")

    # Template Download Generator
    template_df = pd.DataFrame({
        "Gene_ID": ["StPOT1", "StGA2ox", "StActin"],
        "Sequence": ["ATGC...", "TTGC...", "CCGG..."],
        "Target": ["452", "300-450", ""]
    })
    csv_buffer = io.BytesIO()
    template_df.to_csv(csv_buffer, index=False)
    st.download_button(label="📥 Download Template", data=csv_buffer.getvalue(), file_name="SPUD_Template.csv", mime="text/csv")

    uploaded_file = st.file_uploader("Upload Batch File", type=['csv', 'xlsx'])
    if uploaded_file:
        try:
            if uploaded_file.name.endswith('.csv'):
                df_input = pd.read_csv(uploaded_file)
            else:
                df_input = pd.read_excel(uploaded_file)

            for index, row in df_input.iterrows():
                genes_data.append({
                    "Gene_ID": str(row.get("Gene_ID", f"Gene_{index}")),
                    "Sequence": str(row.get("Sequence", "")),
                    "Target": str(row.get("Target", ""))
                })
        except Exception as e:
            st.error(f"Error reading file: {e}")

else:
    # Single Gene Input UI
    col_a, col_b = st.columns([2, 1])
    with col_a:
        s_seq = st.text_area("Target Sequence (5' to 3'):", height=150)
    with col_b:
        s_id = st.text_input("Gene ID:", "MyGene")
        s_target = st.text_input("Target (Smart Input):", placeholder="e.g., 452 OR 300-450")
    if s_seq:
        genes_data.append({"Gene_ID": s_id, "Sequence": s_seq, "Target": s_target})

# --- Execution Engine ---
if st.button("🚀 Run SPUD Engine", type="primary") and genes_data:
    all_results = []
    gene_map_data = {}  # gene_id -> {"seq_len", "junction", "candidates"} for the visual map

    global_args = {
        'PRIMER_OPT_SIZE': 20, 'PRIMER_MIN_SIZE': 18, 'PRIMER_MAX_SIZE': 25,
        'PRIMER_OPT_TM': target_tm, 'PRIMER_MIN_TM': target_tm - 5.0, 'PRIMER_MAX_TM': target_tm + 5.0,
        'PRIMER_DNA_CONC': primer_conc, 'PRIMER_SALT_DIVALENT': mg_conc,
        'PRIMER_SALT_MONOVALENT': 50.0, 'PRIMER_DNTP_CONC': 0.8,
        'PRIMER_TM_FORMULA': 1, 'PRIMER_SALT_CORRECTIONS': 1, 'PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT': 1,
        'PRIMER_NUM_RETURN': max(20, num_returns * 5)
    }

    if min_amp.isdigit() and max_amp.isdigit():
        global_args['PRIMER_PRODUCT_SIZE_RANGE'] = [[int(min_amp), int(max_amp)]]

    progress_bar = st.progress(0)
    status_text = st.empty()

    for idx, gene in enumerate(genes_data):
        status_text.text(f"Processing {gene['Gene_ID']} ({idx+1}/{len(genes_data)})...")
        clean_seq = "".join(gene['Sequence'].split()).upper()

        if len(clean_seq) < 50:
            all_results.append({"Gene_ID": gene['Gene_ID'], "Status": "Error: Sequence too short"})
            continue

        target_type, target_val = parse_target(gene['Target'])
        seq_args = {'SEQUENCE_ID': gene['Gene_ID'], 'SEQUENCE_TEMPLATE': clean_seq}
        junction_0based = None

        if target_type == "Junction":
            if target_val < 20 or target_val > len(clean_seq) - 20:
                all_results.append({"Gene_ID": gene['Gene_ID'], "Status": "Error: Junction too close to edge"})
                continue
            # Convert 1-based user input to 0-based position for primer3
            junction_0based = target_val - 1
            seq_args['SEQUENCE_OVERLAP_JUNCTION_LIST'] = [junction_0based]

        elif target_type == "ROI":
            roi_start = max(0, target_val[0] - 1)  # 1-based -> 0-based
            roi_len = target_val[1] - target_val[0] + 1
            if roi_len > 0 and roi_start + roi_len <= len(clean_seq):
                seq_args['SEQUENCE_INCLUDED_REGION'] = [roi_start, roi_len]
            else:
                all_results.append({"Gene_ID": gene['Gene_ID'], "Status": "Error: Invalid ROI range"})
                continue

        elif target_type == "Error":
            all_results.append({"Gene_ID": gene['Gene_ID'], "Status": "Error: Invalid Target format"})
            continue

        try:
            raw_res = primer3.bindings.designPrimers(seq_args, global_args)
            valid_candidates = []

            for i in range(raw_res.get('PRIMER_PAIR_NUM_RETURNED', 0)):
                cand = extract_candidate(raw_res, i, gene['Gene_ID'])

                # Check 7bp rule if in Junction mode (using 0-based junction)
                if target_type == "Junction":
                    if not (validate_7bp_anchor(cand['f_start'], cand['f_len'], junction_0based, False)
                            or validate_7bp_anchor(cand['r_start'], cand['r_len'], junction_0based, True)):
                        continue

                # Secondary Structure: hairpin of the PRIMER ITSELF (not its rev-comp)
                f_hairpin_tm = primer3.calc_hairpin(cand['Forward_Seq']).tm
                r_hairpin_tm = primer3.calc_hairpin(cand['Reverse_Seq']).tm
                sec_struct_flag = "!" if (f_hairpin_tm >= cand['F_Tm'] - 3.0 or r_hairpin_tm >= cand['R_Tm'] - 3.0) else "v"
                if sec_struct_flag == "!":
                    cand['Penalty'] += 50.0

                # Keep coordinates (f_start/r_start) so the map can use them
                valid_candidates.append({
                    "Gene_ID": gene['Gene_ID'],
                    "Forward_Seq": cand['Forward_Seq'],
                    "Reverse_Seq": cand['Reverse_Seq'],
                    "F_Tm": round(cand['F_Tm'], 1),
                    "R_Tm": round(cand['R_Tm'], 1),
                    "Amp_Len": cand['Amp_Len'],
                    "Primer_Fold": sec_struct_flag,
                    "Penalty": cand['Penalty'],
                    "f_start": cand['f_start'],
                    "r_start": cand['r_start'],
                })

            valid_candidates.sort(key=lambda x: x['Penalty'])

            top_candidates = valid_candidates[:num_returns]
            for rank, cand in enumerate(top_candidates):
                cand['Rank'] = rank + 1
                cand['Status'] = "Pending"
                all_results.append(cand)

            # Store map data for this gene (only if it produced candidates)
            if top_candidates:
                gene_map_data[gene['Gene_ID']] = {
                    "seq_len": len(clean_seq),
                    "junction": junction_0based,
                    "candidates": [dict(c) for c in top_candidates],
                }

            if not valid_candidates:
                all_results.append({"Gene_ID": gene['Gene_ID'], "Status": "No candidates found under constraints"})

        except Exception as e:
            all_results.append({"Gene_ID": gene['Gene_ID'], "Status": f"Error: {e}"})

        progress_bar.progress((idx + 1) / len(genes_data))

    # --- Final Audit & Formatting ---
    status_text.text("Finalizing Quality Checks & BLAST...")

    for cand in all_results:
        if cand.get('Status') == "Pending":
            # Outlier is now measured against the user's Target Tm, not the skewed batch mean
            mean_pair_tm = (cand['F_Tm'] + cand['R_Tm']) / 2.0
            cand['Tm_Flag'] = "Outlier" if abs(mean_pair_tm - target_tm) > 1.5 else "OK"

            should_blast = (
                blast_mode == "All Candidates"
                or (blast_mode == "Rank 1 Only" and cand['Rank'] == 1)
            )
            if should_blast:
                status_text.text(f"BLASTing {cand['Gene_ID']} (Rank {cand['Rank']})...")
                cand['Status'] = blast_primer_pair(cand['Forward_Seq'], cand['Reverse_Seq'], species_name)
            else:
                cand['Status'] = "Skipped BLAST"

    status_text.text("Complete!")

    # Organize columns explicitly (hide internal coordinate fields from the table)
    df_results = pd.DataFrame(all_results)
    desired_cols = ['Gene_ID', 'Rank', 'Forward_Seq', 'Reverse_Seq', 'F_Tm', 'R_Tm',
                    'Tm_Flag', 'Amp_Len', 'Primer_Fold', 'Penalty', 'Status']
    final_cols = [c for c in desired_cols if c in df_results.columns]

    # Add any leftover columns (like error messages), but skip internal coords
    for c in df_results.columns:
        if c not in final_cols and c not in ('f_start', 'r_start'):
            final_cols.append(c)

    df_results = df_results[final_cols]

    # Persist results + map data across reruns
    st.session_state['spud_results'] = df_results
    st.session_state['spud_map_data'] = gene_map_data

# --- Results Display (reads from session_state so it survives reruns) ---
if 'spud_results' in st.session_state:
    df_results = st.session_state['spud_results']
    gene_map_data = st.session_state.get('spud_map_data', {})

    st.subheader("📋 Results Table")
    st.dataframe(df_results, use_container_width=True)
    st.download_button(
        label="📥 Download Results (CSV)",
        data=df_results.to_csv(index=False).encode('utf-8'),
        file_name="SPUD_Results.csv",
        mime="text/csv",
        type="primary"
    )

    # --- Visual Amplicon Maps (one per gene, stacked vertically) ---
    if gene_map_data:
        st.divider()
        st.subheader("🗺️ Amplicon Maps")
        for gene_id, data in gene_map_data.items():
            fig = build_gene_map(
                gene_id, data["seq_len"], data["candidates"], data["junction"]
            )
            st.plotly_chart(fig, use_container_width=True)
