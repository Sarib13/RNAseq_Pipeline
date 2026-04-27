#!/usr/bin/env python3
"""
ChEMBL Drug Information Pipeline
Fetches Drug ID, Name, Approval Year, Max Phase, Mechanism of Action, Action Type, Disease.
Uses direct drug-target interactions only.
"""

import time
import os
import sys
import argparse
import pandas as pd
from chembl_webresource_client.new_client import new_client


RNA_SEQ_GENE_SHEET = "Pathway_Related_DEGs"

NEGATIVE_MODULATORS = ["ANTAGONIST", "INHIBITOR", "BLOCKER", "INVERSE AGONIST", "NEGATIVE ALLOSTERIC MODULATOR", "ALLOSTERIC ANTAGONIST", "ANTISENSE INHIBITOR", "RNAI INHIBITOR", "DEGRADER", "RELEASING AGENT", "GENE EDITING NEGATIVE MODULATOR"]
POSITIVE_MODULATORS = ["ACTIVATOR", "AGONIST", "OPENER", "PARTIAL AGONIST", "POSITIVE ALLOSTERIC MODULATOR", "EXOGENOUS GENE", "EXOGENOUS PROTEIN"]
OTHER_ACTION_TYPES = ["BINDING AGENT", "CHELATING AGENT", "CROSS-LINKING AGENT", "DISRUPTING AGENT", "HYDROLYTIC ENZYME", "METHYLATING AGENT", "MODULATOR", "OTHER", "OXIDATIVE ENZYME", "PROTEOLYTIC ENZYME", "REDUCING AGENT", "SEQUESTERING AGENT", "STABILISER", "SUBSTRATE", "VACCINE ANTIGEN"]

LOGGING_ENABLED = False
LOG_FILE_HANDLE = None
ORIGINAL_STDOUT = None
ORIGINAL_STDERR = None


class TeeStream:
    """Mirror stream output to both terminal and log file."""
    def __init__(self, primary_stream, secondary_stream):
        self.primary_stream = primary_stream
        self.secondary_stream = secondary_stream
    def write(self, message):
        self.primary_stream.write(message)
        self.secondary_stream.write(message)
    def flush(self):
        self.primary_stream.flush()
        self.secondary_stream.flush()


def setup_run_logging(log_file_path=None):
    """Initialize run logging that captures all print output."""
    global LOG_FILE_HANDLE, ORIGINAL_STDOUT, ORIGINAL_STDERR
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    final_path = log_file_path or f"chembl_pipeline_run_{timestamp}.log"
    final_path = os.path.abspath(final_path)
    log_dir = os.path.dirname(final_path)
    if log_dir and not os.path.exists(log_dir):
        os.makedirs(log_dir, exist_ok=True)
    ORIGINAL_STDOUT = sys.stdout
    ORIGINAL_STDERR = sys.stderr
    LOG_FILE_HANDLE = open(final_path, "a", encoding="utf-8", buffering=1)
    sys.stdout = TeeStream(ORIGINAL_STDOUT, LOG_FILE_HANDLE)
    sys.stderr = TeeStream(ORIGINAL_STDERR, LOG_FILE_HANDLE)
    return final_path


def teardown_run_logging():
    """Restore stdout/stderr and close log handle."""
    global LOG_FILE_HANDLE, ORIGINAL_STDOUT, ORIGINAL_STDERR
    if ORIGINAL_STDOUT is not None:
        sys.stdout = ORIGINAL_STDOUT
    if ORIGINAL_STDERR is not None:
        sys.stderr = ORIGINAL_STDERR
    if LOG_FILE_HANDLE is not None:
        LOG_FILE_HANDLE.flush()
        LOG_FILE_HANDLE.close()
    LOG_FILE_HANDLE = None
    ORIGINAL_STDOUT = None
    ORIGINAL_STDERR = None

def log_message(msg, level="INFO"):
    """Write logs to file always, and to console when --logs is enabled."""
    timestamp = time.strftime("%H:%M:%S")
    formatted = f"[{timestamp}] [{level}] {msg}"
    if LOG_FILE_HANDLE is not None:
        LOG_FILE_HANDLE.write(formatted + "\n")
        LOG_FILE_HANDLE.flush()
    if LOGGING_ENABLED:
        stream = ORIGINAL_STDOUT if ORIGINAL_STDOUT is not None else sys.stdout
        stream.write(formatted + "\n")
        stream.flush()


def check_api_connection():
    """Return True if ChEMBL API is reachable, otherwise print a short message."""
    try:
        target_api = getattr(new_client, "target", None)
        if target_api is None:
            print("ERROR: ChEMBL API client is not initialized correctly.")
            return False
        probe = target_api.filter(target_chembl_id="CHEMBL1824")
        _ = next(iter(probe), None)
        return True
    except Exception:
        print("ERROR: ChEMBL API server is down or unreachable. Please try again later.")
        return False

def find_target(gene_symbol):
    """Search for target by gene symbol, return ChEMBL target ID."""
    log_message(f"Searching for target: {gene_symbol}")
    target_api = new_client.target
    results = target_api.search(gene_symbol)
    if not results:
        log_message(f"No target found for: {gene_symbol}", "WARNING")
        return None
    for item in results:
        if item.get("target_type") == "SINGLE PROTEIN" and "Homo sapiens" in item.get("organism", ""):
            chembl_id = item["target_chembl_id"]
            log_message(f"Found human target: {item.get('pref_name', 'Unknown')} ({chembl_id})")
            return chembl_id
    chembl_id = results[0]["target_chembl_id"]
    log_message(f"Using fallback target: {results[0].get('pref_name', 'Unknown')} ({chembl_id})", "WARNING")
    return chembl_id

def get_mechanisms_for_target(target_id):
    """Get all mechanisms (drug-target interactions) for a target."""
    log_message(f"Retrieving mechanisms for target: {target_id}")
    mech_api = new_client.mechanism
    mechanisms = mech_api.filter(target_chembl_id=target_id)
    mech_list = list(mechanisms)
    log_message(f"Found {len(mech_list)} mechanisms")
    return mech_list

def get_approval_year(mol_data):
    """Extract approval year from molecule data."""
    max_phase = float(mol_data.get("max_phase", -1))
    if max_phase < 4:
        return "Not Approved"
    if max_phase == 4:
        first_approval = mol_data.get("first_approval")
        if first_approval:
            year_str = str(first_approval)
            if "-" in year_str:
                year_str = year_str.split("-")[0]
            log_message(f"    Found approval year: {year_str}", "DEBUG")
            return year_str
        usan_year = mol_data.get("usan_year")
        if usan_year:
            log_message(f"    Found USAN year: {usan_year}", "DEBUG")
            return str(usan_year)
        year_of_approval = mol_data.get("year_of_approval")
        if year_of_approval:
            log_message(f"    Found year_of_approval: {year_of_approval}", "DEBUG")
            return str(year_of_approval)
        log_message(f"    No approval year found for approved drug", "WARNING")
        return "Approved (Year Unknown)"
    return "Not Approved"

def _contains_disease_text(source_text, disease_text):
    """Case-insensitive substring match for disease text."""
    if not source_text or not disease_text:
        return False
    return str(disease_text).strip().lower() in str(source_text).strip().lower()


def get_drug_details(molecule_id, api_delay=0.1, disease_filter=None):
    """Fetch name, approval year, max phase, and disease details."""
    time.sleep(api_delay)
    molecule_api = new_client.molecule
    indication_api = new_client.drug_indication
    mol = molecule_api.get(molecule_id)
    if not mol:
        log_message(f"Failed to get molecule: {molecule_id}", "ERROR")
        return None
    name = mol.get("pref_name", "N/A")
    if name == "N/A":
        synonyms = mol.get("molecule_synonyms", [])
        if synonyms:
            name = synonyms[0].get("molecule_synonym", "N/A")
    max_phase = mol.get("max_phase", "N/A")
    approval_year = get_approval_year(mol)
    disease_str = "N/A"
    disease_matched = disease_filter is None

    try:
        indications = indication_api.filter(molecule_chembl_id=molecule_id)
        ind_list = list(indications)
        log_message(f"  indication_api returned {len(ind_list)} records for {molecule_id}", "DEBUG")
        
        if len(ind_list) > 0:
            for ind in ind_list[:15]:
                mesh_heading = ind.get("mesh_heading")
                efo_term = ind.get("efo_term")
                if disease_filter:
                    if mesh_heading and _contains_disease_text(mesh_heading, disease_filter):
                        disease_str = efo_term or mesh_heading or "N/A"
                        disease_matched = True
                        log_message(
                            f"    Disease matched mesh_heading='{mesh_heading}', storing efo_term='{disease_str}'",
                            "DEBUG"
                        )
                        break

                    if efo_term and _contains_disease_text(efo_term, disease_filter):
                        disease_str = mesh_heading or efo_term or "N/A"
                        disease_matched = True
                        log_message(
                            f"    Disease matched efo_term='{efo_term}', storing mesh_heading='{disease_str}'",
                            "DEBUG"
                        )
                        break
                else:
                    if mesh_heading:
                        disease_str = mesh_heading
                    elif efo_term:
                        disease_str = efo_term
                    elif ind.get("disease_mesh_name"):
                        disease_str = ind["disease_mesh_name"]
                    elif ind.get("disease_efo_term"):
                        disease_str = ind["disease_efo_term"]
            
            if disease_filter and disease_str != "N/A":
                log_message(f"  ✓ Disease matched for {name}: {disease_str}", "DEBUG")
            elif not disease_filter and disease_str != "N/A":
                log_message(f"  ✓ Disease captured for {name}: {disease_str}", "DEBUG")
        
    except Exception as e:
        log_message(f"Error fetching disease for {molecule_id}: {str(e)}", "ERROR")

    if disease_filter and not disease_matched:
        return None
    return {
        "Drug_ID": molecule_id,
        "Name": name,
        "Approval_Year": approval_year,
        "Max_Phase": max_phase,
        "Disease": disease_str
    }


def _normalize_column_name(col_name):
    """Normalize a column name for robust matching."""
    return "".join(ch for ch in str(col_name).lower() if ch.isalnum())


def normalize_regulation_value(raw_value):
    """Convert regulation text into a stable up/down category."""
    if pd.isna(raw_value):
        return "Unknown"

    value = str(raw_value).strip().lower()
    if not value:
        return "Unknown"

    if value.startswith("up") or "upreg" in value or "overexp" in value or "induc" in value:
        return "Upregulated"
    if value.startswith("down") or "downreg" in value or "underexp" in value or "suppress" in value:
        return "Downregulated"

    return str(raw_value).strip()


def get_relevant_action_types(condition):
    """Return allowed ChEMBL action types for a target regulation state."""
    if condition == "upregulated":
        return NEGATIVE_MODULATORS
    if condition == "downregulated":
        return POSITIVE_MODULATORS
    return OTHER_ACTION_TYPES


def _normalize_condition(regulation):
    """Normalize condition into upregulated/downregulated/other."""
    regulation_norm = str(regulation or "").strip().lower()

    if regulation_norm.startswith("up") or "upreg" in regulation_norm:
        return "upregulated"
    if regulation_norm.startswith("down") or "downreg" in regulation_norm:
        return "downregulated"
    return "other"


def _normalize_action_type(action_type):
    """Normalize ChEMBL action type text for robust matching."""
    return " ".join(str(action_type or "").upper().split())


def action_type_allowed_for_regulation(regulation, action_type):
    """Filter by allowed action types given regulation status."""
    allowed_types = get_relevant_action_types(_normalize_condition(regulation))
    normalized_action_type = _normalize_action_type(action_type)

    if not normalized_action_type:
        return False

    # Match exact action type and common compound values separated by delimiters.
    if normalized_action_type in allowed_types:
        return True

    for delimiter in [";", "|", ","]:
        if delimiter in normalized_action_type:
            parts = [part.strip() for part in normalized_action_type.split(delimiter)]
            if any(part in allowed_types for part in parts):
                return True

    return False


def is_phase_4_approved(drug_info):
    """Return True when a drug is approved (phase 4)."""
    try:
        return float(drug_info.get("Max_Phase", -1)) == 4.0
    except (TypeError, ValueError):
        return False


def _sanitize_sheet_name(name):
    """Return an Excel-compatible sheet name (max 31 chars, no invalid chars)."""
    invalid_chars = set('[]:*?/\\')
    cleaned = "".join(ch for ch in str(name) if ch not in invalid_chars).strip()
    if not cleaned:
        cleaned = "Target"
    return cleaned[:31]


def build_unique_sheet_name(base_name, used_names):
    """Generate a unique Excel sheet name within 31 chars."""
    base = _sanitize_sheet_name(base_name)
    if base not in used_names:
        used_names.add(base)
        return base

    counter = 2
    while True:
        suffix = f"_{counter}"
        candidate = _sanitize_sheet_name(base[: 31 - len(suffix)] + suffix)
        if candidate not in used_names:
            used_names.add(candidate)
            return candidate
        counter += 1


def load_genes_from_excel(excel_path, gene_column=None, sheet_name="0"):
    """Load gene symbols and regulation values from an RNA-seq Excel sheet."""
    try:
        df = pd.read_excel(excel_path, sheet_name=sheet_name)
    except ValueError as exc:
        workbook = pd.ExcelFile(excel_path)
        raise ValueError(
            f"Sheet '{sheet_name}' not found. Available sheets: {workbook.sheet_names}"
        ) from exc

    if df.empty:
        raise ValueError("Excel sheet is empty.")

    regulation_column = None
    normalized_to_original = {
        _normalize_column_name(col): col for col in df.columns
    }
    if "regulation" in normalized_to_original:
        regulation_column = normalized_to_original["regulation"]
    else:
        raise ValueError(
            "Could not find a REGULATION column in the Excel sheet. "
            f"Available columns: {list(df.columns)}"
        )

    if gene_column:
        if gene_column not in df.columns:
            raise ValueError(
                f"Column '{gene_column}' not found. Available columns: {list(df.columns)}"
            )
        selected_column = gene_column
    else:
        candidates = [
            "gene", "genesymbol", "symbol", "hgncsymbol", "genename", "externalgenename"
        ]
        selected_column = None
        for candidate in candidates:
            if candidate in normalized_to_original:
                selected_column = normalized_to_original[candidate]
                break

        if selected_column is None:
            raise ValueError(
                "Could not auto-detect gene column. Please provide --gene-column. "
                f"Available columns: {list(df.columns)}"
            )

    gene_records = []
    seen = set()
    for _, row in df[[selected_column, regulation_column]].dropna(subset=[selected_column]).iterrows():
        gene = str(row[selected_column]).strip()
        if not gene:
            continue
        if gene.lower() in {"na", "none", "nan"}:
            continue

        gene = gene.upper()
        regulation = normalize_regulation_value(row[regulation_column])
        record_key = (gene, regulation)
        if record_key not in seen:
            seen.add(record_key)
            gene_records.append({"gene": gene, "regulation": regulation})

    if not gene_records:
        raise ValueError(f"No valid genes found in column '{selected_column}'.")
    return gene_records, selected_column, regulation_column

def main():
    global LOGGING_ENABLED
    parser = argparse.ArgumentParser(
        description="Fetch drug information from ChEMBL",
        epilog=(
            "Examples: python chembl_drugs.py --gene EGFR TP53\n"
            "          python chembl_drugs.py --gene-excel rna_seq_output.xlsx --gene-column Gene"
        )
    )
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("--gene", "-g", nargs="+", help="One or more target gene symbols (e.g., EGFR TP53)")
    input_group.add_argument("--gene-excel", type=str, help="Path to RNA-seq output Excel containing gene symbols")
    parser.add_argument("--gene-column", type=str, default=None, help="Gene column name in --gene-excel (optional; auto-detect if omitted)")
    parser.add_argument("--disease", "-d", type=str, default=None, help="Disease name for filtering (optional)")
    parser.add_argument("--output", "-o", type=str, default="chembl_drugs.xlsx", help="Output Excel file name")
    parser.add_argument("--delay", type=float, default=0.1, help="Seconds between API calls (default: 0.1)")
    parser.add_argument("--logs", action="store_true", help="Show detailed debug logging output")
    parser.add_argument("--log-file", type=str, default=None, help="Path to run log file (default: chembl_pipeline_run_YYYYMMDD_HHMMSS.log)")
    args = parser.parse_args()
    LOGGING_ENABLED = args.logs
    log_file_path = setup_run_logging(args.log_file)

    try:
        print(f"Run log file: {log_file_path}")

        if not check_api_connection():
            return

        if args.gene:
            gene_records = []
            seen = set()
            for gene in args.gene:
                cleaned = gene.strip().upper()
                if cleaned and cleaned not in seen:
                    seen.add(cleaned)
                    gene_records.append({"gene": cleaned, "regulation": "Unknown"})
            source_info = "CLI --gene"
            selected_sheet = None
            print(f"[STEP] Using CLI gene input... done ({len(gene_records)} unique genes)")
        else:
            try:
                print(f"[STEP] Reading Excel gene input from: {args.gene_excel}")
                gene_records, detected_column, regulation_column = load_genes_from_excel(
                    args.gene_excel,
                    gene_column=args.gene_column,
                    sheet_name=RNA_SEQ_GENE_SHEET
                )
                source_info = (
                    f"Excel {args.gene_excel} "
                    f"(sheet: {RNA_SEQ_GENE_SHEET}, gene column: {detected_column}, regulation column: {regulation_column})"
                )
                selected_sheet = RNA_SEQ_GENE_SHEET
                print(
                    "[STEP] Excel read complete: "
                    f"sheet='{RNA_SEQ_GENE_SHEET}', gene column='{detected_column}', "
                    f"regulation column='{regulation_column}', genes detected={len(gene_records)}"
                )
            except Exception as exc:
                print(f"ERROR: Failed to load genes from Excel: {exc}")
                return

        genes = gene_records

        print("=" * 70)
        print("ChEMBL Drug Pipeline")
        print("=" * 70)
        print(f"Gene source: {source_info}")
        if selected_sheet is not None:
            print(f"Gene sheet: {selected_sheet}")
        print(f"Gene count : {len(genes)}")
        print(f"Disease    : {args.disease if args.disease else 'No filter'}")
        print("Columns    : Target_Gene, Regulation, Target_ChEMBL_ID, Drug_ID, Name, Approval_Year, Max_Phase, Mechanism_of_Action, Action_Type, Disease")
        if args.logs:
            print("Logging    : Enabled")
        else:
            print("Logging    : Console debug disabled (file logging always enabled)")
        print("=" * 70)
        start_time = time.time()
        results = []
        results_by_target = {}
        genes_with_targets = 0

        for gene_idx, gene_record in enumerate(genes, start=1):
            gene = gene_record["gene"]
            regulation = gene_record.get("regulation", "Unknown")
            print(f"\n[{gene_idx}/{len(genes)}] Processing gene: {gene} ({regulation})")
            print("  [STEP] Resolving ChEMBL target...")
            target_id = find_target(gene)
            if not target_id:
                print(f"  WARNING: No target found for {gene}. Skipping.")
                continue
            print(f"  [DONE] Target resolved: {target_id}")
            genes_with_targets += 1
            print(f"  [STEP] Fetching mechanisms for {target_id}...")
            mechanisms = get_mechanisms_for_target(target_id)
            if not mechanisms:
                print(f"  WARNING: No mechanisms found for {gene} ({target_id}).")
                continue
            print(f"  [DONE] Mechanisms fetched: {len(mechanisms)}")
            print("  [STEP] Evaluating mechanism candidates and fetching drug details...")
            target_candidates = []
            for idx, mech in enumerate(mechanisms):
                mol_id = mech.get("molecule_chembl_id")
                if not mol_id:
                    continue

                log_message(f"Processing {idx+1}/{len(mechanisms)}: {mol_id}", "INFO")
                moa = mech.get("mechanism_of_action", "N/A")
                action_type = mech.get("action_type", "N/A")
                if not action_type_allowed_for_regulation(regulation, action_type):
                    continue
                drug_info = get_drug_details(mol_id, api_delay=args.delay, disease_filter=args.disease)

                if drug_info:
                    drug_info["Target_Gene"] = gene
                    drug_info["Regulation"] = regulation
                    drug_info["Target_ChEMBL_ID"] = target_id
                    drug_info["Mechanism_of_Action"] = moa
                    drug_info["Action_Type"] = action_type
                    target_candidates.append(drug_info)

            print(f"  [DONE] Candidate drugs passing filters: {len(target_candidates)}")
            approved_count = sum(1 for drug in target_candidates if is_phase_4_approved(drug))
            if approved_count > 3:
                target_candidates = [drug for drug in target_candidates if is_phase_4_approved(drug)]
                print(f"  Phase 4 gate applied: {approved_count} approved drugs found, keeping approved only.")
            else:
                print(f"  Phase 4 gate not applied: {approved_count} approved drugs found, keeping all phases.")

            if target_candidates:
                target_key = (gene, regulation, target_id)
                results_by_target[target_key] = target_candidates
                results.extend(target_candidates)
        elapsed_time = time.time() - start_time
        if results:
            columns = [
                "Target_Gene", "Regulation", "Target_ChEMBL_ID", "Drug_ID", "Name",
                "Approval_Year", "Max_Phase", "Mechanism_of_Action", "Action_Type", "Disease"
            ]
            print(f"\n[STEP] Writing results to Excel: {args.output}")
            used_sheet_names = set()
            with pd.ExcelWriter(args.output) as writer:
                overview_df = pd.DataFrame(results)[columns]
                overview_df = overview_df.drop_duplicates(
                    subset=["Target_Gene", "Regulation", "Target_ChEMBL_ID", "Drug_ID", "Action_Type"]
                )
                overview_df.to_excel(writer, sheet_name="All_Selected_Drugs", index=False)

                for (gene, regulation, target_id), target_rows in results_by_target.items():
                    target_df = pd.DataFrame(target_rows)[columns]
                    target_df = target_df.drop_duplicates(
                        subset=["Target_Gene", "Regulation", "Target_ChEMBL_ID", "Drug_ID", "Action_Type"]
                    )

                    base_sheet = f"{gene}_{target_id}"
                    sheet_name = build_unique_sheet_name(base_sheet, used_sheet_names)
                    target_df.to_excel(writer, sheet_name=sheet_name, index=False)
            print(f"[DONE] Excel export complete: {len(results)} rows, {len(results_by_target)} target sheets + 1 overview")
            print("\n" + "=" * 70)
            print(f"SUCCESS! Completed in {elapsed_time:.1f} seconds")
            print(f"   Saved {len(results)} drugs to {args.output}")
            print("=" * 70)
            approved_with_year = sum(1 for d in results if d["Approval_Year"] and d["Approval_Year"].isdigit())
            approved_no_year = sum(1 for d in results if d["Approval_Year"] == "Approved (Year Unknown)")
            not_approved = sum(1 for d in results if d["Approval_Year"] == "Not Approved")
            print(f"\nSummary:")
            print(f"  Genes requested: {len(genes)}")
            print(f"  Genes with target match: {genes_with_targets}")
            print(f"  Total drugs: {len(results)}")
            print(f"  Target sheets created: {len(results_by_target)} (+1 overview sheet)")
            print(f"  Approved with year: {approved_with_year}")
            print(f"  Approved (year unknown): {approved_no_year}")
            print(f"  Not approved / in trials: {not_approved}")
            print(f"  Drugs with disease info: {sum(1 for d in results if d['Disease'] != 'N/A')}/{len(results)}")
        else:
            print("\nNo drugs found matching criteria.")
    finally:
        teardown_run_logging()

if __name__ == "__main__":
    main()