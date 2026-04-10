# ChEMBL_target-to-drugs_query_tool
A command-line Python tool that queries the ChEMBL database to discover drugs targeting specific genes/proteins. Returns comprehensive drug information including approval status, mechanism of action, and disease indications [mesh] and [EFO terms].

## 📋 Table of Contents

- [Overview](#overview)
- [Features](#features)
- [What This Tool Does](#what-this-tool-does)
- [Installation](#installation)
- [Usage](#usage)
- [Command Line Arguments](#command-line-arguments)
- [Output Format](#output-format)
- [Examples](#examples)
- [Recommended Targets](#recommended-targets)
- [Logging & Debugging](#logging--debugging)
- [How It Works](#how-it-works)
- [Limitations](#limitations)
- [Requirements](#requirements)
- [License](#license)

## Overview

This tool provides a simple, reproducible way to connect a target gene to potential drugs using ChEMBL's public API. It's designed for bioinformatics researchers, drug discovery scientists, and students who want to quickly identify:

- What drugs target a specific gene/protein
- Whether those drugs are FDA-approved
- How they work (mechanism of action)
- What diseases they treat

**No API key required** - ChEMBL is freely accessible.


## Features

| Feature | Description |
|---------|-------------|
| **7 Key Data Fields** | Drug ID, Name, Approval Year, Max Phase, Mechanism of Action, Action Type, Disease |
| **Disease Filtering** | Filter results by specific disease (MeSH or EFO terms) |
| **Dual Disease Search** | Searches both `mesh_heading` AND `efo_term` for comprehensive coverage |
| **Approval Year Extraction** | Retrieves approval year for FDA-approved drugs |
| **Excel Output** | Clean, formatted Excel file with multiple columns |
| **Verbose Logging** | Optional `--logs` flag for debugging and transparency |
| **Rate-Limited** | Built-in delays to respect ChEMBL API limits |

---

## What This Tool Does

```
Input:  Gene symbol (e.g., "EGFR", "BRCA1", "HER2")
          ↓
Step 1: Find ChEMBL target ID for the gene
          ↓
Step 2: Retrieve all drug mechanisms for that target
          ↓
Step 3: Fetch drug details (name, approval year, disease indications)
          ↓
Step 4: Filter by disease (optional)
          ↓
Output: Excel file with 7 columns of drug information
```

**It answers questions like:**
- "What drugs target EGFR?"
- "Is there an approved drug for BRCA1-related breast cancer?"
- "What is the mechanism of action for Trastuzumab?"
- "Which diseases is Osimertinib approved for?"

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/Sarib13/ChEMBL_target-to-drugs_query_tool.git
cd ChEMBL_target-to-drugs_query_tool
```

### 2. Install required packages

```bash
pip install chembl_webresource_client pandas openpyxl
```

Or using requirements.txt:

```bash
pip install -r requirements.txt
```

---

## Usage

### Basic Usage

```bash
python ChEMBL_target-to-drugs_query_tool.py --gene EGFR
```

### With Disease Filter

```bash
python ChEMBL_target-to-drugs_query_tool.py --gene EGFR --disease "Colorectal Neoplasms"
```

### Custom Output File

```bash
python ChEMBL_target-to-drugs_query_tool.py --gene BRCA1 --output brca1_drugs.xlsx
```

### Limit Number of Drugs

```bash
python ChEMBL_target-to-drugs_query_tool.py --gene HER2 --max-drugs 20
```

### Enable Debug Logging

```bash
python ChEMBL_target-to-drugs_query_tool.py --gene EGFR --logs
```

---

## Command Line Arguments

| Argument | Short | Required | Default | Description |
|----------|-------|----------|---------|-------------|
| `--gene` | `-g` | Yes | None | Target gene symbol (e.g., EGFR, BRCA1, HER2) |
| `--disease` | `-d` | No | None | Disease name for filtering results |
| `--output` | `-o` | No | `chembl_drugs.xlsx` | Output Excel file name |
| `--max-drugs` | - | No | 30 | Maximum number of drugs to fetch |
| `--delay` | - | No | 0.1 | Seconds between API calls |
| `--logs` | - | No | False | Show detailed debug logging |

---

## Output Format

The tool generates an Excel file with the following 7 columns:

| Column | Description | Example |
|--------|-------------|---------|
| **Drug_ID** | ChEMBL unique identifier | `CHEMBL1201827` |
| **Name** | Drug common name | `Panitumumab` |
| **Approval_Year** | FDA approval year or status | `2006` or `Not Approved` |
| **Max_Phase** | Highest development phase (4=approved) | `4` |
| **Mechanism_of_Action** | How the drug affects the target | `Epidermal growth factor receptor antagonist` |
| **Action_Type** | Pharmacological action | `ANTAGONIST`, `INHIBITOR` |
| **Disease** | Disease indications (MeSH + EFO terms) | `Colorectal Neoplasms; Breast Neoplasms` |

---

## Examples

### Example 1: Basic EGFR Query

```bash
python ChEMBL_target-to-drugs_query_tool.py --gene EGFR --max-drugs 10
```

**Output sample:**

| Drug_ID | Name | Approval_Year | Max_Phase | Mechanism_of_Action | Action_Type | Disease |
|---------|------|---------------|-----------|---------------------|-------------|---------|
| CHEMBL1201827 | Panitumumab | 2006 | 4 | EGFR antagonist | ANTAGONIST | Colorectal Neoplasms |
| CHEMBL1201577 | Cetuximab | 2004 | 4 | EGFR antagonist | ANTAGONIST | Colorectal Neoplasms |
| CHEMBL1079742 | Erlotinib | 2004 | 4 | EGFR inhibitor | INHIBITOR | Non-Small Cell Lung Cancer |

### Example 2: BRCA1 with Disease Filter

```bash
python ChEMBL_target-to-drugs_query_tool.py --gene BRCA1 --disease "ovarian cancer" --max-drugs 15
```

**What happens:** Only drugs with "ovarian cancer" in their Disease column are retained.

### Example 3: HER2 with Debug Logging

```bash
python ChEMBL_target-to-drugs_query_tool.py --gene HER2 --logs --output her2_results.xlsx
```

**Logs will show:**
- Which target ID was found
- How many mechanisms retrieved
- Disease fetching progress for each drug
- Filtering results

---

## Recommended Targets

Targets that have direct drug interactions well-documented in ChEMBL:

| Gene | Disease Area |
|------|--------------|
| `EGFR` | Lung, colorectal cancer |
| `HER2` (ERBB2) | Breast cancer |
| `BRCA1` | Ovarian, breast cancer |
| `VEGFA` | Multiple cancers |
| `BCR-ABL` | Leukemia |

### Targets to Avoid

| Gene | Reason |
|------|--------|
| `TP53` | Tumor suppressor - no direct drugs |
| `KRAS` | Historically undruggable - few records |
| `MYC` | Transcription factor - difficult to target |

---

## Logging & Debugging

Use the `--logs` flag to see detailed execution information:

```
[18:01:44] [INFO] Searching for target: EGFR
[18:01:44] [INFO] Found human target: Epidermal growth factor receptor (CHEMBL203)
[18:01:44] [INFO] Retrieving mechanisms for target: CHEMBL203
[18:01:44] [INFO] Found 30 mechanisms

[18:01:44] [INFO] Processing 1/30: CHEMBL1201827
[18:01:44] [DEBUG]   indication_api returned 38 records
[18:01:44] [DEBUG]     Found mesh_heading: Colorectal Neoplasms
[18:01:44] [DEBUG]   ✓ Disease captured: Colorectal Neoplasms
```

This helps you understand:
- Which target ChEMBL selected
- How many drugs were found
- Whether disease information was retrieved

---

## How It Works

### Step-by-Step Technical Flow

1. **Target Search** (`find_target()`)
   - Queries ChEMBL target endpoint with gene symbol
   - Filters for human single proteins
   - Returns ChEMBL target ID (e.g., `CHEMBL203`)

2. **Mechanism Retrieval** (`get_mechanisms_for_target()`)
   - Queries ChEMBL mechanism endpoint
   - Retrieves all drug-target interactions
   - Returns list of mechanism records containing drug IDs, MOA, action type

3. **Drug Details Fetch** (`get_drug_details()`)
   - For each drug ID, queries:
     - Molecule endpoint (name, approval year, max phase)
     - Drug indication endpoint (disease terms)
   - Searches both `mesh_heading` and `efo_term` for diseases

4. **Approval Year Logic** (`get_approval_year()`)
   - If `max_phase == 4`: Drug is approved
   - Extracts `first_approval` year if available
   - Otherwise returns "Approved (Year Unknown)"

5. **Excel Generation**
   - Creates pandas DataFrame
   - Saves to multi-column Excel file

### API Endpoints Used

| Endpoint | Purpose |
|----------|---------|
| `/target` | Find target by gene symbol |
| `/mechanism` | Get drug-target mechanisms |
| `/molecule` | Get drug name and approval data |
| `/drug_indication` | Get disease indications |

---

## Limitations

| Limitation | Explanation |
|------------|-------------|
| **Only direct drug-target interactions** | Does not capture drugs working on upstream/downstream regulators |
| **ChEMBL coverage** | Not all drugs are in ChEMBL; newer drugs may be missing |
| **Approval years** | Not all approved drugs have `first_approval` populated |
| **Disease annotations** | Some drugs lack MeSH or EFO disease terms |
| **Rate limiting** | API delays make large queries (100+ drugs) slower |

---

## Requirements

- Python 3.7+
- Internet connection (for ChEMBL API calls)

### Python Packages

```
chembl_webresource_client>=0.10.0
pandas>=1.3.0
openpyxl>=3.0.0
```

---

## Troubleshooting

### "No mechanisms found for this target"

**Possible causes:**
- The gene has no direct drugs (e.g., TP53, KRAS)
- Gene symbol is misspelled
- Try a different gene

### "No disease information retrieved"

**Solution:**
- Many research compounds lack disease annotations
- Use `--logs` flag to see if API returned records
- Try approved drugs (max_phase=4) for better disease coverage

### Slow performance

**Solution:**
- Reduce `--max-drugs` value (default 30 is reasonable)
- Increase `--delay` if getting rate-limited (rare)

---

## License

MIT License - Free for academic and commercial use.

---

## Author
[Sarib13](https://github.com/Sarib13)
---

---

## Acknowledgments

- [ChEMBL](https://www.ebi.ac.uk/chembl/) for providing the free API
- European Bioinformatics Institute (EBI) for maintaining the database
