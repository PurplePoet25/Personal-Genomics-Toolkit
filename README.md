# Personal Genomics Toolkit

> A desktop app (Tkinter + Polars) for loading, exploring, filtering, and visualizing consumer genotyping raw data (e.g., 23andMe).

---

## ✨ Features

- **Fast loader (Polars backend):** Handles large raw-data text files with comment headers.
- **Clean tabbed UI:** Dedicated tabs for:
  - **Load** (choose File A/B)
  - **Search** (lookup by rsID, chromosome, or position range; live preview)
  - **Filter** (facet filters for zygosity, indels, missing calls; preview & export)
  - **Summary** (key counts, interactive tables)
  - **SNPs per Chromosome** (bar chart)
  - **Genotype Composition** (bar chart)
  - **Chromosome Painting** (two-column layout, long chromosome tracks with legend)
  - **ROH Scan** (runs-of-homozygosity; adjustable window/threshold; result table)
  - **Compare A vs B** (overlap/discordance using shared rsIDs)
  - **Export** (CSV exports from Search/Filter results; HTML Report)
  - **Report** (self-contained HTML with all charts + ROH + metadata)
- **Non-blocking UI:** Heavy work is performed off the main thread; UI stays responsive.
- **Robust plotting:** Matplotlib figures are isolated per tab to prevent cross-updates.
- **23andMe compatible format:** Reads official raw data `.txt` with comment header lines.
- **Sample data included:** Two realistic sample files with overlapping rsIDs:
  - `sample_raw_data_A.txt`
  - `sample_raw_data_B.txt`

> If you are using a previous build, replace your `main.py` with the latest one from this repository.

---

## 🧰 Installation

1. **Python 3.10+ recommended (works with 3.13).**  
2. Create a virtual environment and install requirements:

```bash
python -m venv .venv
# Windows
.venv\Scripts\activate
# macOS/Linux
source .venv/bin/activate

pip install -r requirements.txt
```

3. Run the app:

```bash
python main.py
```

> On first launch, the app will create a working directory (`./output`) for exports and reports.

---

## 📦 Requirements

See [`requirements.txt`](requirements.txt). Core stack:

- **polars** for fast data handling
- **pyarrow** (enables zero-copy conversions where needed)
- **matplotlib** for charts
- **jinja2** for HTML report templating
- **pillow** for saving figure images

---

## 📑 Supported Input Formats

### 1) 23andMe Raw Data (recommended)
- Plain text `.txt`
- **Header lines start with `#`** (metadata, comments)
- **Data lines are TAB-separated** with 4 columns:

```
rsid    chromosome  position  genotype
rs123   1           123456    AG
```

- **Chromosome values:** `1..22`, `X`, `Y`, `MT` (case-insensitive).
- **Genotype values:** two-letter bases (`AA`, `AG`, `CT`, …), indels (`I`, `D` combos), or missing (`--`).

### 2) Generic TSV/CSV
- If you export/convert to CSV/TSV, ensure **columns are named exactly**:  
  `rsid, chromosome, position, genotype`  
- Comment lines (starting with `#`) are ignored automatically.

---

## 🖥️ How to Use

1. **Load Tab**
   - Click **Load File A** and/or **Load File B** and choose a `.txt` or `.csv` file.
   - The app parses headers (lines with `#`) and only reads data rows.

2. **Search Tab**
   - Search by **rsID**, **chromosome**, or **position range** (`start-end`).
   - Results preview instantly; use **Export CSV** to save them.

3. **Filter Tab**
   - Choose zygosity (**Homozygous**, **Heterozygous**, **Indel (I/D)**, **Missing (--)**), chromosome subsets, or genotype masks.
   - **Always shows a table preview** (not just counts). Export results to CSV.

4. **Summary / Charts Tabs**
   - **SNPs per Chromosome** and **Genotype Composition** are **independent figures** (isolated canvases). They don’t mutate each other or other tabs.
   - Chromosome Painting renders in **two columns** with **long track lengths** plus a **legend/key** explaining colors for genotype classes.

5. **ROH Scan Tab**
   - Tune **window size** (e.g., 50–200 SNPs) and **homozygosity threshold** (e.g., ≥95% homozygous within window).
   - Outputs table: chromosome, start–end position, length (bp), #SNPs.
   - If you see “No ROH segments found,” try lowering the threshold or the minimum length/SNPs.

6. **Compare A vs B Tab**
   - Uses **shared rsIDs** to compute overlap, concordance/discordance, and per-chromosome summary.

7. **Export / Report**
   - **Export CSV** from Search/Filter/Compare tabs.
   - **Generate Report (HTML)** creates a **self-contained** file with charts and summary tables without altering any in-app figures.

---

## 🧪 Sample Data

Two synthetic-but-realistic files are included with **overlapping rsIDs** and **varied distributions**:

- `sample_raw_data_A.txt`
- `sample_raw_data_B.txt`

> Use these to test the compare & ROH features and the report generator.

---

## 🛠️ Troubleshooting

- **“the truth value of a DataFrame is ambiguous”**  
  Use `.is_empty()` or explicit length checks; the app does this internally. If you modify code, avoid `if df:` patterns.
- **Freezing during load**  
  Large file parsing runs in a worker thread; if you still see stutters, disable antivirus for the project directory temporarily or move files to a local (non-cloud-synced) folder.
- **PyArrow missing**  
  `pip install pyarrow` (already in `requirements.txt`).
- **Charts collapse to a single bar or disappear**  
  This build **isolates figure state per tab** and **never reuses Matplotlib objects** across tabs or for the report.
- **“No ROH segments found”**  
  Adjust thresholds (lower homozygosity, shorter min length/window). Some modern arrays have sparse coverage in regions.
- **Windows long paths**  
  Keep the project path short (avoid deeply nested OneDrive paths).

---

## 📄 Project Structure

```
.
├── main.py
├── requirements.txt
├── sample_raw_data_A.txt
├── sample_raw_data_B.txt
└── output/
    ├── exports/         # CSV exports
    └── reports/         # HTML reports + embedded images
```

---

## 📷 Screenshots (placeholders)

_Add screenshots/GIFs here to demonstrate each tab:_

- Load & Summary  
- Search & Filter (with preview)  
- SNPs/Chromosome & Genotype Composition  
- Chromosome Painting (two-column)  
- ROH Scan (with segments)  
- Compare A vs B  
- Report (HTML preview)

---

## 🤝 Contributing

PRs and issues welcome. Please include steps to reproduce bugs and attach small data snippets if possible.

---

## ⚠️ Disclaimer

This tool is for **research, educational, and informational** use only. It is **not** a medical device and should not be used for clinical decisions.

---

## 📬 Contact

Open an issue with your OS/Python version and a minimal data excerpt if you need help.
