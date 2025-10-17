# Personal Genomics Kit
#
# Requirements (Python 3.13):
#   py -3.13 -m pip install "polars==1.34.0" pandas>=2 numpy>=1.24 matplotlib>=3.7

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import threading, os, webbrowser, tempfile, shutil, math
import polars as pl
import pandas as pd
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

APP_TITLE = "Personal Genomics Toolkit"
DBSNP_URL = "https://www.ncbi.nlm.nih.gov/snp/{rsid}"
EXPECTED_COLS = ["rsid","chromosome","position","genotype"]
GENO_CLASSES = ["Homozygous", "Heterozygous", "Indel(I/D)", "Missing(--)", "Other"]

HELP_TEXT = """\
Load one or two 23andMe adjacent raw_data.txt files (Build 37 typical) and explore.

Tabs:
• Summary — key counts and quality stats (pure Polars).
• Visualize — per-chromosome bars or genotype composition (respects current filter).
• Chromosome Painting — heterozygosity density, long bars, **two columns**, scrollable.
• Search — rsID full/partial; double-click opens dbSNP; export hits.
• Filter & Export — filter by chromosome / genotype class / rsID contains; preview + export CSV.
• Table — paged (1k rows); optional annotations.csv (rsid,label) shows in “ann”.
• ROH (educational) — simple runs of homozygosity; can use whole data or current filter.
• Compare A vs B — overlap & conflicts across shared rsIDs.

IMPORTANT: Research/education only — not medical use.
"""

# ──────────────────────────── Loader (Polars + IPC cache) ─────────────────────────

def read_23andme_pl(path: str) -> pl.DataFrame:
    """Read 23andMe file with Polars 1.34 (Windows). Skips leading comment lines."""
    comment_rows = 0
    has_real = False
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                comment_rows += 1
                continue
            if line.lower().lstrip().startswith("rsid\tchromosome\tposition\tgenotype"):
                has_real = True
            break

    if has_real:
        df = pl.read_csv(path, separator="\t", has_header=True, skip_rows=comment_rows,
                         infer_schema_length=0, ignore_errors=True, rechunk=True)
    else:
        df = pl.read_csv(path, separator="\t", has_header=False, skip_rows=comment_rows,
                         new_columns=EXPECTED_COLS, infer_schema_length=0, ignore_errors=True, rechunk=True)

    df = df.rename({c: str(c).strip().lower() for c in df.columns})
    for col in EXPECTED_COLS:
        if col not in df.columns:
            raise ValueError(f"Missing expected column '{col}'")

    df = df.with_columns([
        pl.col("rsid").cast(pl.Utf8).str.strip_chars(),
        pl.col("chromosome").cast(pl.Utf8).str.replace("^chr", "", literal=False).str.to_uppercase(),
        pl.col("position").cast(pl.Int64, strict=False),
        pl.col("genotype").cast(pl.Utf8).str.to_uppercase().str.strip_chars(),
    ]).drop_nulls(["rsid","chromosome","position","genotype"])
    return df.rechunk()

def feather_cache_path(path:str)->str:
    return path + ".feather"

def load_with_cache(path:str)->pl.DataFrame:
    fp = feather_cache_path(path)
    try:
        src_mtime = os.path.getmtime(path)
        if os.path.exists(fp) and os.path.getmtime(fp) >= src_mtime:
            return pl.read_ipc(fp, memory_map=True)
    except Exception:
        pass
    df = read_23andme_pl(path)
    try:
        df.write_ipc(fp)
    except Exception:
        pass
    return df

# ───────────────────────────── Utilities / Stats (Polars) ─────────────────────────

def classify_genotype_expr(expr: pl.Expr) -> pl.Expr:
    return (
        pl.when(expr == "--").then(pl.lit("Missing(--)"))
         .when(expr.str.contains("^[ID]{1,2}$")).then(pl.lit("Indel(I/D)"))
         .when((expr.str.len_bytes()==2) & (expr.str.contains("^[ACGT]{2}$")) &
               (expr.str.slice(0,1)==expr.str.slice(1,1))).then(pl.lit("Homozygous"))
         .when((expr.str.len_bytes()==2) & (expr.str.contains("^[ACGT]{2}$"))).then(pl.lit("Heterozygous"))
         .otherwise(pl.lit("Other"))
    )

def chr_sort_key_str(ch:str)->int:
    mapping = {str(i): i for i in range(1, 23)}
    mapping.update({"X": 23, "Y": 24, "MT": 25, "M": 25, "MTDNA": 25})
    return mapping.get(str(ch), 1000)

def summarize_pl(df: pl.DataFrame) -> dict:
    total = df.height
    ch_counts = (
        df.group_by("chromosome").len()
          .with_columns(pl.col("chromosome").map_elements(chr_sort_key_str, return_dtype=pl.Int64).alias("ord"))
          .sort("ord").select(["chromosome","len"])
    )
    ch_dict = {row[0]: int(row[1]) for row in ch_counts.iter_rows()}

    classes_pl = df.with_columns(geno_class=classify_genotype_expr(pl.col("genotype"))).group_by("geno_class").len()
    classes = {row[0]: int(row[1]) for row in classes_pl.iter_rows()}
    geno_counts = {k: int(classes.get(k, 0)) for k in GENO_CLASSES}
    miss_rate = (classes.get("Missing(--)", 0) / total) if total else 0.0
    indel_rate = (classes.get("Indel(I/D)", 0) / total) if total else 0.0

    # GC content via Polars list ops
    gc_stats = (
        df.select(pl.col("genotype"))
          .with_columns(pl.when(pl.col("genotype")!="--").then(pl.col("genotype")).otherwise(None).alias("gt"))
          .drop_nulls("gt")
          .with_columns(pl.col("gt").str.split("").alias("chars"))
          .with_columns(pl.col("chars").list.slice(1, 2).alias("alleles"))
          .select(pl.col("alleles").list.explode().alias("allele"))
          .filter(pl.col("allele").is_in(["A","C","G","T"]))
          .with_columns((pl.col("allele").is_in(["G","C"]).cast(pl.Int64)).alias("is_gc"))
          .select([pl.len().alias("n_all"), pl.col("is_gc").sum().alias("n_gc")])
    )
    if gc_stats.height == 0:
        gc_content = 0.0
    else:
        n_all = int(gc_stats["n_all"][0]); n_gc = int(gc_stats["n_gc"][0])
        gc_content = (n_gc / n_all) if n_all else 0.0

    return {"total_snps": total, "per_chromosome": ch_dict, "geno_counts": geno_counts,
            "missing_rate": miss_rate, "indel_rate": indel_rate, "gc_content": gc_content}

def col_contains_expr(col: str, term: str) -> pl.Expr:
    try:
        return pl.col(col).str.contains(term, literal=False)
    except TypeError:
        return pl.col(col).str.contains(term, literal=True)

def naive_roh_pl(df: pl.DataFrame, min_run=50, max_gap_bp=300000) -> pd.DataFrame:
    tmp = df.with_columns(geno_class=classify_genotype_expr(pl.col("genotype"))) \
            .select(["chromosome","position","geno_class"]).to_pandas()
    out = []
    for ch, sub in tmp.sort_values("position").groupby("chromosome"):
        start = None; prev = None; n = 0
        for _, row in sub.iterrows():
            hom = (row["geno_class"]=="Homozygous")
            if hom:
                if start is None:
                    start = row["position"]; n = 1
                else:
                    if prev is not None and (row["position"]-prev) > max_gap_bp:
                        if n >= min_run: out.append((ch, int(start), int(prev), int(n)))
                        start = row["position"]; n = 1
                    else:
                        n += 1
            else:
                if start is not None and n >= min_run: out.append((ch, int(start), int(prev), int(n)))
                start = None; n = 0
            prev = row["position"]
        if start is not None and n >= min_run: out.append((ch, int(start), int(prev), int(n)))
    res = pd.DataFrame(out, columns=["chromosome","start_bp","end_bp","n_snps"])
    if not res.empty: res["length_bp"] = res["end_bp"] - res["start_bp"]
    return res

# ────────────────────────────────────── GUI ───────────────────────────────────────

class ExplorerApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title(APP_TITLE); self.geometry("1360x940"); self.minsize(1200,820)
        self.dfA: pl.DataFrame|None = None
        self.dfB: pl.DataFrame|None = None
        self.filteredA: pl.DataFrame|None = None
        self.annotations = None  # rsid -> label
        self.page_index = 0; self.page_size = 1000
        self._build_ui()

    def _build_ui(self):
        top = ttk.Frame(self); top.pack(side=tk.TOP, fill=tk.X, padx=8, pady=6)
        ttk.Button(top, text="Open File A", command=lambda:self.open_file("A")).pack(side=tk.LEFT)
        ttk.Button(top, text="Open File B (optional)", command=lambda:self.open_file("B")).pack(side=tk.LEFT, padx=(6,0))
        ttk.Button(top, text="Help", command=lambda:messagebox.showinfo("Help", HELP_TEXT)).pack(side=tk.LEFT, padx=(6,0))
        ttk.Button(top, text="Report (HTML)", command=self.export_report).pack(side=tk.LEFT, padx=(6,0))
        self.status = tk.StringVar(value="Load a 23andMe file to begin…")
        ttk.Label(top, textvariable=self.status, foreground="#4b5563").pack(side=tk.LEFT, padx=12)

        self.nb = ttk.Notebook(self); self.nb.pack(fill=tk.BOTH, expand=True, padx=8, pady=8)

        self.tab_summary = ttk.Frame(self.nb); self.nb.add(self.tab_summary, text="Summary"); self._build_summary()
        self.tab_viz = ttk.Frame(self.nb); self.nb.add(self.tab_viz, text="Visualize"); self._build_visualize()
        self.tab_paint = ttk.Frame(self.nb); self.nb.add(self.tab_paint, text="Chromosome Painting"); self._build_paint()
        self.tab_search = ttk.Frame(self.nb); self.nb.add(self.tab_search, text="Search"); self._build_search()
        self.tab_filter = ttk.Frame(self.nb); self.nb.add(self.tab_filter, text="Filter & Export"); self._build_filter_export()
        self.tab_table = ttk.Frame(self.nb); self.nb.add(self.tab_table, text="Table"); self._build_table()
        self.tab_roh = ttk.Frame(self.nb); self.nb.add(self.tab_roh, text="ROH (educational)"); self._build_roh()
        self.tab_cmp = ttk.Frame(self.nb); self.nb.add(self.tab_cmp, text="Compare A vs B"); self._build_compare()

    # ── Summary
    def _build_summary(self):
        frm = ttk.Frame(self.tab_summary); frm.pack(fill=tk.BOTH, expand=True, padx=12, pady=12)
        self.txt_summary = tk.Text(frm, height=20, wrap="word")
        self.txt_summary.pack(fill=tk.BOTH, expand=True)
        self.txt_summary.insert("1.0", HELP_TEXT); self.txt_summary.configure(state="disabled")

    def refresh_summary(self):
        if self.dfA is None or self.dfA.is_empty(): return
        s = summarize_pl(self.filteredA if self.filteredA is not None else self.dfA)
        lines = [f"Total SNP rows (A): {s['total_snps']:,}\n", "Per-chromosome SNP counts:"]
        for ch,cnt in sorted(s["per_chromosome"].items(), key=lambda kv: chr_sort_key_str(kv[0])):
            lines.append(f"  chr{ch}: {cnt:,}")
        lines.append("\nGenotype class counts:")
        for k in GENO_CLASSES:
            lines.append(f"  {k}: {s['geno_counts'].get(k,0):,}")
        lines.append("")
        lines.append(f"Missing rate: {s['missing_rate']*100:.2f}%")
        lines.append(f"Indel(I/D) rate: {s['indel_rate']*100:.2f}%")
        lines.append(f"GC content (A/C/G/T alleles): {s['gc_content']*100:.2f}%")
        if self.annotations:
            lines.append(f"\nAnnotation set loaded: {len(self.annotations):,} rsIDs flagged")
        self.txt_summary.configure(state="normal"); self.txt_summary.delete("1.0", tk.END)
        self.txt_summary.insert("1.0", "\n".join(lines)); self.txt_summary.configure(state="disabled")

    # ── Visualize
    def _build_visualize(self):
        frm = ttk.Frame(self.tab_viz); frm.pack(fill=tk.BOTH, expand=True, padx=12, pady=12)
        top = ttk.Frame(frm); top.pack(side=tk.TOP, fill=tk.X)
        ttk.Label(top, text="Chart").pack(side=tk.LEFT)
        self.cmb_chart = ttk.Combobox(top, values=["SNPs per chromosome", "Genotype composition"],
                                      state="readonly", width=28)
        self.cmb_chart.set("SNPs per chromosome"); self.cmb_chart.pack(side=tk.LEFT, padx=6)
        ttk.Button(top, text="Render", command=self.render_chart).pack(side=tk.LEFT)
        self.fig = Figure(figsize=(8.8,4.8), dpi=100); self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master=frm); self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True, pady=(8,0))

    def _bar_with_labels(self, ax, labels, values, rotate=False, title="", xlabel="", ylabel=""):
        ax.clear()
        x = np.arange(len(labels))
        ax.bar(x, values)
        ax.set_xticks(x)
        if rotate:
            ax.set_xticklabels(labels, rotation=20, ha="right")
        else:
            ax.set_xticklabels(labels)
        if xlabel: ax.set_xlabel(xlabel)
        if ylabel: ax.set_ylabel(ylabel)
        if title: ax.set_title(title)

    def render_chart(self):
        if self.dfA is None or self.dfA.is_empty(): return
        df = self.filteredA if self.filteredA is not None else self.dfA
        s = summarize_pl(df)
        if self.cmb_chart.get()=="SNPs per chromosome":
            items = sorted(s["per_chromosome"].items(), key=lambda kv: chr_sort_key_str(kv[0]))
            labels = [k for k,_ in items]; vals = [v for _,v in items]
            self._bar_with_labels(self.ax, labels, vals, rotate=True,
                                  title="SNPs per chromosome", xlabel="Chromosome", ylabel="Count")
        else:
            labels = list(GENO_CLASSES); vals = [s["geno_counts"].get(k,0) for k in labels]
            self._bar_with_labels(self.ax, labels, vals, rotate=True,
                                  title="Genotype composition", ylabel="Count")
        self.fig.tight_layout(); self.canvas.draw_idle()

    # ── Chromosome painting (2 columns, scrollable)
    def _build_paint(self):
        outer = ttk.Frame(self.tab_paint); outer.pack(fill=tk.BOTH, expand=True, padx=12, pady=12)
        # Scrollable canvas
        self.paint_canvas = tk.Canvas(outer, background="#f8fafc")
        vbar = ttk.Scrollbar(outer, orient=tk.VERTICAL, command=self.paint_canvas.yview)
        self.paint_canvas.configure(yscrollcommand=vbar.set)
        self.paint_canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        vbar.pack(side=tk.RIGHT, fill=tk.Y)

        # Redraw on size changes
        self.paint_canvas.bind("<Configure>", lambda e: self._schedule_paint_redraw())
        ttk.Label(self.tab_paint, text="Shade = heterozygosity density (0 = dark, 1 = light)").pack(pady=(0,8))

        self._paint_redraw_pending = False

    def _schedule_paint_redraw(self):
        if self._paint_redraw_pending: return
        self._paint_redraw_pending = True
        self.after(150, self._do_paint_redraw)

    def _do_paint_redraw(self):
        self._paint_redraw_pending = False
        self.render_paint()

    def render_paint(self):
        if self.dfA is None or self.dfA.is_empty(): 
            self.paint_canvas.delete("all"); 
            self.paint_canvas.configure(scrollregion=(0,0,0,0))
            return

        df = self.filteredA if self.filteredA is not None else self.dfA

        chs = sorted(df.select(pl.col("chromosome")).unique().to_series().to_list(), key=chr_sort_key_str)
        if not chs:
            self.paint_canvas.delete("all"); 
            self.paint_canvas.configure(scrollregion=(0,0,0,0))
            return

        # Layout
        self.paint_canvas.update_idletasks()
        W = self.paint_canvas.winfo_width() or 1200
        margin = 18
        gutter = 80
        cols = 2
        col_width = (W - (margin*2) - gutter) // cols
        lefts = [margin, margin + col_width + gutter]
        rights = [l + col_width for l in lefts]

        bar_h = 24
        gap = 16
        rows_per_col = math.ceil(len(chs) / cols)
        total_height = 70 + rows_per_col*(bar_h+gap)
        self.paint_canvas.delete("all")

        # Legend (top)
        lg_left = margin + 100
        for i in range(200):
            val = int(255*(i/199))
            color = f"#{val:02x}{val:02x}{val:02x}"
            x0 = lg_left + i*2
            self.paint_canvas.create_rectangle(x0, 18, x0+2, 30, outline="", fill=color)
        self.paint_canvas.create_text(lg_left-8, 24, text="0", anchor="e")
        self.paint_canvas.create_text(lg_left+400+8, 24, text="1", anchor="w")
        self.paint_canvas.create_text(margin, 24, text="Heterozygosity density", anchor="w", font=("Segoe UI", 10, "bold"))

        pandas = df.select(["chromosome","position","genotype"]).to_pandas()
        # Draw per chromosome in 2 columns
        for idx, ch in enumerate(chs):
            col = idx // rows_per_col   # 0 or 1
            row = idx % rows_per_col
            left = lefts[col]; right = rights[col]
            y = 60 + row*(bar_h+gap)

            sub = pandas[pandas["chromosome"]==ch]
            self.paint_canvas.create_text(left-60, y+bar_h/2, text=f"chr{ch}", anchor="w")
            if len(sub) < 5:
                # Draw empty bar
                self.paint_canvas.create_rectangle(left, y, right, y+bar_h, outline="#cbd5e1", fill="#e5e7eb")
                continue

            sub = sub.sort_values("position")
            het = (sub["genotype"].str.len()==2) & (sub["genotype"].str[0]!=sub["genotype"].str[1])
            x0 = sub["position"].iloc[0]; x1 = sub["position"].iloc[-1]; span = max(1, x1-x0)

            bins = 300
            xs = ((sub["position"]-x0)/span).to_numpy()
            het_arr = het.to_numpy().astype(float)
            counts = np.zeros(bins); hsum = np.zeros(bins)
            idxs = np.minimum((xs*bins).astype(int), bins-1)
            for i,b in enumerate(idxs):
                counts[b]+=1; hsum[b]+=het_arr[i]
            dens = np.divide(hsum, np.maximum(counts,1))

            for i,d in enumerate(dens):
                xseg0 = left + (right-left)*i/bins
                xseg1 = left + (right-left)*(i+1)/bins
                val = int(255*(d))  # higher het -> lighter
                color = f"#{val:02x}{val:02x}{val:02x}"
                self.paint_canvas.create_rectangle(xseg0, y, xseg1, y+bar_h, outline="", fill=color)

        # Scroll region
        self.paint_canvas.configure(scrollregion=(0, 0, W, total_height))

    # ── Search
    def _build_search(self):
        frm = ttk.Frame(self.tab_search); frm.pack(fill=tk.BOTH, expand=True, padx=12, pady=12)
        top = ttk.Frame(frm); top.pack(side=tk.TOP, fill=tk.X)
        ttk.Label(top, text="Search rsID (full or contains):").pack(side=tk.LEFT)
        self.ent_search = ttk.Entry(top, width=30); self.ent_search.pack(side=tk.LEFT, padx=6)
        ttk.Button(top, text="Find", command=self.run_search).pack(side=tk.LEFT)
        ttk.Button(top, text="Export hits CSV", command=self.export_search_hits).pack(side=tk.LEFT, padx=8)

        self.search_table = ttk.Treeview(frm, columns=("rsid","chromosome","position","genotype"),
                                         show="headings", height=22)
        for col,w in zip(("rsid","chromosome","position","genotype"), (200,110,120,120)):
            self.search_table.heading(col, text=col); self.search_table.column(col, width=w, anchor=tk.CENTER)
        self.search_table.pack(fill=tk.BOTH, expand=True, pady=8)
        self.search_table.bind("<Double-1>", lambda e:self.open_selected_from(self.search_table))

    def run_search(self):
        if self.dfA is None or self.dfA.is_empty(): return
        term = self.ent_search.get().strip()
        if not term: return
        df = self.dfA.filter(col_contains_expr("rsid", term)).select(EXPECTED_COLS)
        dfp = df.head(1000).to_pandas()
        for i in self.search_table.get_children(): self.search_table.delete(i)
        for _,r in dfp.iterrows():
            self.search_table.insert("", tk.END, values=(r["rsid"], r["chromosome"], int(r["position"]), r["genotype"]))
        self.status.set(f"Search hits: {df.height:,} (showing first {min(len(dfp),1000)})")

    def export_search_hits(self):
        if self.dfA is None or self.dfA.is_empty(): return
        term = self.ent_search.get().strip()
        if not term: return
        df = self.dfA.filter(col_contains_expr("rsid", term)).select(EXPECTED_COLS)
        out = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV","*.csv")], title="Save search hits CSV")
        if not out: return
        df.to_pandas().to_csv(out, index=False)
        messagebox.showinfo("Saved", f"Exported {df.height:,} matching rows.")

    def open_selected_from(self, tree):
        sel = tree.selection()
        if not sel: return
        rsid = tree.item(sel[0], "values")[0]
        webbrowser.open_new(DBSNP_URL.format(rsid=rsid))

    # ── Filter & Export (with preview)
    def _build_filter_export(self):
        frm = ttk.Frame(self.tab_filter); frm.pack(fill=tk.BOTH, expand=True, padx=12, pady=12)
        row = ttk.Frame(frm); row.pack(side=tk.TOP, fill=tk.X, pady=(0,8))
        ttk.Label(row, text="Chromosome").pack(side=tk.LEFT)
        self.flt_chr = ttk.Combobox(row, values=["(any)"], width=10, state="readonly")
        self.flt_chr.set("(any)"); self.flt_chr.pack(side=tk.LEFT, padx=6)
        ttk.Label(row, text="Genotype class").pack(side=tk.LEFT, padx=(10,0))
        self.flt_class = ttk.Combobox(row, values=["(any)"]+GENO_CLASSES, width=16, state="readonly")
        self.flt_class.set("(any)"); self.flt_class.pack(side=tk.LEFT, padx=6)
        ttk.Label(row, text="rsID contains").pack(side=tk.LEFT, padx=(10,0))
        self.flt_rsid = ttk.Entry(row, width=20); self.flt_rsid.pack(side=tk.LEFT, padx=6)
        ttk.Button(row, text="Apply & Preview", command=self.apply_filters_to_state).pack(side=tk.LEFT, padx=8)
        ttk.Button(row, text="Export CSV", command=self.export_filtered).pack(side=tk.LEFT)
        self.flt_info = ttk.Label(frm, text="No filter applied."); self.flt_info.pack(anchor="w", pady=(4,6))
        self.preview_table = ttk.Treeview(frm, columns=("rsid","chromosome","position","genotype"),
                                          show="headings", height=18)
        for col,w in zip(("rsid","chromosome","position","genotype"), (200,110,120,120)):
            self.preview_table.heading(col, text=col); self.preview_table.column(col, width=w, anchor=tk.CENTER)
        self.preview_table.pack(fill=tk.BOTH, expand=True)
        self.preview_table.bind("<Double-1>", lambda e:self.open_selected_from(self.preview_table))

    def apply_filters_to_state(self):
        if self.dfA is None or self.dfA.is_empty(): return
        df = self.dfA
        if self.flt_chr.get() != "(any)":
            df = df.filter(pl.col("chromosome") == self.flt_chr.get())
        if self.flt_class.get() != "(any)":
            clazz = self.flt_class.get()
            df = df.with_columns(tmp=classify_genotype_expr(pl.col("genotype"))).filter(pl.col("tmp")==clazz).drop("tmp")
        term = self.flt_rsid.get().strip()
        if term:
            df = df.filter(col_contains_expr("rsid", term))
        self.filteredA = df
        self.page_index = 0
        self.flt_info.configure(text=f"Filter active: {df.height:,} rows (preview first 1,000 below)")
        for i in self.preview_table.get_children(): self.preview_table.delete(i)
        for r in df.select(EXPECTED_COLS).head(1000).to_pandas().itertuples(index=False):
            self.preview_table.insert("", tk.END, values=(r.rsid, r.chromosome, int(r.position), r.genotype))
        self.refresh_summary(); self.render_chart(); self.render_paint(); self.refresh_table()

    # ── Table (paged)
    def _build_table(self):
        top = ttk.Frame(self.tab_table); top.pack(side=tk.TOP, fill=tk.X, padx=12, pady=8)
        ttk.Button(top, text="Prev", command=lambda:self.change_page(-1)).pack(side=tk.LEFT)
        ttk.Button(top, text="Next", command=lambda:self.change_page(1)).pack(side=tk.LEFT, padx=6)
        self.lbl_page = ttk.Label(top, text="Page 1"); self.lbl_page.pack(side=tk.LEFT, padx=10)
        ttk.Button(top, text="Open dbSNP for selected",
                   command=lambda:self.open_selected_from(self.table)).pack(side=tk.LEFT, padx=6)
        frm = ttk.Frame(self.tab_table); frm.pack(fill=tk.BOTH, expand=True, padx=12, pady=6)
        self.table = ttk.Treeview(frm, columns=("rsid","chromosome","position","genotype","ann"),
                                  show="headings", height=23)
        for col,w in zip(("rsid","chromosome","position","genotype","ann"), (200,110,110,120,240)):
            self.table.heading(col, text=col); self.table.column(col, width=w, anchor=tk.CENTER)
        self.table.pack(fill=tk.BOTH, expand=True)

    def change_page(self, delta:int):
        self.page_index = max(0, self.page_index + delta)
        self.refresh_table()

    def refresh_table(self):
        df = self.filteredA if self.filteredA is not None else self.dfA
        if df is None or df.is_empty(): return
        total = df.height
        start = self.page_index*self.page_size
        if start >= total and total>0:
            self.page_index = max(0, (total-1)//self.page_size); start = self.page_index*self.page_size
        end = min(total, start+self.page_size)
        sub = df.select(["rsid","chromosome","position","genotype"]).slice(start, end-start).to_pandas()
        for i in self.table.get_children(): self.table.delete(i)
        if self.annotations:
            ann = self.annotations
            for _,r in sub.iterrows():
                label = ann.get(r["rsid"], "")
                self.table.insert("", tk.END, values=(r["rsid"], r["chromosome"], int(r["position"]), r["genotype"], label))
        else:
            for _,r in sub.iterrows():
                self.table.insert("", tk.END, values=(r["rsid"], r["chromosome"], int(r["position"]), r["genotype"], ""))
        self.lbl_page.configure(text=f"Page {self.page_index+1} (rows {start+1}-{end} of {total:,})")

    # ── ROH
    def _build_roh(self):
        top = ttk.Frame(self.tab_roh); top.pack(side=tk.TOP, fill=tk.X, padx=12, pady=8)
        ttk.Label(top, text="Min SNPs").pack(side=tk.LEFT)
        self.e_min = ttk.Entry(top, width=6); self.e_min.insert(0,"50"); self.e_min.pack(side=tk.LEFT, padx=6)
        ttk.Label(top, text="Max gap (bp)").pack(side=tk.LEFT)
        self.e_gap = ttk.Entry(top, width=12); self.e_gap.insert(0,"300000"); self.e_gap.pack(side=tk.LEFT, padx=6)
        self.roh_use_filter = tk.BooleanVar(value=False)
        ttk.Checkbutton(top, text="Use current filter", variable=self.roh_use_filter).pack(side=tk.LEFT, padx=10)
        ttk.Button(top, text="Find ROH", command=self.start_roh).pack(side=tk.LEFT, padx=8)
        self.txt_roh = tk.Text(self.tab_roh, height=24, wrap="word")
        self.txt_roh.pack(fill=tk.BOTH, expand=True, padx=12, pady=6)

    def start_roh(self):
        if self.dfA is None or self.dfA.is_empty(): return
        try:
            min_snp = int(self.e_min.get()); max_gap = int(self.e_gap.get())
        except Exception:
            messagebox.showerror("Bad input","Min SNPs and Max gap must be integers."); return
        self.txt_roh.delete("1.0", tk.END); self.txt_roh.insert("1.0","Scanning (threaded)…")

        def work():
            scope = (self.filteredA if (self.filteredA is not None and self.roh_use_filter.get()) else self.dfA)
            res = naive_roh_pl(scope, min_run=min_snp, max_gap_bp=max_gap)
            hom_count = scope.with_columns(gc=classify_genotype_expr(pl.col("genotype"))) \
                             .filter(pl.col("gc")=="Homozygous").height
            return res, hom_count

        def done(res, err):
            # DO NOT touch charts/painting here.
            self.txt_roh.delete("1.0", tk.END)
            if err:
                self.txt_roh.insert("1.0", f"ROH failed: {err}"); return
            df_res, hom_count = res
            if df_res.empty:
                self.txt_roh.insert("1.0",
                    f"No ROH segments found.\n"
                    f"(Homozygous SNPs in scope: {hom_count:,}. "
                    f"Try Min SNPs: 20–30 and/or Max gap: 1,000,000 bp.)"
                ); return
            self.txt_roh.insert("1.0", f"Found {len(df_res)} ROH segments:\n\n")
            for _,r in df_res.sort_values(["chromosome","start_bp"]).iterrows():
                self.txt_roh.insert(tk.END, f"chr{r['chromosome']}: {int(r['start_bp']):,}–{int(r['end_bp']):,} "
                                            f"({int(r['n_snps'])} SNPs; ~{int(r['length_bp']):,} bp)\n")
            # Optional safety redraw of painting (no harm if not visible)
            self.render_paint()
        self.run_bg(work, done)

    # ── Compare
    def _build_compare(self):
        frm = ttk.Frame(self.tab_cmp); frm.pack(fill=tk.BOTH, expand=True, padx=12, pady=12)
        self.txt_cmp = tk.Text(frm, height=24, wrap="word"); self.txt_cmp.pack(fill=tk.BOTH, expand=True)

    def refresh_compare(self):
        self.txt_cmp.delete("1.0", tk.END)
        if self.dfA is None or self.dfB is None or self.dfA.is_empty() or self.dfB.is_empty():
            self.txt_cmp.insert("1.0","Load File A and File B to compare."); return
        A = self.dfA.select(["rsid","genotype","chromosome","position"]).rename({"genotype":"genoA"})
        B = self.dfB.select(["rsid","genotype"]).rename({"genotype":"genoB"})
        joined = A.join(B, on="rsid", how="inner")
        total_shared = joined.height
        if total_shared==0:
            self.txt_cmp.insert("1.0","No overlapping rsIDs found."); return
        conflicts = joined.filter(pl.col("genoA")!=pl.col("genoB"))
        conflict_n = conflicts.height
        self.txt_cmp.insert("1.0", f"Shared rsIDs: {total_shared:,}\n"
                                   f"Conflicting genotypes: {conflict_n:,} ({(conflict_n/total_shared*100):.2f}%)\n\n"
                                   f"Preview of conflicts (up to 50):\n")
        for row in conflicts.head(50).to_pandas().itertuples(index=False):
            self.txt_cmp.insert(tk.END, f"{row.rsid}: A={row.genoA}, B={row.genoB} (chr{row.chromosome}:{row.position})\n")

    # ── Export helpers
    def export_filtered(self):
        if self.dfA is None or self.dfA.is_empty(): return
        df = self.filteredA if self.filteredA is not None else self.dfA
        out = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV","*.csv")], title="Save filtered CSV")
        if not out: return
        df.select(EXPECTED_COLS).to_pandas().to_csv(out, index=False)
        messagebox.showinfo("Saved", f"Exported {df.height:,} rows to:\n{out}")

    def export_report(self):
        if self.dfA is None or self.dfA.is_empty(): return
        df = self.filteredA if self.filteredA is not None else self.dfA

        def work():
            s = summarize_pl(df)
            tmpdir = tempfile.mkdtemp()

            # Off-screen figures with explicit ticks/labels (no UI state touched)
            fig1 = Figure(figsize=(8.8,4.8), dpi=120); ax1 = fig1.add_subplot(111)
            items = sorted(s["per_chromosome"].items(), key=lambda kv: chr_sort_key_str(kv[0]))
            labels = [k for k,_ in items]; vals = [v for _,v in items]
            x = np.arange(len(labels)); ax1.bar(x, vals)
            ax1.set_xticks(x); ax1.set_xticklabels(labels, rotation=20, ha="right")
            ax1.set_title("SNPs per chromosome"); ax1.set_xlabel("Chromosome"); ax1.set_ylabel("Count")
            fig1.tight_layout(); p1 = os.path.join(tmpdir, "chart_chr.png"); fig1.savefig(p1)

            fig2 = Figure(figsize=(8.8,4.8), dpi=120); ax2 = fig2.add_subplot(111)
            labels2 = list(GENO_CLASSES); vals2 = [s["geno_counts"].get(k,0) for k in labels2]
            x2 = np.arange(len(labels2)); ax2.bar(x2, vals2)
            ax2.set_xticks(x2); ax2.set_xticklabels(labels2, rotation=20, ha="right")
            ax2.set_title("Genotype composition"); ax2.set_ylabel("Count")
            fig2.tight_layout(); p2 = os.path.join(tmpdir, "chart_comp.png"); fig2.savefig(p2)

            html = f"""
            <html><head><meta charset='utf-8'><title>23andMe Explorer Report</title>
            <style>body{{font-family:system-ui,Segoe UI,Arial,sans-serif;margin:32px}}h1,h2{{margin:0.2em 0}}</style></head><body>
            <h1>23andMe Explorer Report</h1>
            <p><i>Research/education only — not medical.</i></p>
            <h2>Summary</h2>
            <p>Total SNP rows: {s['total_snps']:,}<br>
            Missing rate: {s['missing_rate']*100:.2f}% &nbsp; Indel rate: {s['indel_rate']*100:.2f}% &nbsp; GC content: {s['gc_content']*100:.2f}%</p>
            <h2>Charts</h2>
            <img src='chart_chr.png' width='880'><br>
            <img src='chart_comp.png' width='880'><br>
            <p>(Chromosome painting is interactive in the app.)</p>
            </body></html>
            """
            return tmpdir, p1, p2, html

        def done(res, err):
            if err:
                messagebox.showerror("Report failed", str(err)); return
            tmpdir, p1, p2, html = res
            out = filedialog.asksaveasfilename(defaultextension=".html", filetypes=[("HTML","*.html")], title="Save HTML report")
            if not out:
                shutil.rmtree(tmpdir, ignore_errors=True)
                return
            out_dir = os.path.dirname(out)
            shutil.copy(p1, os.path.join(out_dir, "chart_chr.png"))
            shutil.copy(p2, os.path.join(out_dir, "chart_comp.png"))
            with open(out, "w", encoding="utf-8") as f: f.write(html)
            shutil.rmtree(tmpdir, ignore_errors=True)
            # Safety: repaint to ensure nothing disappears visually
            self.render_paint()
            messagebox.showinfo("Saved", f"Report saved:\n{out}\n\nImages copied alongside the HTML.")

        self.run_bg(work, done)

    # ── File open & threading
    def open_file(self, which:str):
        path = filedialog.askopenfilename(title=f"Open 23andMe raw_data.txt for {which}",
                                          filetypes=[("23andMe text","*.txt"),("All files","*.*")])
        if not path: return
        self.status.set(f"Loading {which}…")

        def work():
            df = load_with_cache(path)
            ann = None
            if which=="A":
                ann_path = os.path.join(os.path.dirname(path), "annotations.csv")
                if os.path.exists(ann_path):
                    try:
                        ann_df = pd.read_csv(ann_path)
                        if {"rsid","label"}.issubset(ann_df.columns):
                            ann = {str(r["rsid"]): str(r["label"]) for _,r in ann_df.iterrows()}
                    except Exception:
                        ann = None
            return df, ann, path

        def done(res, err):
            if err:
                messagebox.showerror("Load failed", str(err)); self.status.set("Load failed."); return
            df, ann, path2 = res
            if which=="A":
                self.dfA = df; self.annotations = ann
                chs = sorted(df.select(pl.col("chromosome")).unique().to_series().to_list(), key=chr_sort_key_str)
                self.flt_chr.configure(values=["(any)"]+chs); self.flt_chr.set("(any)")
                self.flt_class.set("(any)"); self.flt_rsid.delete(0, tk.END)
                self.filteredA = None; self.page_index=0
                self.refresh_summary(); self.render_chart(); self.render_paint(); self.refresh_table()
            else:
                self.dfB = df; self.refresh_compare()
            self.status.set(f"Loaded {which}: {os.path.basename(path2)} — {df.height:,} rows")

        self.run_bg(work, done)

    def run_bg(self, fn, on_done):
        def wrap():
            try:
                res = fn()
                self.after(0, lambda res=res: on_done(res, None))
            except Exception as exc:
                self.after(0, lambda err=exc: on_done(None, err))
        threading.Thread(target=wrap, daemon=True).start()

if __name__ == "__main__":
    ExplorerApp().mainloop()
