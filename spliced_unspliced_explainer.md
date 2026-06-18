# From FASTQ to Spliced/Unspliced: An End-to-End Example

**Sample:** GSM8253792 — CD1, Control, 0 weeks recovery  
**Organism:** *Mus musculus* (GRCm39)  
**Protocol:** 10x Genomics Chromium v3 (single-nucleus RNA-seq)  
**Pipeline:** kallisto|bustools (kb-python) with lamanno chemistry  

---

## 1. What is in the raw FASTQ?

When you download this sample from SRA with `fasterq-dump`, you get two FASTQ files per run:

```
SRR_ID_2.fastq   # Read 1 — 28 bp
SRR_ID_3.fastq   # Read 2 — 91 bp
```

A single sequenced fragment looks like this across the two files:

**R1 (28 bp) — cell identity, NOT genome sequence:**
```
@SRR_ID.1 1/1
CTGTGGGAGGTCACCCGATTACAGTCGC
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
```
- Bases 1–16: **cell barcode** → `CTGTGGGAGGTCACCC`
- Bases 17–28: **UMI** → `GATTACAGTCGC` (random 12-mer; 4¹² ≈ 16.7M possible sequences)

> Note: a poly-T UMI like `TTTTTTTTTTTT` is an artifact of internal priming — the
> poly-dT capture primer on the gel bead mis-anneals to an internal A-rich region of
> the mRNA rather than the 3' poly-A tail, so the bead's own poly-T sequence gets
> sequenced in the UMI slot. These are filtered by kb-python and cellranger before
> molecules are counted.

**R2 (91 bp) — the actual cDNA read that maps to the genome:**
```
@SRR_ID.1 1/2
GCAGCAGTGACTGGAGAAGCCATCAAGCAGCTGCAGAATGAGCTGGAGAAGCAGAAGAAGCAGCAGCAGCAGCAGCAGCAG
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
```

R1 tells you **which cell** and **which molecule** (barcode + UMI).  
R2 tells you **what was transcribed** — this is what gets aligned to the genome.

Together, a barcode+UMI pair defines one unique mRNA **molecule**. Multiple reads can share the same barcode+UMI — those reads form the **supporting read set** for that molecule.

---

## 2. The Reference: How Exons and Introns Are Defined

Before any alignment, you need a genome assembly and gene annotation:

```
Mus_musculus.GRCm39.dna.primary_assembly.fa   # genome sequence
Mus_musculus.GRCm39.110.gtf                   # gene annotation
```

The GTF file defines **exon coordinates** for every transcript. Introns are never listed explicitly — they are the gaps between consecutive exons of the same transcript.

### How a read's coordinates are classified against the GTF map

Once the aligner returns `[read_start, read_end]` on a chromosome, the tool works through this decision tree against the GTF exon intervals:

```
Given: read maps to [read_start, read_end] on chrX
       GTF defines exon intervals [exon_start, exon_end] for each transcript

Is [read_start, read_end] entirely            YES ──→ EXONIC
within one exon [exon_start, exon_end]?               → SPLICED vote
                    │
                    NO
                    ↓
Does the read overlap the gap between         YES ──→ spans exon-intron
two consecutive exons of the same                     boundary
transcript? (read_start < exon_end            → UNSPLICED vote
AND read_end > next_exon_start)
                    │
                    NO
                    ↓
Is the read entirely between two exons        YES ──→ purely intronic
(read_start > exon_end                                (no exon overlap)
AND read_end < next_exon_start)?              → UNSPLICED vote
                    │
                    NO
                    ↓
Does the CIGAR string contain an N (skip)     YES ──→ exon-exon junction
whose length matches the annotated                    (intron already removed)
intron length in the GTF?                     → SPLICED vote
                    │
                    NO
                    ↓
                AMBIGUOUS → discard
```

After collecting all votes across every read in a molecule's supporting read set:

```
All votes → SPLICED                   → molecule counted in spliced.mtx
All votes → UNSPLICED (for all 
  compatible transcript models)       → molecule counted in unspliced.mtx
Mixed votes                           → AMBIGUOUS → discarded
```

---

### Real example: *Becn2* (ENSMUSG00000104158.3) — our spliced-only gene

```
GTF entries for Becn2 on chromosome 1, + strand:

Exon 1:  chr1  175,747,895 – 175,749,254   (1,360 bp)
                          ↕↕↕
Intron 1: chr1  175,749,255 – 175,749,362   (108 bp)  ← inferred gap
                          ↕↕↕
Exon 2:  chr1  175,749,363 – 175,749,791   (429 bp)
```

### Real example: *9330185C12Rik* (ENSMUSG00000097648.2) — our unspliced-only gene

This is a lncRNA on the **minus strand** with 3 exons:

```
Genomic coordinates (left = lower position on chr1):

[Exon 3]              [Exon 2]          [Exon 1]
113,819,115-113,819,778  113,837,065-113,837,093  113,891,549-113,891,864
          └──────── Intron 2 ────────┘└──────────────── Intron 1 (54,455 bp) ───────────────┘
```

Because it is on the minus strand, transcription runs **right to left** on the chromosome. Reads mapping to intron 1 (chr1:113,837,094–113,891,548) support unspliced status.

### Real example: *Pcmtd1* (ENSMUSG00000051285.18) — our both-spliced-and-unspliced gene

Transcript Pcmtd1-201 has 6 exons on chromosome 1, + strand:

```
Exon 1: 7,159,144 – 7,159,440   (297 bp)
Intron 1: 7,159,441 – 7,190,417 (30,977 bp)  ← reads mapping here = unspliced
Exon 2: 7,190,418 – 7,190,839   (422 bp)
Intron 2: 7,190,840 – 7,217,860
Exon 3: 7,217,861 – 7,217,963
... (3 more exons)
Exon 6: 7,239,739 – 7,243,852
```

---

## 3. Building the kb-python Index

kb-python (kallisto|bustools) avoids traditional genome alignment entirely. Instead, it performs **pseudo-alignment** — matching reads to a k-mer index built from transcript sequences. For spliced/unspliced quantification, it builds **two separate indices**:

```bash
kb ref \
  -i index.idx \
  -g t2g.txt \
  -f1 cdna.fa \
  -f2 intron.fa \
  --workflow lamanno \
  Mus_musculus.GRCm39.dna.primary_assembly.fa \
  Mus_musculus.GRCm39.110.gtf
```

This produces:

| File | Contents | Purpose |
|------|----------|---------|
| `cdna.fa` | Spliced transcript sequences (exons only, joined) | Detecting mature mRNA |
| `intron.fa` | Pre-mRNA sequences (full gene body incl. introns) | Detecting nascent mRNA |
| `cdna_t2c.txt` | Transcript → gene for spliced | Spliced counting |
| `intron_t2c.txt` | Transcript-I → gene for intronic | Unspliced counting |
| `index.idx` | Combined k-mer index (both cdna + intron) | Pseudo-alignment target |
| `t2g.txt` | Full transcript-to-gene mapping with coordinates | Gene-level summarisation |

From your actual files:

```
# cdna_t2c.txt (one entry per spliced transcript)
ENSMUST00000180837.2          → 9330185C12Rik

# intron_t2c.txt (suffix -I.1 = intronic capture region)
ENSMUST00000180837.2-I.1      → 9330185C12Rik
```

The `-I.1` suffix marks a capture region that overlaps the intronic sequence. If a read pseudo-aligns to this entry, it is voting for the unspliced count.

---

## 4. Running kb count

```bash
kb count \
  -i sra_runs/index.idx \
  -g sra_runs/t2g.txt \
  -x 10xv3 \
  --lamanno \
  -c1 sra_runs/cdna_t2c.txt \
  -c2 sra_runs/intron_t2c.txt \
  -o sra_runs/kb_output_GSM8253792/ \
  SRR_ID_2.fastq SRR_ID_3.fastq
```

What happens internally:

1. **Barcode + UMI extraction:** R1 is parsed — first 16 bp = cell barcode, next 12 bp = UMI. The barcode is corrected against the 10x v3 whitelist.
2. **Pseudo-alignment:** R2 (91 bp) is broken into k-mers (k=31) and matched against the combined index.
3. **BUS file:** Every (barcode, UMI, equivalence class) tuple is written to a `.bus` file.
4. **Sorting and correction:** `bustools sort` + `bustools correct` refine barcode assignments.
5. **Counting:** `bustools count` using `-c1` (cdna) and `-c2` (intron) capture lists separates molecules into spliced vs unspliced.

Output:
```
kb_output_GSM8253792/counts_unfiltered/
├── spliced.barcodes.txt   # 816,643 barcodes with ≥1 spliced molecule
├── spliced.genes.txt      # 56,941 genes
├── spliced.mtx            # 13,411,965 nonzero (barcode × gene) entries
├── unspliced.barcodes.txt # 712,804 barcodes with ≥1 unspliced molecule
├── unspliced.genes.txt    # 56,941 genes
└── unspliced.mtx          # 11,422,859 nonzero entries
```

The ratio of 11.4M unspliced / 13.4M spliced ≈ **0.85** — consistent with the La Manno et al. observation that 15–25% of 10x reads carry intronic sequence.

---

## 5. How a Single Read Gets Classified

This is where the La Manno et al. rules are applied in practice.

### The classification logic

For every unique molecule (barcode + UMI), kb collects all pseudo-alignment hits across both the cdna and intron indices. The rule is:

```
All hits → cdna only?          → SPLICED   → spliced.mtx
All hits → intron for all 
  compatible transcripts?      → UNSPLICED → unspliced.mtx
Mixed (some cdna, some intron)?→ AMBIGUOUS → discarded
```

---

## 6. Concrete Read-Level Examples

All examples below are from **cell barcode `CTGTGGGAGGTCACCC`** — the highest-UMI cell in this sample (45,354 spliced + 41,369 unspliced molecules, 7,579 genes detected).

---

### Example A: *Becn2* — SPLICED (Rule 1)

**Observed:** 1 spliced UMI, 0 unspliced UMIs

```
Gene structure (chr1, + strand):
[====== Exon 1 ======]          [= Exon 2 =]
175,747,895    175,749,254  175,749,363   175,749,791
                           ↑
                      108 bp intron (175,749,255–175,749,362)

Read R2 maps to:
[=====READ (91bp)=====]
175,747,950       175,748,040

→ Entirely within Exon 1 coordinates.
→ No bases fall in the 108 bp intron or Exon 2.
→ Pseudo-aligns to cdna index ONLY.
→ Annotated: SPLICED → counted in spliced.mtx
```

**Biological interpretation:** This is a **mature, fully processed mRNA** for Becn2. The intron has already been spliced out. The cell has this gene's protein being translated right now.

---

### Example B: *9330185C12Rik* — UNSPLICED (Rule 2)

**Observed:** 0 spliced UMIs, 18 unspliced UMIs

```
Gene structure (chr1, − strand, transcription flows right→left):

[Exon 3]          [Exon 2]     [====== Exon 1 ======]
113,819,115       113,837,065  113,891,549    113,891,864
         └─── Intron 2 ───┘└──────── Intron 1 (54,455 bp) ─────────┘

Scenario i — read spans exon–intron boundary:
                              [=====READ (91bp)=====]
                             113,891,500         113,891,590
→ Read starts in Exon 1 (≥113,891,549) but 50 bp of it falls into 
  Intron 1 (≤113,891,548). Spans the exon-intron boundary.
→ Matches intron index (ENSMUST00000180837.2-I.1), NOT cdna index.

Scenario ii — read maps entirely in intron:
              [=====READ (91bp)=====]
             113,850,000         113,850,090
→ Entirely within intron 1 (113,837,094–113,891,548).
→ This arises from secondary priming at intronic polyT sequences (common 
  in 10x Chromium, as described in La Manno et al. Fig. 1, Extended Data).
→ Matches intron index only.

All 18 molecules fall into scenario i or ii.
→ For the single compatible transcript (9330185C12Rik-201), every molecule 
  has at least one intronic read.
→ Annotated: UNSPLICED → counted in unspliced.mtx (count = 18)
```

**Biological interpretation:** This is **nascent pre-mRNA** — the gene is being actively transcribed but the introns have not been removed yet. High intronic read counts indicate active transcription.

---

### Example C: *Pcmtd1* — BOTH (Rule 1 and Rule 2 for different molecules)

**Observed:** 13 spliced UMIs, 8 unspliced UMIs

These are **21 distinct molecules** (different UMIs) for the same gene in the same cell. Each molecule is classified independently.

```
The 13 SPLICED molecules:
  Each had supporting reads landing exclusively on the joined exon sequence:
  [Exon1]--[Exon2]--[Exon3]--[Exon4]--[Exon5]--[Exon6]
  e.g., a read spanning the Exon1-Exon2 junction:
  
  Genomic:  ...7,159,300]  GAP (30,977 bp intron1)  [7,190,418...
  mRNA:     ...7,159,300]──────────────────────────── [7,190,418...
                         ↑ intron is ABSENT in mRNA ↑
  
  A read that aligns with a CIGAR like 60M30977N31M crossed the splice
  junction → genome gap matches annotated intron length → exon-exon 
  junction read → SPLICED.

The 8 UNSPLICED molecules:
  Supporting reads land within the 30,977 bp intron 1 
  (chr1:7,159,441–7,190,417), OR span the Exon1→Intron1 boundary.
  → UNSPLICED.
```

**Biological interpretation:** This gene is in **mid-production**. Some copies of Pcmtd1 mRNA have been fully spliced (contributing to translation), while new copies are still being transcribed (unspliced pre-mRNA). The ratio 13:8 ≈ 1.6 is close to steady state, meaning transcription rate ≈ splicing/degradation rate for this gene in this cell.

---

### Example D: *Sntg1* — Actively Upregulating (high unspliced:spliced ratio)

**Observed:** 1 spliced UMI, 47 unspliced UMIs (ratio = 47:1)

```
Sntg1 has 9 transcript isoforms (Sntg1-201 through Sntg1-209).
The gene body spans chr1:8,431,699–9,370,103 (minus strand) — nearly 1 Mb.

47 molecules had reads in intronic regions across this gene body.
Only 1 molecule produced a fully exonic read set.

In RNA velocity terms (La Manno et al. Fig. 1c-d):

  unspliced (u)                    
       ↑                           
       │    * this cell (u=47, s=1)← far above equilibrium line
       │   /
       │  / slope γ = steady-state u/s ratio
       │ /
       │/_________________________ spliced (s)

When u >> γ·s, the cell is ABOVE the equilibrium line.
→ RNA velocity predicts spliced mRNA for Sntg1 is INCREASING.
→ This cell is actively upregulating Sntg1 expression.
```

---

## 7. From Classification to the MTX Files

The `.mtx` files use **Matrix Market format** (sparse matrix). After running kb count, the entries for our demo cell look like:

```
# spliced.mtx (row=barcode_index, col=gene_index, value=UMI_count)
# Cell CTGTGGGAGGTCACCC is barcode index 1 (first row in spliced.barcodes.txt)

row   col    value   → gene
...
1     38     1       → Becn2      (ENSMUSG00000104158.3) — SPLICED only
1     18     13      → Pcmtd1     (ENSMUSG00000051285.18) — spliced molecules
1     81     1       → Sntg1      (ENSMUSG00000025909.17) — spliced molecules
...

# unspliced.mtx
...
1     66     18      → 9330185C12Rik (ENSMUSG00000097648.2) — UNSPLICED only
1     18     8       → Pcmtd1        — unspliced molecules
1     81     47      → Sntg1         — unspliced molecules
...
```

---

## 8. The Full Pipeline in One Diagram

```
SRA Download
    │
    ▼
SRR_ID_2.fastq (R1, 28bp)   SRR_ID_3.fastq (R2, 91bp)
[CTGTGGGAGGTCACCC][UMI]     [GCAGCAGTGACTGGAG...91bp cDNA...]
    │                              │
    │ ←─── kb count ──────────────┘
    │       ↓
    │  1. Extract barcode + UMI from R1
    │  2. Pseudo-align R2 k-mers against combined index
    │  3. Assign to equivalence class:
    │       cdna only     → spliced vote
    │       intron only   → unspliced vote
    │       both / mixed  → ambiguous → discard
    │  4. Collapse by barcode+UMI → one count per molecule
    │  5. Aggregate by gene via t2g.txt
    ▼
counts_unfiltered/
├── spliced.mtx    ← s matrix (mature mRNA counts per cell per gene)
└── unspliced.mtx  ← u matrix (nascent pre-mRNA counts per cell per gene)
    │
    ▼
RNA Velocity (scVelo / velocyto)
  ds/dt = β·u − γ·s
  
  For each gene in each cell, estimate whether s is going UP or DOWN
  based on the current u/s ratio vs the equilibrium γ.
  
  Sntg1 in CTGTGGGAGGTCACCC:  u=47, s=1  →  velocity > 0  (increasing)
  Pcmtd1 in CTGTGGGAGGTCACCC: u=8,  s=13 →  near steady state
```

---

## 9. Why the Distinction Matters

| If you only had spliced counts | If you add unspliced counts |
|---|---|
| Static snapshot of current mRNA levels | Reveals the **direction** mRNA is heading |
| Cannot distinguish actively-transcribing from steady-state cells | Cells with high u relative to s are **upregulating** the gene |
| No sense of developmental time | RNA velocity arrows point toward future cell state |

The key insight of La Manno et al. 2018 is that the unspliced reads — which were previously treated as noise or artifact — are actually a **real-time signal of transcriptional activity**. By reading both the current mRNA level (spliced) and the rate of new production (unspliced), you can estimate the derivative of gene expression and predict where a cell is going, not just where it is.
