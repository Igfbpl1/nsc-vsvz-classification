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

## 2. The Reference: Exon & Intron Coordinates

Two files are needed before any counting can happen:

```
Mus_musculus.GRCm39.dna.primary_assembly.fa   # 2.7B bp genome sequence
Mus_musculus.GRCm39.110.gtf                   # exon coordinates per transcript
```

The GTF lists only **exon** positions. Introns are never written explicitly — they are the gaps between consecutive exons of the same transcript. For example, *Becn2* on chr1:

```
GTF entries (+ strand):

Exon 1:  chr1  175,747,895 ──────────── 175,749,254   (1,360 bp)
                                                  ↕
Intron 1: chr1  175,749,255 ─── 175,749,362   (108 bp, inferred gap)
                                                  ↕
Exon 2:  chr1  175,749,363 ─── 175,749,791   (429 bp)
```

Everything between two exon entries for the same transcript is intronic by elimination. That is the complete basis for classifying reads as exonic or intronic.

---

## 3. Building the Index & Running the Count

### 3a. kb ref — build two separate indices from the GTF

kb-python performs **pseudo-alignment** — it matches read k-mers against a pre-built index rather than doing full genome alignment. For spliced/unspliced quantification it builds two indices from the same GTF:

```bash
kb ref \
  -i index.idx -g t2g.txt -f1 cdna.fa -f2 intron.fa \
  --workflow lamanno \
  Mus_musculus.GRCm39.dna.primary_assembly.fa \
  Mus_musculus.GRCm39.110.gtf
```

| File | What it actually is | Role |
|------|---------------------|------|
| `cdna.fa` | A FASTA of **sequences**: each record is a header (`>ENSMUST… gene_id:… gene_name:… chr/start/end`) followed by the transcript's exons joined together (introns removed) | The **spliced** sequence database — reads matching here are counted spliced |
| `intron.fa` | A FASTA of **sequences**: each record is a header (`-I.N` suffix) followed by an intron sequence of that transcript | The **unspliced** sequence database — reads matching here are counted unspliced |
| `cdna_t2c.txt` | A plain **list of IDs only** — one spliced transcript ID per line (`ENSMUST…`), no sequence | "Capture list" marking which index entries are **spliced** |
| `intron_t2c.txt` | A plain **list of IDs only** — one intron-capture ID per line (`ENSMUST…-I.N`), no sequence | "Capture list" marking which index entries are **unspliced** |
| `t2g.txt` | One line per index entry → its gene (both `ENSMUST…` and `ENSMUST…-I.N` point to the same gene) | Maps a matched ID back to its **gene** |
| `index.idx` | `cdna.fa` + `intron.fa` compiled into one searchable structure | The single pseudo-alignment target reads are matched against |

Two distinctions worth keeping straight: (1) the **`.fa` files hold the sequences** (with an ID + coordinates in each header), whereas the **`t2c` files hold only the ID lists** — they are the table of contents of each `.fa`, not sequence files; and (2) the **`t2c` files only classify a hit as spliced vs unspliced**, while **`t2g.txt` is what maps a hit to a gene** — these are separate jobs.

The `-I.N` suffix is the signal that ties them together: `ENSMUST00000180837.2` is the spliced entry for *9330185C12Rik* (listed in `cdna_t2c.txt`), while `ENSMUST00000180837.2-I.1` is its intronic capture entry (listed in `intron_t2c.txt`). Both map to the same gene via `t2g.txt`, but a read matching the `-I.1` entry is counted unspliced.

### 3b. kb count — align, deduplicate, and produce count matrices

```bash
kb count \
  -i sra_runs/index.idx -g sra_runs/t2g.txt \
  -x 10xv3 --lamanno \
  -c1 sra_runs/cdna_t2c.txt \
  -c2 sra_runs/intron_t2c.txt \
  -o sra_runs/kb_output_GSM8253792/ \
  SRR_ID_2.fastq SRR_ID_3.fastq
```

Internally:

1. **Extract:** First 16 bp of R1 = cell barcode; next 12 bp = UMI. Barcode corrected against 10x v3 whitelist.
2. **Pseudo-align:** R2 broken into 31-mers, matched against `index.idx`. Each match is a vote for `cdna` or `intron`.
3. **BUS file:** Every (barcode, UMI, equivalence class) tuple written to disk.
4. **Collapse by barcode+UMI:** All reads sharing a barcode+UMI = one molecule. Votes aggregated → spliced / unspliced / ambiguous.
5. **Count:** Molecules tallied per gene via `t2g.txt`, written to `.mtx`.

```
kb_output_GSM8253792/counts_unfiltered/
├── spliced.mtx    — 13,411,965 nonzero entries across 816,643 barcodes × 56,941 genes
└── unspliced.mtx  — 11,422,859 nonzero entries across 712,804 barcodes × 56,941 genes
```

> The 816,643 barcodes are **unfiltered** — they include ~800k empty droplets that captured
> ambient RNA. Real cells (~4,600) sit above the UMI count knee and are selected in the
> next analysis step. The ratio 11.4M / 13.4M ≈ 0.85 unspliced:spliced is consistent with
> La Manno et al.'s observation that 15–25% of 10x reads carry intronic sequence.

---

## 4. Classification & Examples

### The decision tree (applied per read, then aggregated per molecule)

Once a read's coordinates are returned by the aligner, the tool checks them against the GTF exon intervals:

```
Given: read maps to [read_start, read_end] on chrX

Is read entirely within one exon          YES ──→ SPLICED vote
[exon_start, exon_end]?                           (exonic)
                    │ NO
                    ↓
Does read overlap the gap between         YES ──→ UNSPLICED vote
two consecutive exons?                            (spans exon-intron boundary)
                    │ NO
                    ↓
Is read entirely between two exons        YES ──→ UNSPLICED vote
with no exon overlap?                             (purely intronic)
                    │ NO
                    ↓
Does CIGAR N-skip match annotated         YES ──→ SPLICED vote
intron length in GTF?                             (exon-exon junction read)
                    │ NO
                    ↓
                AMBIGUOUS → discard

After all reads in a molecule's supporting set are voted:

  All votes SPLICED                              → spliced.mtx
  All votes UNSPLICED (across all isoforms)      → unspliced.mtx
  Mixed                                          → discarded
```

All examples below are from cell **`CTGTGGGAGGTCACCC`** — the highest-UMI cell in this sample (45,354 spliced + 41,369 unspliced molecules across 7,579 genes).

---

### Example A: *Becn2* — SPLICED only

**Counts:** 1 spliced UMI, 0 unspliced UMIs → `spliced.mtx` row 1, col 38, value 1

```
chr1 (+ strand):
[══════════ Exon 1 ══════════]          [══ Exon 2 ══]
175,747,895              175,749,254  175,749,363  175,749,791
                                    ↑
                           108 bp intron (175,749,255–175,749,362)

Read R2:  [═══════ 91 bp ═══════]
          175,747,950         175,748,040

Decision: entirely within Exon 1 → SPLICED vote
          cdna index hit only → molecule = SPLICED
```

*Mature, fully processed mRNA. The intron has been removed; protein is being translated.*

---

### Example B: *9330185C12Rik* — UNSPLICED only

**Counts:** 0 spliced UMIs, 18 unspliced UMIs → `unspliced.mtx` row 1, col 66, value 18

```
chr1 (− strand, transcription runs right → left):

[Exon 3]              [Exon 2]         [══════ Exon 1 ══════]
113,819,115-113,819,778  113,837,065-113,837,093  113,891,549-113,891,864
           └──── Intron 2 ────┘└───────── Intron 1 (54,455 bp) ──────────┘

Scenario i — boundary read:
                                      [═══ 91 bp ═══]
                                     113,891,500  113,891,590
  → 50 bp in Exon 1, 41 bp in Intron 1 → spans boundary → UNSPLICED vote

Scenario ii — purely intronic read:
              [═══ 91 bp ═══]
             113,850,000  113,850,090
  → Entirely within Intron 1 (54,455 bp gap) → UNSPLICED vote
  → Arises from secondary priming at intronic polyT sequences

All 18 molecules fall into scenario i or ii → molecule = UNSPLICED
```

*Nascent pre-mRNA. Introns present. Gene is being actively transcribed.*

---

### Example C: *Pcmtd1* — BOTH (different molecules, same gene)

**Counts:** 13 spliced UMIs, 8 unspliced UMIs → both matrices, row 1, col 18

These are 21 distinct molecules (distinct UMIs). Each is classified independently.

```
Pcmtd1-201, chr1 (+ strand), 6 exons:

[Exon 1]      ←── Intron 1 (30,977 bp) ───→      [Exon 2] ... [Exon 6]
7,159,144-7,159,440                              7,190,418-7,190,839

The 13 SPLICED molecules — reads hit joined exon sequence:
  CIGAR 60M30977N31M: 60 bp in Exon 1, skip 30,977 bp, 31 bp in Exon 2
  Genome gap matches annotated intron → exon-exon junction → SPLICED

The 8 UNSPLICED molecules — reads hit Intron 1:
  Reads land in chr1:7,159,441–7,190,417 or span Exon1→Intron1 boundary
  → UNSPLICED
```

*Mid-production: some copies are mature and translating (spliced), new copies are still being transcribed (unspliced). Ratio 13:8 ≈ 1.6 — near steady state.*

---

### Example D: *Sntg1* — actively upregulating

**Counts:** 1 spliced UMI, 47 unspliced UMIs → `spliced.mtx` value 1, `unspliced.mtx` value 47

```
Sntg1, chr1 (− strand), 9 isoforms, gene body ~938 kb

47 molecules: reads in intronic regions → UNSPLICED
 1 molecule:  reads fully exonic        → SPLICED

Phase portrait for this cell:

  u │    * (u=47, s=1) ← far above steady-state line
    │   /
    │  / slope = γ/β
    │ /
    │/──────────────── s

u >> γ·s  →  ds/dt > 0  →  spliced mRNA is INCREASING
```

*This gene is being rapidly upregulated in this cell right now.*

---

### MTX output summary for cell `CTGTGGGAGGTCACCC`

The `.mtx` files store only nonzero entries as `barcode_index gene_index count`. This cell is index 1 in `spliced.barcodes.txt`:

```
spliced.mtx                          unspliced.mtx
─────────────────────────────────    ─────────────────────────────────────
row  col  val  gene                  row  col  val  gene
1    38   1    Becn2     (spliced    1    66   18   9330185C12Rik (unspliced
                          only)                          only)
1    18   13   Pcmtd1               1    18   8    Pcmtd1
1    81   1    Sntg1                1    81   47   Sntg1
```

---

## 8. The Full Pipeline in One Diagram

```
╔══════════════════════════════════════════════════════════════════════╗
║  STEP 1 — SEQUENCING                                                 ║
║                                                                      ║
║  SRA Download: fasterq-dump SRR_ID --split-files                     ║
║                                                                      ║
║  SRR_ID_2.fastq (R1, 28bp)      SRR_ID_3.fastq (R2, 91bp)            ║
║  ┌──────────────────┬──────────┐  ┌───────────────────────────────┐  ║
║  │ cell barcode     │  UMI     │  │ cDNA sequence (91 bp)         │  ║
║  │ CTGTGGGAGGTCACCC │GATTACAGT │  │ GCAGCAGTGACTGGAGAAGCCAT...    │  ║
║  │ (16 bp)          │ (12 bp)  │  │                               │  ║
║  └──────────────────┴──────────┘  └───────────────────────────────┘  ║
║       ↑ WHO + WHICH MOLECULE            ↑ WHAT WAS TRANSCRIBED       ║
╚══════════════════════════════════════════════════════════════════════╝
                              │
                              ▼
╔══════════════════════════════════════════════════════════════════════╗
║  STEP 2 — REFERENCE BUILDING (done once, reused for all samples)     ║
║                                                                      ║
║  Inputs:                                                             ║
║    Mus_musculus.GRCm39.dna.primary_assembly.fa  (2.7B bp genome)     ║
║    Mus_musculus.GRCm39.110.gtf                  (exon coordinates)   ║
║                                                                      ║
║  GTF defines exon boundaries — introns are inferred gaps:            ║
║                                                                      ║
║    Becn2 chr1:  [═══ Exon1 ═══]──intron──[═Exon2═]                   ║
║                  175,747,895                       175,749,791       ║
║                                                                      ║
║  kb ref builds two indices:                                          ║
║    cdna.fa    → exon sequences joined (spliced transcripts)          ║
║    intron.fa  → full gene body incl. introns (pre-mRNA)              ║
║                                                                      ║
║  Output: index.idx  t2g.txt  cdna_t2c.txt  intron_t2c.txt            ║
╚══════════════════════════════════════════════════════════════════════╝
                              │
                              ▼
╔══════════════════════════════════════════════════════════════════════╗
║  STEP 3 — PSEUDO-ALIGNMENT & COUNTING  (kb count)                    ║
║                                                                      ║
║  For every read pair:                                                ║
║    1. Extract barcode + UMI from R1                                  ║
║    2. Correct barcode against 10x v3 whitelist                       ║
║    3. Break R2 into 31-mers, match against index.idx                 ║
║                                                                      ║
║  Classify each read by coordinate vs GTF interval:                   ║
║                                                                      ║
║   read coords vs GTF exon table                                      ║
║       │                                                              ║
║       ├─ entirely within one exon?          → EXONIC (spliced vote)  ║
║       ├─ overlaps exon-intron boundary?     → INTRONIC (unspl. vote) ║
║       ├─ entirely within intron gap?        → INTRONIC (unspl. vote) ║
║       ├─ CIGAR N matches annotated intron?  → EXONIC (spliced vote)  ║
║       └─ mixed across isoforms?             → AMBIGUOUS → discard    ║
║                                                                      ║
║  Collapse by barcode+UMI (one count per unique molecule):            ║
║    all reads exonic       → 1 spliced molecule                       ║
║    all reads intronic     → 1 unspliced molecule                     ║
║    mixed                  → discarded                                ║
║                                                                      ║
║  Aggregate by gene via t2g.txt                                       ║
║                                                                      ║
║  Output:  counts_unfiltered/                                         ║
║    spliced.mtx    13,411,965 nonzero entries  ← s  (mature mRNA)     ║
║    unspliced.mtx  11,422,859 nonzero entries  ← u  (nascent pre-mRNA)║
║                                                                      ║
║  NOTE: 816,643 barcodes total — mostly empty droplets.               ║
║  Real cells (~4,600) are identified by UMI count knee plot.          ║
╚══════════════════════════════════════════════════════════════════════╝
                              │
                              ▼
╔══════════════════════════════════════════════════════════════════════╗
║  STEP 4 — CELL FILTERING & NORMALISATION  (scanpy / scVelo)          ║
║                                                                      ║
║  Filter to real cells above UMI knee (discard empty droplets)        ║
║  Normalise counts by cell size                                       ║
║  Select highly variable genes                                        ║
║                                                                      ║
║  Result: two dense matrices for ~4,600 cells × ~3,000 HVGs           ║
║    s[cells × genes]   u[cells × genes]                               ║
╚══════════════════════════════════════════════════════════════════════╝
                              │
                              ▼
╔══════════════════════════════════════════════════════════════════════╗
║  STEP 5 — DEGRADATION RATE ESTIMATION  (scVelo)                      ║
║                                                                      ║
║  The model:   du/dt = α − β·u                                        ║
║               ds/dt = β·u − γ·s                                      ║
║                                                                      ║
║    α = transcription rate  ─┐                                        ║
║    β = splicing rate        ├─ never directly measured               ║
║    γ = degradation rate    ─┘                                        ║
║                                                                      ║
║  Degradation cannot be observed from reads — a degraded mRNA         ║
║  leaves no trace in the sequencing data.                             ║
║                                                                      ║
║  Instead, γ is INFERRED from the steady-state constraint:            ║
║    At steady state:  du/dt = 0  and  ds/dt = 0                       ║
║    → β·u = γ·s   →   u = (γ/β)·s                                     ║
║    → u vs s across cells forms a LINE with slope γ/β                 ║
║                                                                      ║
║    u │          · ·                                                  ║
║      │       · · /← steady-state cells                               ║
║      │     · · /   cluster on this line                              ║
║      │   · · /     slope = γ/β  ← estimated by regression            ║
║      │  · · /      on upper quantile of scatter                      ║
║      │────────────────── s                                           ║
║                                                                      ║
║  scVelo fits this slope gene-by-gene across all cells.               ║
║  Output: one γ value per gene                                        ║
╚══════════════════════════════════════════════════════════════════════╝
                              │
                              ▼
╔══════════════════════════════════════════════════════════════════════╗
║  STEP 6 — VELOCITY COMPUTATION  (scVelo)                             ║
║                                                                      ║
║  For each gene in each cell:                                         ║
║    velocity = ds/dt = β·u − γ·s                                      ║
║                                                                      ║
║  Interpretation — where is each cell relative to steady-state line?  ║
║                                                                      ║
║    u │    ·  ← ABOVE line: u too high for current s                  ║
║      │   /     β·u > γ·s  →  ds/dt > 0  →  s is INCREASING           ║
║      │  /── steady-state line                                        ║
║      │ /   ·  ← BELOW line: u too low for current s                  ║
║      │/      β·u < γ·s  →  ds/dt < 0  →  s is DECREASING             ║
║      │──────────────── s                                             ║
║                                                                      ║
║  Real examples from cell CTGTGGGAGGTCACCC:                           ║
║                                                                      ║
║    Gene      s     u    position        velocity                     ║
║    ──────────────────────────────────────────────────────────        ║
║    Sntg1     1    47    far above line  ds/dt >> 0  (rising fast)    ║
║    Pcmtd1   13     8    near line       ds/dt ≈ 0   (steady state)   ║
║    Becn2     1     0    below line      ds/dt < 0   (degrading)      ║
║                                                                      ║
║  Output: velocity matrix [cells × genes]                             ║
╚══════════════════════════════════════════════════════════════════════╝
                              │
                              ▼
╔══════════════════════════════════════════════════════════════════════╗
║  STEP 7 — PROJECTION ONTO EMBEDDING                                  ║
║                                                                      ║
║  Velocity vectors (one per cell, across thousands of genes) are      ║
║  projected onto a low-dimensional embedding (UMAP / t-SNE / PCA).    ║
║                                                                      ║
║  Each arrow on the UMAP shows:                                       ║
║    direction → which cell state this cell is moving toward           ║
║    length    → how fast the transition is occurring                  ║
║                                                                      ║
║  OPC ──────────────────────────────────────────────────────────      ║
║   ·  → → →                                                           ║
║   · · → → → →   ← velocity arrows pointing toward COP fate           ║
║   · · · → → → →                                                      ║
║  TAP          COP ──→ OL                                             ║
║   ·  ·  ·                                                            ║
║   · · · (no arrows = cells not actively transitioning)               ║
╚══════════════════════════════════════════════════════════════════════╝
```