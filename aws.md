# V-SVZ Pipeline — Run Reference

---

## Running Locally (Mac)

### Prerequisites
- `uv` for Python
- kb-python installed via uv (already in `pyproject.toml` — `uv sync` installs it)
- kallisto installed via Homebrew; bustools compiled via `kb compile` (no Homebrew formula)
- Run all kb commands as `uv run kb ...`

```bash
uv sync                  # installs kb-python and all deps
brew install cmake       # required for kb compile
brew install kallisto    # pre-built ARM binary for Apple Silicon
uv run kb compile bustools --cmake-arguments="-DCMAKE_POLICY_VERSION_MINIMUM=3.5"
                         # compiles bustools from source into the venv
uv run kb info           # verify kb is available
```

> kb-python bundles binaries for Linux but not Mac. `bustools` has no Homebrew formula and
> `kb compile all` fails on cmake 4.x due to a bifrost sub-dependency. Solution: use Homebrew
> kallisto + compile bustools separately with the cmake policy flag.
>
> Compiled bustools lands at:
> `.venv/lib/python3.12/site-packages/kb_python/bins/compiled/bustools/bustools`
> Pass this full path with `--bustools` on every `kb ref` and `kb count` call.
> kallisto is at `/opt/homebrew/bin/kallisto` on Apple Silicon (`/usr/local/bin/` on Intel).

### Steps in order

**0. Build the reference (one-time, ~30 min, ~16GB RAM)**

If `sra_runs/index.idx` doesn't exist, you have two options:

**Option A — copy from EC2 (faster if EC2 has it):**
```bash
scp -i ~/.ssh/vsvz-key.pem \
  ubuntu@ec2-3-80-204-115.compute-1.amazonaws.com:~/project/sra_runs/index.idx \
  ubuntu@ec2-3-80-204-115.compute-1.amazonaws.com:~/project/sra_runs/t2g.txt \
  ubuntu@ec2-3-80-204-115.compute-1.amazonaws.com:~/project/sra_runs/cdna_t2c.txt \
  ubuntu@ec2-3-80-204-115.compute-1.amazonaws.com:~/project/sra_runs/intron_t2c.txt \
  sra_runs/
```

**Option B — build from scratch:**
```bash
cd /Users/viv/development/nsc-vsvz-classification/sra_runs

# download genome + annotation (~800MB + ~40MB)
curl -O https://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
curl -O https://ftp.ensembl.org/pub/release-110/gtf/mus_musculus/Mus_musculus.GRCm39.110.gtf.gz

cd /Users/viv/development/nsc-vsvz-classification

uv run kb ref \
  -i sra_runs/index.idx -g sra_runs/t2g.txt \
  -f1 sra_runs/cdna.fa -f2 sra_runs/intron.fa \
  -c1 sra_runs/cdna_t2c.txt -c2 sra_runs/intron_t2c.txt \
  --workflow lamanno \
  --kallisto /opt/homebrew/bin/kallisto \
  --bustools /Users/viv/development/nsc-vsvz-classification/.venv/lib/python3.12/site-packages/kb_python/bins/compiled/bustools/bustools \
  sra_runs/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz \
  sra_runs/Mus_musculus.GRCm39.110.gtf.gz
```

**1. Download FASTQs**

```bash
cd /Users/viv/development/nsc-vsvz-classification/sra_runs

nohup fasterq-dump --split-files --include-technical SRR_ID1 --outdir . > SRR_ID1.log 2>&1 &
nohup fasterq-dump --split-files --include-technical SRR_ID2 --outdir . > SRR_ID2.log 2>&1 &
```

**Identify which files are R1 and R2** — use `tail -2` to see the last read in each file and check its length:

```bash
tail -2 SRR_ID_1.fastq   # shows sequence + quality of last read
tail -2 SRR_ID_2.fastq
tail -2 SRR_ID_3.fastq
tail -2 SRR_ID_4.fastq   # if dual-indexed (4 files)
```

| Length | What it is | Use? |
|--------|-----------|------|
| 10 bp | I1 or I2 — sequencer index reads for lane demultiplexing | No |
| 28 bp | R1 — cell barcode (16bp) + UMI (12bp) | **Yes** |
| 90–91 bp | R2 — cDNA read, maps to genome | **Yes** |

Single-indexed (3 files): `_1`=I1, `_2`=R1, `_3`=R2 → use `_2` + `_3`
Dual-indexed (4 files): `_1`=I1, `_2`=I2, `_3`=R1, `_4`=R2 → use `_3` + `_4`

**2. kb count for each sample**

```bash
cd /Users/viv/development/nsc-vsvz-classification

uv run kb count \
  -i sra_runs/index.idx \
  -g sra_runs/t2g.txt \
  -x 10xv3 --workflow lamanno \
  --kallisto /opt/homebrew/bin/kallisto \
  --bustools /Users/viv/development/nsc-vsvz-classification/.venv/lib/python3.12/site-packages/kb_python/bins/compiled/bustools/bustools \
  -c1 sra_runs/cdna_t2c.txt \
  -c2 sra_runs/intron_t2c.txt \
  -o sra_runs/kb_output_<GSM> \
  sra_runs/<SRR>_3.fastq sra_runs/<SRR>_4.fastq   # dual-indexed layout
  # OR sra_runs/<SRR>_2.fastq sra_runs/<SRR>_3.fastq  # single-indexed layout
```

**3. Run the full pipeline**

```bash
cd /Users/viv/development/nsc-vsvz-classification

# first run (builds everything)
uv run python run_pipeline.py

# force rebuild if inputs changed
uv run python run_pipeline.py --rebuild
```

This runs in order:
- `build_processed()` → `outputs/processed.h5ad`
- `train_ol_classifier` → `outputs/trigger_genes.csv`, `outputs/ol_commitment.csv`
- `tap_analysis` → `outputs/out_of_sample_tap_comparison_test2.csv`
- `velocity_build.build_velocity()` → `outputs/velocity/velocity_combined.h5ad` (~30 min)
- `rna_velocity_pipeline` → stream plots + driver CSVs in `outputs/velocity/`
- `compare_tap_fate_methods` → `outputs/tap_fate_*.csv`

**4. Individual steps** (if you only need to re-run part)

```bash
uv run python -c "import velocity_build; velocity_build.build_velocity()"
uv run python rna_velocity_pipeline.py
uv run python compare_tap_fate_methods.py
```

---

## EC2 — Connecting & Monitoring

### Current instance
- **Public DNS**: `ec2-3-80-204-115.compute-1.amazonaws.com`
- **Key**: `~/.ssh/vsvz-key.pem`
- **Type**: `r6i.2xlarge` (64GB RAM, 8 vCPUs)
- **Disk**: 800 GiB gp3
- **AMI**: Ubuntu 26.04 LTS

### Connect
```bash
chmod 400 ~/.ssh/vsvz-key.pem   # one-time
ssh -i ~/.ssh/vsvz-key.pem ubuntu@ec2-3-80-204-115.compute-1.amazonaws.com
```

### Copy files Mac → EC2
```bash
scp -i ~/.ssh/vsvz-key.pem \
  /Users/viv/development/nsc-vsvz-classification/<file> \
  ubuntu@ec2-3-80-204-115.compute-1.amazonaws.com:~/project/

# multiple files
scp -i ~/.ssh/vsvz-key.pem \
  velocity_build.py rna_velocity_pipeline.py preprocess.py markers.py \
  ubuntu@ec2-3-80-204-115.compute-1.amazonaws.com:~/project/
```

### Copy files EC2 → Mac
```bash
# kb_output tarballs
scp -i ~/.ssh/vsvz-key.pem \
  "ubuntu@ec2-3-80-204-115.compute-1.amazonaws.com:~/project/sra_runs/GSM*_counts.tgz" \
  /Users/viv/development/nsc-vsvz-classification/sra_runs/

# velocity outputs
scp -r -i ~/.ssh/vsvz-key.pem \
  ubuntu@ec2-3-80-204-115.compute-1.amazonaws.com:~/project/outputs/velocity/ \
  /Users/viv/development/nsc-vsvz-classification/outputs/
```

### Extract tarballs (on Mac after download)
```bash
cd /Users/viv/development/nsc-vsvz-classification
for tgz in sra_runs/GSM*_counts.tgz; do
    GSM=$(basename "$tgz" _counts.tgz)
    mkdir -p sra_runs/kb_output_${GSM}
    tar xzf "$tgz" -C sra_runs/kb_output_${GSM}
done
```

### Monitoring (on EC2)
```bash
df -h                                        # disk usage
du -sh ~/project/sra_runs/kb_output_*/      # per-sample size
ps aux | grep -E "kb|fasterq"               # running processes
tail -f ~/project/sra_runs/kb_count.log     # kb count progress
free -h                                      # memory
```

### Background jobs (survive SSH disconnect)
```bash
nohup <command> > output.log 2>&1 &
tail -f output.log     # follow (Ctrl+C to stop following)
ps aux | grep <name>   # check still running
kill <PID>             # stop if needed
```

### Shutdown / cost control
```bash
sudo shutdown -h now   # stop instance (EBS still bills)
# To terminate fully: AWS Console → EC2 → Terminate Instance
```

**Billing:** Running ~$0.40/hr + ~$0.05/hr EBS. Stopped: EBS still ~$1.50/day. Terminated: all billing stops.

---

## EC2 — One-Time Setup

```bash
# system packages
sudo apt update && sudo apt install -y git curl tmux htop sra-toolkit kallisto entrez-direct

# uv
curl -Ls https://astral.sh/uv/install.sh | sh && source ~/.bashrc

# python venv
sudo apt install python3.14-venv -y
python3 -m venv ~/venv
source ~/venv/bin/activate
echo "source ~/venv/bin/activate" >> ~/.bashrc

# python packages
uv pip install kb-python scvelo cellrank xgboost shap scanpy

# project dirs
mkdir -p ~/project/sra_runs ~/project/outputs/velocity
```

### Build kb reference (one-time, ~30 min)
```bash
cd ~/project/sra_runs

wget https://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/release-110/gtf/mus_musculus/Mus_musculus.GRCm39.110.gtf.gz

nohup kb ref \
  -i index.idx -g t2g.txt -f1 cdna.fa -f2 intron.fa \
  -c1 cdna_t2c.txt -c2 intron_t2c.txt \
  --workflow lamanno \
  Mus_musculus.GRCm39.dna.primary_assembly.fa.gz \
  Mus_musculus.GRCm39.110.gtf.gz > ref.log 2>&1 &

tail -f ref.log
```

---

## EC2 — Processing a Sample

### 1. Download FASTQs
```bash
cd ~/project/sra_runs

# run in parallel for speed
nohup fasterq-dump --split-files --include-technical SRR_ID1 --outdir . > SRR_ID1.log 2>&1 &
nohup fasterq-dump --split-files --include-technical SRR_ID2 --outdir . > SRR_ID2.log 2>&1 &

tail -f SRR_ID1.log SRR_ID2.log
ps aux | grep fasterq   # confirm running
```

**Check layout** after download (single-indexed → use `_2`+`_3`; dual-indexed → use `_3`+`_4`):
```bash
awk 'NR==2{print length($0); exit}' SRR_ID1_1.fastq
awk 'NR==2{print length($0); exit}' SRR_ID1_2.fastq
awk 'NR==2{print length($0); exit}' SRR_ID1_3.fastq
```

### 2. Run kb count
```bash
# single-indexed (_2=R1 barcode, _3=R2 cDNA)
nohup kb count \
  -i index.idx -g t2g.txt -x 10xv3 \
  -o kb_output_<GSM> \
  -c1 cdna_t2c.txt -c2 intron_t2c.txt \
  --workflow lamanno \
  SRR_ID1_2.fastq SRR_ID1_3.fastq \
  SRR_ID2_2.fastq SRR_ID2_3.fastq \
  > kb_count.log 2>&1 &

tail -f kb_count.log
```

### 3. Package and clean up
```bash
tar czf <GSM>_counts.tgz -C kb_output_<GSM> counts_unfiltered

# free disk: delete FASTQs and intermediates
rm SRR*.fastq
cd kb_output_<GSM>
rm -f *.bus matrix.ec transcripts.txt inspect*.json 10x_version3_whitelist.txt
rm -rf tmp
```

---

## Tips

- **Always use `nohup`** — processes die when SSH disconnects without it.
- **FASTQ layout matters** — dual-indexed (4 files) uses `_3`+`_4`; single-indexed (3 files) uses `_2`+`_3`. Verify with `awk` before running kb count.
- **`df -h` regularly** during multi-sample runs.
- **Terminate, don't just stop** — stopped instances still bill for EBS (~$1.50/day for 800GB).
- **Local kb count** requires `--kallisto /opt/homebrew/bin/kallisto --bustools /opt/homebrew/bin/bustools` on Mac (EC2 doesn't need these as they're in PATH).
