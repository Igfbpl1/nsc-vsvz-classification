# AWS EC2 Reference — RNA Velocity Pipeline

Quick reference for running the V-SVZ velocity pipeline on AWS EC2.

---

## Quick Reference — Connecting & Monitoring

### Current instance
- **Public DNS**: `ec2-3-80-204-115.compute-1.amazonaws.com`
- **Key**: `~/.ssh/vsvz-key.pem`
- **Type**: `r6i.2xlarge` (64GB RAM, 8 vCPUs)
- **Disk**: 800 GiB gp3
- **AMI**: Ubuntu 26.04 LTS

### Connect

```bash
# fix key permissions (one-time)
chmod 400 ~/.ssh/vsvz-key.pem

# SSH in
ssh -i ~/.ssh/vsvz-key.pem ubuntu@ec2-3-80-204-115.compute-1.amazonaws.com
```

### Copy files Mac → EC2

```bash
# single file
scp -i ~/.ssh/vsvz-key.pem \
  /Users/chandra/development/nsc-vsvz-classification/<file> \
  ubuntu@ec2-3-80-204-115.compute-1.amazonaws.com:~/project/

# directory recursive
scp -r -i ~/.ssh/vsvz-key.pem \
  /path/to/dir \
  ubuntu@ec2-3-80-204-115.compute-1.amazonaws.com:~/project/

# multiple files
scp -i ~/.ssh/vsvz-key.pem \
  process_sample.sh process_all_remaining.sh rna_velocity_pipeline.py \
  ubuntu@ec2-3-80-204-115.compute-1.amazonaws.com:~/project/
```

### Copy files EC2 → Mac

```bash
# from your Mac
scp -i ~/.ssh/vsvz-key.pem \
  ubuntu@ec2-3-80-204-115.compute-1.amazonaws.com:~/project/sra_runs/GSM*_counts.tgz \
  /Users/chandra/development/nsc-vsvz-classification/sra_runs/

# directory recursive
scp -r -i ~/.ssh/vsvz-key.pem \
  ubuntu@ec2-3-80-204-115.compute-1.amazonaws.com:~/project/outputs/velocity/ \
  /Users/chandra/development/nsc-vsvz-classification/outputs/velocity_remote/
```

### Monitoring commands (run on EC2)

```bash
# overall disk usage
df -h

# folder-level usage (find what's eating space)
du -sh ~/project/sra_runs/*
du -sh ~/project/sra_runs/kb_output_*/

# running processes
ps aux | grep -E "kb|kallisto|bustools|fasterq"
jobs                    # background jobs in current shell
top                     # interactive process viewer (q to quit)
htop                    # nicer if installed (sudo apt install htop)

# memory usage
free -h

# CPU info
nproc
cat /proc/cpuinfo | grep "model name" | head -1

# follow log files
tail -f ~/project/sra_runs/ref.log
tail -f ~/project/sra_runs/kb_count.log
tail -f ~/project/process_all.log

# check if specific file exists & size
ls -lh ~/project/sra_runs/index.idx
ls -lh ~/project/sra_runs/GSM*_counts.tgz
```

### Background jobs that survive SSH disconnect

```bash
# launch in background, redirect logs
nohup bash some_script.sh > some_script.log 2>&1 &

# check if it's still running
jobs
ps aux | grep some_script

# follow log (Ctrl+C to stop following, process keeps running)
tail -f some_script.log

# kill a background job by PID
kill <PID>
```

### Shutdown / cost control

```bash
# auto-shutdown on completion (saves money)
some_long_command && sudo shutdown -h now

# manual stop from EC2 (instance can be restarted, EBS still charges)
sudo shutdown -h now

# permanent stop — must terminate from AWS Console:
# Console → EC2 → select instance → "Terminate Instance"
```

**Billing reminder:**
- **Running** = ~$0.40/hr (r6i.2xlarge on-demand) + ~$0.05/hr EBS
- **Stopped** = compute charge stops; **EBS still bills** at ~$1.50/day for 800GB
- **Terminated** = all billing stops; volume deleted

---

## One-Time Setup (when launching a new instance)

### Launch instance via AWS Console

1. Go to [console.aws.amazon.com](https://console.aws.amazon.com) → **EC2** → **Launch Instance**
2. Set:
   - **Name**: `vsvz-velocity`
   - **AMI**: Ubuntu 26.04 LTS
   - **Instance type**: `r6i.2xlarge` (64GB RAM, 8 vCPUs)
   - **Key pair**: `vsvz-key` (or create new and save `.pem` to `~/.ssh/`)
   - **Storage**: 800 GiB gp3 root volume
   - **Purchasing option**: On-Demand (Spot is cheaper but can be interrupted)
3. Click **Launch Instance**

### Install dependencies (one-time per instance)

```bash
# system packages
sudo apt update && sudo apt install -y git curl tmux htop sra-toolkit kallisto entrez-direct

# uv (Python package manager)
curl -Ls https://astral.sh/uv/install.sh | sh
source ~/.bashrc

# create venv (python3-venv required on Ubuntu 26.04)
sudo apt install python3.14-venv -y
python3 -m venv ~/venv
source ~/venv/bin/activate
echo "source ~/venv/bin/activate" >> ~/.bashrc

# install Python packages
uv pip install bustools kb-python

# create project dir
mkdir -p ~/project/sra_runs ~/project/outputs
```

### Build kb reference (one-time per instance, ~30 min)

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

tail -f ref.log    # Ctrl+C to stop tailing
```

---

## Processing One Sample (small/test runs)

### Look up SRR sizes first

Before downloading, check sizes at:
```
https://www.ncbi.nlm.nih.gov/Traces/study/?acc=<GSM_ID>
```

### Download FASTQs — in parallel for speed

For a sample with multiple SRR runs, download all in parallel to maximize throughput:

```bash
cd ~/project/sra_runs

# parallel prefetch + fasterq-dump per SRR
nohup fasterq-dump --split-files --include-technical SRR_ID1 --outdir . > fastq_1.log 2>&1 &
nohup fasterq-dump --split-files --include-technical SRR_ID2 --outdir . > fastq_2.log 2>&1 &
nohup fasterq-dump --split-files --include-technical SRR_ID3 --outdir . > fastq_3.log 2>&1 &

# follow all three
tail -f fastq_1.log fastq_2.log fastq_3.log
```

**Important:** `--include-technical` is required to get the barcode reads (`_3` files). Without it, fasterq-dump only writes cDNA reads (`_4` files) and kb count fails.

### Run kb count directly

```bash
nohup kb count \
  -i index.idx \
  -g t2g.txt \
  -x 10xv3 \
  -o kb_output_GSM8647353 \
  -c1 cdna_t2c.txt \
  -c2 intron_t2c.txt \
  --workflow lamanno \
  SRR31443695_3.fastq SRR31443695_4.fastq \
  SRR31443696_3.fastq SRR31443696_4.fastq \
  SRR31443697_3.fastq SRR31443697_4.fastq \
  SRR31443700_3.fastq SRR31443700_4.fastq \
  > kb_count.log 2>&1 &

tail -f kb_count.log
```

### Clean up to free disk after kb count finishes

```bash
cd ~/project/sra_runs

# delete FASTQs
rm SRR*.fastq

# delete BUS files and intermediates (keep only counts_unfiltered/)
cd kb_output_<GSM>
rm -f *.bus matrix.ec transcripts.txt inspect*.json 10x_version3_whitelist.txt
rm -rf tmp
cd ..

# package counts_unfiltered for download
tar czf <GSM>_counts.tgz -C kb_output_<GSM> counts_unfiltered
```

---

## Processing All Remaining Samples (batch mode)

Use the helper scripts `process_sample.sh` and `process_all_remaining.sh` in the project root.

### Upload scripts to EC2

```bash
# from Mac
scp -i ~/.ssh/vsvz-key.pem \
  process_sample.sh process_all_remaining.sh \
  ubuntu@ec2-3-80-204-115.compute-1.amazonaws.com:~/project/
```

### Test with one sample first

```bash
# on EC2
cd ~/project
bash process_sample.sh GSM8253792
```

This auto-looks-up SRR IDs via NCBI Entrez, downloads + runs kb count + cleans up. Output: `~/project/sra_runs/GSM8253792_counts.tgz`.

### Run all 7 remaining samples sequentially

```bash
# on EC2
cd ~/project
nohup bash process_all_remaining.sh > process_all.log 2>&1 &

# monitor
tail -f process_all.log

# check disk between samples
df -h
```

Estimated runtime: ~3-4 hours per sample × 7 ≈ 24-28 hours unattended.

### Download tarballs back to Mac

```bash
# from Mac
scp -i ~/.ssh/vsvz-key.pem \
  ubuntu@ec2-3-80-204-115.compute-1.amazonaws.com:~/project/sra_runs/GSM*_counts.tgz \
  /Users/chandra/development/nsc-vsvz-classification/sra_runs/

# extract each
cd /Users/chandra/development/nsc-vsvz-classification
for tgz in sra_runs/GSM*_counts.tgz; do
    GSM=$(basename "$tgz" _counts.tgz)
    mkdir -p sra_runs/kb_output_${GSM}
    tar xzf "$tgz" -C sra_runs/kb_output_${GSM}
done
```

---

## Running the Velocity Pipeline on EC2

If you want to run the full velocity pipeline on EC2 (avoids Mac RAM/CPU bottleneck for 10-sample analysis):

```bash
# upload pipeline + processed.h5ad
scp -i ~/.ssh/vsvz-key.pem \
  /Users/chandra/development/nsc-vsvz-classification/rna_velocity_pipeline.py \
  /Users/chandra/development/nsc-vsvz-classification/markers.py \
  ubuntu@ec2-3-80-204-115.compute-1.amazonaws.com:~/project/

scp -i ~/.ssh/vsvz-key.pem \
  /Users/chandra/development/nsc-vsvz-classification/outputs/processed.h5ad \
  ubuntu@ec2-3-80-204-115.compute-1.amazonaws.com:~/project/outputs/

# on EC2: run with auto-shutdown
cd ~/project
nohup bash -c "python rna_velocity_pipeline.py && sudo shutdown -h now" > velocity.log 2>&1 &
tail -f velocity.log
```

Instance auto-stops when done — billing pauses (EBS still bills).

---

## Cost Estimate for Full 10-Sample Run

- **kb count per sample**: ~3-4 hrs × 7 samples = ~25 hrs
- **Velocity pipeline**: ~5 hrs
- **Total compute**: ~30 hrs × $0.40/hr = **$12**
- **EBS storage** during run (800GB × ~2 days): **$5**
- **Data transfer out** (~10GB results to download): **$0.90**
- **Total: ~$18**

Spot pricing on r6i.2xlarge can reduce this 60-70%, but interruption risk during long runs makes On-Demand safer.

---

## Tips & Gotchas

- **Always run long steps with `nohup` or `tmux`** — otherwise they die when SSH disconnects.
- **`sudo apt install entrez-direct`** provides `esearch`/`elink`/`efetch` — needed for the auto-SRR-lookup in `process_sample.sh`.
- **The old `_4` files from a non-`--include-technical` run** may be incompatible with the new `_3` files. If kb count complains about mismatched read counts, delete and regenerate.
- **EBS volume size cannot be decreased**, only increased. Start conservative; if you run out, expand via Console.
- **Spot interruption mid-kb-count** loses ~3 hours of work for that sample. The cleanup logic in `process_sample.sh` handles partial state by re-checking for `counts_unfiltered/spliced.mtx`.
- **`df -h` regularly** during multi-sample runs to make sure you're not approaching disk full.
- **Terminate the instance** after downloading results — a stopped instance still bills for EBS.
