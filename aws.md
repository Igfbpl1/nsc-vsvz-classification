# AWS EC2 Setup — RNA Velocity Pipeline (GSM8253799)

## Step 1 — Check GSM8253799 file sizes first

Before launching, go to this URL in your browser to see the exact sizes:
```
https://www.ncbi.nlm.nih.gov/Traces/study/?acc=GSM8253799
```
Note the total size of all SRR runs — this determines your EBS volume size.

---

## Step 2 — Launch EC2 Instance (browser)

1. Go to [console.aws.amazon.com](https://console.aws.amazon.com) → **EC2** → **Launch Instance**
2. Set these options:
   - **Name**: `vsvz-velocity`
   - **AMI**: Ubuntu 26.04 LTS (`ami-0b6d9d3d33ba97d99`)
   - **Instance type**: `r6i.2xlarge` (64GB RAM, 8 vCPUs) — On-Demand, not Spot
   - **Key pair**: `vsvz-key` (already created, saved to `~/.ssh/vsvz-key.pem`)
   - **Storage**: 800 GiB gp3 root volume
3. Click **Launch Instance**

**Instance details (already launched):**
- Public DNS: `ec2-3-80-204-115.compute-1.amazonaws.com`
- Key: `~/.ssh/vsvz-key.pem`

---

## Step 3 — Connect via Terminal (no extra software needed)

Open Terminal on your Mac:

```bash
# fix key permissions (already done)
chmod 400 ~/.ssh/vsvz-key.pem

# SSH in
ssh -i ~/.ssh/vsvz-key.pem ubuntu@ec2-3-80-204-115.compute-1.amazonaws.com
```

---

## Step 4 — Install tools on EC2

```bash
# system packages
sudo apt update && sudo apt install -y git curl tmux

# sra-tools (for downloading FASTQs)
sudo apt install -y sra-toolkit

# kallisto
sudo apt install -y kallisto

# uv (Python package manager)
curl -Ls https://astral.sh/uv/install.sh | sh
source ~/.bashrc

# create virtual environment (python3-venv required on Ubuntu 26.04)
sudo apt install python3.14-venv -y
python3 -m venv ~/venv
source ~/venv/bin/activate
echo "source ~/venv/bin/activate" >> ~/.bashrc

# install bustools + kb-python into venv
uv pip install bustools kb-python
```

---

## Step 5 — Rebuild the reference on EC2

This avoids uploading the large `index.idx` from your Mac — download the genome directly from Ensembl instead:

```bash
mkdir -p ~/project/sra_runs && cd ~/project/sra_runs

# download genome + GTF (same files used locally)
wget https://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/release-110/gtf/mus_musculus/Mus_musculus.GRCm39.110.gtf.gz

# build kb reference (takes ~30 min)
kb ref \
  -i index.idx \
  -g t2g.txt \
  -f1 cdna.fa \
  -f2 intron.fa \
  -c1 cdna_t2c.txt \
  -c2 intron_t2c.txt \
  --workflow lamanno \
  Mus_musculus.GRCm39.dna.primary_assembly.fa.gz \
  Mus_musculus.GRCm39.110.gtf.gz
```

Run with `nohup` to keep it running if you close Terminal:
```bash
nohup kb ref \
  -i index.idx \
  -g t2g.txt \
  -f1 cdna.fa \
  -f2 intron.fa \
  -c1 cdna_t2c.txt \
  -c2 intron_t2c.txt \
  --workflow lamanno \
  Mus_musculus.GRCm39.dna.primary_assembly.fa.gz \
  Mus_musculus.GRCm39.110.gtf.gz > ref.log 2>&1 &

# check progress
tail -f ref.log

# check if still running
jobs
```

---

## Step 6 — Upload processed.h5ad from your Mac

This is the only file you need to upload (your Mac → EC2). Open a **new Terminal tab** on your Mac:

```bash
scp -i ~/.ssh/vsvz-key.pem \
  /Users/chandra/development/nsc-vsvz-classification/outputs/processed.h5ad \
  ubuntu@ec2-3-80-204-115.compute-1.amazonaws.com:~/project/outputs/
```

---

## Step 7 — Download GSM8253799 FASTQs on EC2

SRR runs for GSM8253799 (total ~77GB):
- SRR28912869 (32GB)
- SRR28912870 (32GB)
- SRR28912871 (13GB)

```bash
cd ~/project/sra_runs

nohup bash -c "prefetch SRR28912869 SRR28912870 SRR28912871 && \
  fasterq-dump --split-files --include-technical \
  SRR28912869 SRR28912870 SRR28912871 --outdir ." \
  > fastq.log 2>&1 &

# check progress
tail -f fastq.log
```

**Note:** `--include-technical` is required to get barcode reads (`_3` files). Without it, fasterq-dump only writes cDNA reads (`_4` files) and kb count will fail.

---

## Step 8 — Run kb count

Add GSM8253799 to `run_kb_count.sh` (fill in SRR IDs from Step 7), then:

```bash
bash run_kb_count.sh GSM8253799

# clean up immediately after to free space
rm sra_runs/*.fastq
rm sra_runs/kb_output_GSM8253799/*.bus
rm sra_runs/kb_output_GSM8253799/matrix.ec
rm sra_runs/kb_output_GSM8253799/transcripts.txt
```

---

## Step 9 — Run the velocity pipeline

Upload the Python code, then run:

```bash
scp -i ~/.ssh/vsvz-key.pem \
  /Users/chandra/development/nsc-vsvz-classification/rna_velocity_pipeline.py \
  ubuntu@ec2-3-80-204-115.compute-1.amazonaws.com:~/project/

cd ~/project
tmux new -s velocity
# auto-shutdown when pipeline completes so billing stops
uv run python rna_velocity_pipeline.py && sudo shutdown -h now
```

**Billing note:**
- Instance auto-stops when pipeline finishes
- Stop = compute billing pauses but EBS disk still charges (~$1.50/day for 800GB)
- Terminate = everything deleted, all billing stops
- Download results first, then terminate

---

## Step 10 — Download results + terminate instance

```bash
# from your Mac
scp -r -i ~/.ssh/vsvz-key.pem \
  ubuntu@ec2-3-80-204-115.compute-1.amazonaws.com:~/project/outputs/velocity/ \
  /Users/chandra/development/nsc-vsvz-classification/outputs/velocity_GSM8253799/
```

Then go to EC2 console → **Terminate Instance** so billing stops.

---

## Notes

- Always run long steps inside `tmux` so they survive SSH disconnects
- Check SRA Run Selector for GSM8253799 SRR IDs and total size before Step 2
- EBS volume = 2x total FASTQ size to account for kb count intermediates
- After kb count completes, delete FASTQs + BUS files immediately to free space
- Terminate the instance as soon as results are downloaded — billing stops immediately
