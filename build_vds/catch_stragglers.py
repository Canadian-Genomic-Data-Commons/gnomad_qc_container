import os
import hail as hl

# Paths
SCRATCH = "/projects/HostSeq/hostseq_jc/scratch"
VDS_OUT = os.path.join(SCRATCH, "stragglers.vds")

# Input lists
with open("sample_list.txt") as fh:
    sample_list = set(line.strip() for line in fh if line.strip())

with open("VDS_content.txt") as fh:
    vds_content = set(line.strip() for line in fh if line.strip())

with open("consent_withdrawn.txt") as fh:
    consent_withdrawn = set(line.strip() for line in fh if line.strip())

# Determine stragglers
stragglers = sorted(sample_list - vds_content - consent_withdrawn)

print(f"Found {len(stragglers)} stragglers eligible for import.")

if not stragglers:
    print("No stragglers to process, exiting.")
    exit(0)

# Map samples -> gVCF path
def get_gvcf_path(sample_id):
    return f"/projects/gvcfs/{sample_id}.g.vcf.gz"

gvcf_paths = {s: get_gvcf_path(s) for s in stragglers}

# Import and combine
print("Starting Hail import...")
hl.init(tmp_dir=os.path.join(SCRATCH, "hail_tmp"))

combiner = hl.vds.new_combiner(
    output_path=f"file://{VDS_OUT}",
    temp_path=os.path.join(SCRATCH, "straggler_tmp"),
    gvcf_paths=gvcf_paths,
    reference_genome="GRCh38",
    import_interval_size=50_000_000,
    gvcf_batch_size=len(gvcf_paths),
    target_records=1_000_000,
    force=True,
)

combiner.run()
print(f"Straggler VDS written to {VDS_OUT}")

# Append to VDS_content.txt so these won't reappear
with open("VDS_content.txt", "a") as fh:
    for s in stragglers:
        fh.write(s + "\n")

print("Updated VDS_content.txt")
