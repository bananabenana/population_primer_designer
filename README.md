# population_primer_designer
This bioinformatics tool designs primers against single genomes or large populations. It identifies shared genomic regions in desired genomes ensuring these regions are absent from unwanted genomes. Then primers are designed on these regions and in silico PCR is run to confirm sizes.
This is done in parallel, so it means you can run 10,000 genomes through this and it will efficiently generate primers.

## Quick start

### Installation
```bash
# Clone github repo
git clone https://github.com/bananabenana/population_primer_designer

# Install dependencies via conda or mamba
cd population_primer_designer; conda env create -f population_primer_designer_env.yaml
```

## Run on example genomes
```bash
# Activate environment
conda activate population_primer_designer_env

# Run autopilot mode
python population_primer_designer.py autopilot \
  --genome_list example_data/genome_list.txt \
  --output_dir example_autopilot_out \
  --amplicon_length 500 \
  --max_amplicon_diff 100 \
  --num_return 20 \
  --threads 4
```

## Dependencies
See `population_primer_designer_env.yaml`
