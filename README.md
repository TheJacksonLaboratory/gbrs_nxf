Run GBRS (https://github.com/churchill-lab/gbrs) on paired-end RNAseq data using
Nextflow pipelines.

This is adapted to Sumner cluster (slurm, singularity) and uses the
jaxreg.jax.org containers repository, follow the User Guide there to add it to
your library.

To test it, install nextflow and run:
```bash
nextflow run TheJacksonLaboratory/gbrs_nxf -profile singularity,slurm \
  -resume --params-file <parameter_yaml_file>
```
Make sure the input files are in the format: `ID_*R{1,2}*.fastq.gz`

The parameter file is a `yaml` formatted file that contains at least the
parameter entries:
```yaml
fastqR1: [path]
fastqR2: [path]
outputDir: [path]
generation: '[G0-G40]'
sex: [M|F]
```
Plus any other optional parameters, which can be shown by passing the `--help`
argument i.e. `nextflow run TheJacksonLaboratory/gbrs_nxf --help`

