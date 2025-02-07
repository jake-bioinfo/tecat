# TECAT

TECAT (Telomere End Chromosome Assaying Tool)

## Features

- Reference sequence analysis for de novo telomere motif discovery
- Highly customizable motif assignment for repettitive sequence matching
- Telomeric read identification with high fidelity
- Accurate telomere length analysis
- All step are parallelized 
- Mapping telomeres back to reference for chromosome end-specific lengths

## Installation

```bash
git clone https://github.com/jake-bioinfo/tecat.git
R CMD build tecat
R CMD INSTALL tecat
```

## Quickstart
### Workflow
Split fastq files 
`auto_split(fastq_file)`
Find organism/sample specific telomere motifs
`telomere_motif(reference_file)`
Find telomeres in samples
`telo_search(split_fastq_file_list)`
Cut potential telomeres into rolling windows
`sliding_window_parallel(telomere_file)`
Calculate frequencies of windows
`frequencies(windows)`
Determine optimal thresholds
`determine_threshold(telomere_frequencies)`
`optimal_thresholds(threshold_dataframe)`
Truncate and map remainder read sequences
`truncate_file(telomere_file)`
`map(truncated_reads, reference_file)`
Plot results!


## License

GPL-3.0 License - see LICENSE file for details