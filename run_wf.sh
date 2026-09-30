nextflow run . -profile test,docker \
  --input /workspaces/test_data/oncotrseq_test/samplesheet.csv \
  --outdir /workspaces/test_data/oncotrseq_test/results \
  -w /workspaces/test_data/oncotrseq_test/work