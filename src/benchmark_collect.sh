#!/bin/bash
mkdir -p ./benchmarks/consumption
cat `find ./benchmarks/star/ -name "*star_align.benchmark"` | sed '1!{/h:m:s/d;}' > benchmarks/consumption/star_align.benchmark
cat `find ./benchmarks/centrifuge/ -name "*centrifuge.benchmark"` | sed '1!{/h:m:s/d;}' > benchmarks/consumption/microbiota.benchmark
cat `find ./benchmarks/optitype/ -name "*optitype.benchmark"` | sed '1!{/h:m:s/d;}' > benchmarks/consumption/optitype.benchmark
cat `find ./benchmarks/rseqc/gene_body_cvg/ -name "*gene_body_cvg_qc.benchmark"` | sed '1!{/h:m:s/d;}' > benchmarks/consumption/rseqc_genebody_cvg.benchmark
cat `find ./benchmarks/rseqc/read_quality/ -name "*read_quality.benchmark"` | sed '1!{/h:m:s/d;}' > benchmarks/consumption/rseqc_read_qual.benchmark
cat `find ./benchmarks/rseqc/read_distrib/ -name "*read_distrib_qc_matrix.benchmark"` | sed '1!{/h:m:s/d;}' > benchmarks/consumption/rseqc_read_distrib.benchmark
cat `find ./benchmarks/rseqc/tin_score/ -name "*tin_score.benchmark"` | sed '1!{/h:m:s/d;}' > benchmarks/consumption/rseqc_tin_score.benchmark
cat `find ./benchmarks/rseqc/insert_size/ -name "*insert_size.benchmark"` | sed '1!{/h:m:s/d;}' > benchmarks/consumption/rseqc_insert_size.benchmark
cat `find ./benchmarks/rseqc/junction_saturation/ -name "*junction_saturation.benchmark"` | sed '1!{/h:m:s/d;}' > benchmarks/consumption/rseqc_junction_saturation.benchmark
cat `find ./benchmarks/salmon/ -name "*salmon.benchmark"` | sed '1!{/h:m:s/d;}' > benchmarks/consumption/quantification.benchmark
cat `find ./benchmarks/varscan/ -name "*somatic.benchmark"` | sed '1!{/h:m:s/d;}' > benchmarks/consumption/varscan.benchmark
cat `find ./benchmarks/vep/ -name "*vep_snp_annot.benchmark"` | sed '1!{/h:m:s/d;}' > benchmarks/consumption/vep.benchmark
cat `find ./benchmarks/pvacseq/ -name "*neoantigen_pvacseq.benchmark"` | sed '1!{/h:m:s/d;}' > benchmarks/consumption/pvacseq.benchmark
cat `find ./benchmarks/fusion/ -name "*star_fusion.benchmark"` | sed '1!{/h:m:s/d;}' > benchmarks/consumption/fusion.benchmark
cat `find ./benchmarks/msisensor/ -name "*msisensor.benchmark"` | sed '1!{/h:m:s/d;}' > benchmarks/consumption/msisensor.benchmark
cat `find ./benchmarks/trust4/ -name "*trust4.benchmark"` | sed '1!{/h:m:s/d;}' > benchmarks/consumption/trust4.benchmark
cat `find ./benchmarks/deseq2/ -name "*deseq2.benchmark"` | sed '1!{/h:m:s/d;}' > benchmarks/consumption/deseq2.benchmark
cat `find ./benchmarks/lisa/ -name "*lisa.benchmark"` | sed '1!{/h:m:s/d;}' > benchmarks/consumption/lisa.benchmark
cp benchmarks/ssgsea/ssgsea.benchmark benchmarks/consumption/
cp benchmarks/wgcna/WGCNA.benchmark benchmarks/consumption/
cp benchmarks/immunedeconv/deconv.benchmark benchmarks/consumption/
cp benchmarks/tide/tide_score.benchmark benchmarks/consumption/
cp benchmarks/msi_est/msi_est_score.benchmark benchmarks/consumption/