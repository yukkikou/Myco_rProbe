bash myco_rprobe_single.sh -s Aspergillus_fumigatus -g example/Aspergillus_fumigatus.ASM265v1.dna.toplevel.fa.gz -t 20 -o single_af -m /media/share/data4/container/Py_env/pan_fungi/pan_fungi3.sif

bash myco_rprobe_multi.sh -g example/fungi_genome.list -t 10 -o pathgen_fungi \
    -m /media/share/data4/container/Py_env/pan_fungi/pan_fungi3.sif \
    -r data/rRNA_list2.tsv > log

bash myco_rprobe_evaluate.sh -g example/batch1_lib_config.tsv -t 10 -o batch1_mapping \
    -m /media/share/data4/container/Py_env/pan_fungi/pan_fungi5.sif \
    -r pathgen_fungi > log

bash myco_rprobe_circ.sh -g example/batch2_lib_config.tsv -t 10 -o batch_mapping \
    -m /media/share/data4/container/CIRIquant/ciriquant-v1.1.3.sif

bash myco_rprobe_circseq.sh -g example/circ_seq.tsv -t 5 -o batch1_mapping/Cryptococcus/4_target -m /media/share/data4/container/Py_env/pan_fungi/pan_fungi5.sif

# plus design
sige /media/share/data4/container/seqkit/seqkit-2.2.0.sif seqkit grep -r -p "sliva_28S_JALPTJ010000291" ../../../pathgen_fungi/sep_28s.fasta > Mucor_28S_sliva_28S.fa
bash ../../../src/plus_design/s01_run_mpile.sh Mucor

sige Py_env/pan_fungi/pan_fungi5.sif python ../../../src/plus_design/s02_mismatch_pileup.py Mucor_28s.pileup Mucor.28S.depth 59

##############################
pan_sig=/media/share/data4/container/Py_env/pan_fungi/pan_fungi3.sif
sige $pan_sig python src/group_probes.py Af/28S_60_0.8_60/cdhit_final.alig.out \
    Af/sep_28s.fasta \
    data/fungi_taxonomy_info.tsv \
    tmp.tsv filter.tsv

sige $pan_sig python src/group_probes_plus.py pathgen_fungi/5S_60_0.85_60/cdhit_final.alig.out \
    pathgen_fungi/sep_5s.fasta \
    data/fungi_taxonomy_info.tsv \
    data/fungi_taxonomy_nodes.tsv \
    group_5s.tsv filter_5s.tsv --num_workers 10


sige $pan_sig python src/grouped_coverage.py filter.tsv \
    0.7 phylum \
    Af/sep_28s.fasta \
    Af/28S_60_0.8_60/cdhit_final.alig.out \
    grouped.tsv \
    Af/28S_60_0.8_60/cdhit_final.fasta \
    grouped_probes

sige $pan_sig python src/decision_tree.py \
    --nodes data/fungi_taxonomy_nodes.tsv \
    --info data/fungi_taxonomy_info.tsv \
    --grouped-dir pathgen_fungi/5S_60_0.85_60/ \
    --probes-fasta pathgen_fungi/5S_60_0.85_60/rrna_probes.fasta \
    --species-map-fasta pathgen_fungi/sep_28s.fasta.merge \
    --fungi-root-name "Fungi" \
    --output-file fungal_tree_report.txt

sige $pan_sig python src/decision_tree.py \
    --nodes data/fungi_taxonomy_nodes.tsv \
    --info data/fungi_taxonomy_info.tsv \
    --grouped-dir pathgen_fungi/28S_59_0.8_59/ \
    --probes-fasta pathgen_fungi/28S_59_0.8_59/rrna_probes.fasta \
    --species-map-fasta pathgen_fungi/sep_28s.fasta.merge \
    --fungi-root-name "Fungi" \
    --output-file 28s_fungal_tree_report.txt


sige $pan_sig python src/draw_tree.py \
    --nodes data/fungi_taxonomy_nodes.tsv \
    --info data/fungi_taxonomy_info.tsv \
    --grouped-dir pathgen_fungi/28S_60_0.8_60/ \
    --probes-fasta pathgen_fungi/28S_60_0.8_60/rrna_probes.fasta \
    --species-map-fasta pathgen_fungi/sep_28s.fasta.merge \
    --fungi-root-name "Fungi" \
    --output-file 28s_fungal_tree_report.txt \
    --output-graph fungal_tree.png

sige $pan_sig Rscript src/stat_coverage.R \
    /media/share/data5/1610305236/Panfungi/script/Myco_rProbe/batch1_mapping \
    Cryptococcus_neoformans \
    28S
