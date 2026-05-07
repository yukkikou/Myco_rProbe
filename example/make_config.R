library(tidyverse)
sam_info = read_tsv("/media/share/data8/volume1/1610305236_PanFungi/batch6/sample_info.tsv", col_names = c("sampleid","sample_name","f1","f2"))
anno_info = read_tsv("fusa_lib_ano.tsv", col_names = c("strain","sampleid","genome","gtf"))

left_join(anno_info, sam_info) %>%
	select(strain,sample_name,f1,f2,genome,gtf) %>%
	write_tsv("fusa_lib_config.tsv", col_names = F)
