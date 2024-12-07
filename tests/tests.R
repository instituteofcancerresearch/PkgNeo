


seq <- PkgStopGain::get_gene_seq("brca1")
mutations <- c("c.5266dupC","c.1016delA","c.1121_1121delC")
results <- PkgStopGain::stop_gain_factors("c.1016delA",seq)
