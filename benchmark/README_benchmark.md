### Validation of the 2-step classification approach

To assess the performance of the 2-step classification approach, 3 metagenomes were simulated from the same simulated community to approximate 2x150bp PE sequencing, 1x80bp (short read sequencing), and 1x40bp reads (aDNA).


For the DB build and classify commands, we tuned the following parameters: maximum DB size (RAM requirements), kraken confidence threshold (implemented in the classify command), RTL threshold (calculated using conifer).


For the coarse DB classification step, we calculated the percentage of reads correctly assigned to each domain (Archaea, Bactera, Eukaryotes, Viruses, random sequences). 
This is shown in the figures ```*_cm_domain.pdf```. 
The precision, recall, and F1 score for each domain are shown in ```*_stats_domain.pdf```. 
Likewise, we also assessed the performance by larger taxonomic groups, i.e. Archaea, Bactera, Viruses, Fungi, Protists (incl. algae), higher plants, Metazoa.
These results are in ```*_cm_taxongroups.pdf``` and ```*_stats_taxongroups.pdf```
The corresponding legend for these plots is in ```sim_legend_coarse.pdf```.


For the second classification step using the high-resolution database only for prokaryotes and viruses, we calculated the number and percentage of correctly classified reads, prokaryotic or viral reads incorrectly assigned to other prokaryotic or viral taxa, eukaryotic or random reads incorrectly assigned to a prokaryotic or viral taxon, and unassigned reads at each taxonomic rank. 
The absolute read numbers are shown in a stacked plot in ```*_taxlevels_abs.pdf``` and the stacked percentages in ```*_taxlevels_rel.pdf```. 


### Conclusion

Setting a confidence threshold in the kraken2 classify command does not give good results: the reduction in false positives comes at a high cost of fals negatives.
Instead, setting a RTL threshold can substantially improve results
Limiting the DB size will lower classification performance to some degree, but for the maximum DB sizes test here, this was acceptable.


Without the coarse classification, a high proportion of non-prokaryotic and non-viral reads would have been assigned to a taxon using only the high-resolution database.
This noise can be reduced by setting a stringest RTL filter, however not to the degree that was achieved with the coarse DB sorting step. 
Converserly, the recall was slightly lower with the 2-step classification approach. 
In the end, the small cost of slightly more false negatives in the 2-step approach is worth the gain of the strongly reduced proportion of false positives.


So far, only the performance of the high-resolution DB for prokaryotes and viruses was validated.


### Assessment of dereplication approach



