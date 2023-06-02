# duplication_detector
During the long-term evolution, ancestors of current species had underdone few to several waves of whole genome duplication (WGD) events, which contributes to the genetic and morphological innovation. Nowadays researchers adopt many accepted approaches, e.g. chromosome-level synteny and Ks distribution of protein genes, to identify and understand the signals of WGD with currently available genomes. Approximate one decade ago, Jiao et al. (2011) developed a phylogenomic  approach and found two ancient WGD event predating the origin of angiosperms and seed plants. This approach takes its unique advantage for studies with a long time span (>300 Myr), because faint chromosome-level synteny and saturated Ks cannot provide precise if any information under such evolutionary scale.

Based on the method adopted in Jiao et al. (2011), I provide a user-friendly script (infer_duplication_node.pl) to fulfil similar phylogenomic analysis.

## Dependencies
- Perl5

- Bioperl

## Usage
Users are supposed to provide two files to the script. The bootstrap_threshold <-b> is used to determine the reliability of nodes in your tree.

```console
Usage: infer_duplication_node.pl -b <bootstrap_threshold> -i <gene_tree_with_bs> -o <output_dir> -l <species_list_with_phylogeny_order>
```

- A gene tree with bootstrap value, e.g. Orthofinder output.
- A species list with phylogeny order as prior information.
For example, 

|         |                                                                  |
| ------- | ---------------------------------------------------------------- |
| Ulva    | 1 |
| Volvox  | 1 |
| Chlamydomonas | 1   |
| Chlorokybus | 2 |
| Mesostigma | 2 |
| Chara | 3 |
| Mesotaenium | 4 |
| Spirogloea | 4 |
| Marchantia | 5 |
| Sphagnum | 5 |
| Physcomitrium | 5 |
| Hypnum | 5 |
| Anthoceros | 5 |
| Adiantum | 6 |
| Cycas | 6 |
| Arabidopsis | 6 |
![image](https://github.com/yurunxian/duplication_detector/assets/48025559/a904648c-07c1-4235-9ce2-f5a58b908cd9)

## Workflow
![image](https://github.com/yurunxian/duplication_detector/assets/48025559/bda6706b-f623-42b2-a47a-89b80382d51c)

## Example
![image](https://github.com/yurunxian/duplication_detector/assets/48025559/580bfda8-d867-4b2d-9b0a-8cffb47fc1b2)
![image](https://github.com/yurunxian/duplication_detector/assets/48025559/bc5a61fa-4812-4994-a3d6-fca33a4ff417)
![image](https://github.com/yurunxian/duplication_detector/assets/48025559/aa21db7b-9a80-4e73-ac1d-fde6795dba68)
