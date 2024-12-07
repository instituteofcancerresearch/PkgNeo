# PkgStopGain
A package for the calculation of up- and down-stream stop gain calculations introduced by a mutation.

**Authors**  
- Rachel Alcraft, Institute of Cancer Research, Research Software Engineering-Scientific Computing, London UK  
- Lorena Magraner, Institute of Cancer Research, Gene Function-Breast Cancer Research, London UK  

## Installation
This can be installed as a pakage from github. 
If it is private the additional private authentication steps need to be followed below.

You will need the package `remotes` installed.  
    
```R
install.packages("remotes")
remotes::install_github("instituteofcancerresearch/PkgStopGain", force=TRUE)
```

## Documentation  

This package has a few gene sequences that can be retrieved, and the functions to calculate the upstream and downstream stop codons for a mutation along with information about the mutation and up-and-down shifts.        

### Find the genes with pre-saved sequences
```R
PkgStopGain::get_pre_saved_genes_list()
>>> [1] "BRCA1" "BRCA2" "RAD51C" "RAD51D" "PALB2" 
```

### Get the gene sequence for a gene (if pre-saved)
```R
gene_seq <- PkgStopGain::get_gene_seq("BRCA1")
>>> [1] "atggatttatctgctcttcgcgttgaagaagtacaaaatgtcattaatgctatgcagaaaatcttagagtgtc..."
```

### Investigate the mutation
```R
gene_seq <- PkgStopGain::get_gene_seq("BRCA1")
mut <- "c.5267_5273delinsAA"
mut_obj <- PkgStopGain::apply_mutation(mut, gene_seq)
>>> 
$type
[1] "delins"
$start
[1] 5267
$end
[1] 5273
$offset
[1] -5
$modulo
[1] 1
$ins
[1] "AA"
$del
[1] "aggacag"
$dna_orig
[1] "atggattt..."
$ok
[1] TRUE
```

### Get the stop gain information
**The downstream stop gain** is the first stop codon after the mutation (downstream of the mutation, to the right).  
The sequence will shift from the wildtype downstream until there is a stop.  
`WWWWWWWWWXXXXXXXXXX*YYYYYYYYYYYYY`  
Where
- W is the wiltype sequence  
- X is shifted sequence after the mutation  
- Y is the shifted sequence after the stop codon  
We want to return the `XXXXXXXXX*` and the equivalent dna sequence.  

**The upstream stop gain** is the first stop codon before the mutation (upstream of the mutation, to the left).  
This assumes that the protein sequence is read from right to left.  
The sequence will shift from the wildtype upstream until there is a stop.  
`YYYYYYYYY*XXXXXXXXXXWWWWWWWWWW`  
Where
- W is the wiltype sequence,
- X is shifted sequence after the mutation (right to left),
- Y is the shifted sequence after the stop codon (right to left)
We want to returnn the `*XXXXXXX` and the equivalent dna sequence.  


```R
gene_seq <- PkgStopGain::get_gene_seq("BRCA1")
mut <- "c.5266dupC"
stop_gains <- PkgStopGain::stop_gain_factors(mut,gene_seq)
>>>
$ok
[1] TRUE
$dna_orig
[1] "atggatttatctgctcttcgcgttgaa..."
$aa_orig
[1] "MDLSALRVEEVQN..."
$type
[1] "dup"
$new_dna
[1] "atggatttatctgc..."
$down_all
[1] "MDLSALRVEEVQNVINAMQKILEC..."
$down_stop_gain
[1] "PGQKDLQGARNLLLWALHQHAHRSTGMDGTAVWCFCGEGAFIIHPWHRCPPNCGCAARCLDRGQWLPCNWADV*"
$down_dna_stop_gain
[1] "cCaggacagaaagatcttcagggggctagaaatctgttgctatgggcccttcaccaacatgcccacagatcaactggaatggatggtacagctgtgtggtgcttctgtggtgaaggagctttcatcattcacccttggcacaggtgtccacccaattgtggttgtgcagccagatgcctggacagaggacaatggcttccatgcaattgggcagatgtgtga"
$up_all
[1] "WIYLLFALKKYKMSLMLCRKS*SV..."
$up_stop_gain
[1] "*MSMILKSEEMWSMEETTKVQSEQENP"
$up_dna_stop_gain
[1] "tgaatgagcatgattttgaagtcagaggagatgtggtcaatggaagaaaccaccaaggtccaaagcgagcaagagaatccc"
      
```
---  

## Additional developer information

#### Private authentication
For a private repo you will need to authenticate with your github credentials. You will need to have been given access to this repo.   
You will need the following R packages installed: `remotes`, `devtools`

1. set config  
`usethis::use_git_config(user.name = "YourName", user.email = "your@mail.com")`

2. Go to github page to generate token  
`usethis::create_github_token() `

3. paste your PAT into pop-up that follows...  
`credentials::set_github_pat()`

4. now remotes::install_github() will work  
`remotes::install_github("instituteofcancerresearch/PkgStopGain")`

---  

### Working locally with the package
To work and develop with the package it can also be cloned locally and you can work with it and install it as a local package.

To do this you will need to clone the repo and then install the package locally.  

```Bash
git clone https://github.com/instituteofcancerresearch/PkgStopGain.git
or 
git clone git@github.com:instituteofcancerresearch/PkgStopGain.git
```

Then you can install the package locally using the following command.  

```R
devtools::install(pkg="path/to/PkgStopGain")
```

You can now make changes in the package locally and then re-install the package to see the changes.  

### Developing and updating the package
When developing the package you can make changes to the functions and then update the package, and test the package from a diufferent R script side by side.  
The `document` function ensures the correct functions are exported, the `build` function builds the package and the `install` function installs the package.  

```R
devtools::document(pkg="path/to/PkgStopGain")
devtools::build(pkg="path/to/PkgStopGain")
devtools::install(pkg="path/to/PkgStopGain",force=TRUE)'
```




