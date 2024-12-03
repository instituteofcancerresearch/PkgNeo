# PkgNeo
A package for the neoantigens calculation at the Institute of Cancer Research.  

**Author** Rachel Alcraft, Institute of Cancer Research, London, UK  

### Installation
This can be installed as a pakage from github. As this is a private repo you will need to authenticate with your github credentials. You will need to have been given access to this repo.   

YOu will need the following R packages installed: `remotes`  

1. remotes::install_github() will work  
`remotes::install_github("instituteofcancerresearch/PkgNeo")`  

### Usage
This package has a few gene sequences that can be retrieved, and the functions to calculate the upstrream and downstram stop codins for a mutation along with information about the mutation and up-and-down shifts.        

```R
library(PkgNeo)
# Get the gene sequence for a gene
gene_seq <- get_gene_seq("brca1")
# Get the up-down factors
up_down_factors <- up_down_factors("c.5266dupC",gene_seq)

# The results include:

# Original seqeunce information
up_down_factors$dna_orig
up_down_factors$aa_orig      

# Mutation information
up_down_factors$ins_piece
up_down_factors$del_piece
up_down_factors$new_dna

# Upstream shift information
up_down_factors$up_aa
up_down_factors$up_stub
up_down_factors$up_stop

# Downstream shift information
up_down_factors$down_aa
up_down_factors$down_stub
up_down_factors$down_stop

# metric information
up_down_factors$distance      
      
```

# Working locally with the package
To work and develop with the package it can also be cloned locally and you can work with it and install it as a local package.

To do this you will need to clone the repo and then install the package locally.  

```Bash
git clone https://github.com/instituteofcancerresearch/PkgNeo.git
or 
git clone git@github.com:instituteofcancerresearch/PkgNeo.git
```

Then you can install the package locally using the following command.  

```R
devtools::install(pkg="path/to/PkgNeo")
```

You can now make changes in the package locally and then re-install the package to see the changes.  

