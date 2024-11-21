
# Creating a dev package in the same folder
```bash
mamba create --name env-r-neo -c bioconda r-base=4.3
mamba activate env-r-neo
mamba install r-devtools
Rscript -e 'devtools::create("../PkgNeo", open=TRUE)'
```
