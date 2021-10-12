conda install -c anaconda biopython
conda install -c conda-forge scikit-bio
conda install -c bioconda htseq pysam pyfaidx
git clone https://github.com/brentp/align.git
#Add a package: git+https://github.com/brentp/align.git

conda env export -f environment.yml
