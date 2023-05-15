# Init
mkdir microbetag_env
cd microbetag_env
PWD=$(pwd)


# Get Julia

# Check if Julia present 
julia -v
retval=$?
# If no, get it
if [ $retval -ne 0 ]; then
    echo "Install Julia"
    wget https://julialang-s3.julialang.org/bin/linux/x64/1.7/julia-1.7.1-linux-x86_64.tar.gz
    tar -zxvf julia-1.7.1-linux-x86_64.tar.gz
    echo "export PATH='$PATH:$PWD/julia-1.7.1/bin'" >> ~/.bashrc
    source ~/.bashrc
    rm *.tar.gz
else
    echo "Julia already installed"
fi

# Get FlashWeave
julia -e 'using Pkg;Pkg.add("FlashWeave")'

# Get FAPROTAX
wget https://pages.uoregon.edu/slouca/LoucaLab/archive/FAPROTAX/SECTION_Download/MODULE_Downloads/CLASS_Latest%20release/UNIT_FAPROTAX_1.2.6/FAPROTAX_1.2.6.zip
unzip FAPROTAX_1.2.6.zip
rm FAPROTAX_1.2.6.zip

# Install BugBase
git clone https://github.com/knights-lab/BugBase.git
echo "export BUGBASE_PATH='/$PWD/microbetag_env/BugBase'" >> ~/.bashrc
echo "export PATH='$BUGBASE_PATH/bin:$PATH'" >> ~/.bashrc
source ~/.bashrc


# BugBase R libraries
mkdir R_libs
Rscript -e 'install.packages("dplyr", repos="https://cran.rstudio.com", lib="R_libs")' 
Rscript -e 'install.packages("RColorBrewer2", repos="https://cran.rstudio.com", lib="R_libs")' 
Rscript -e 'install.packages("beeswarm", repos="https://cran.rstudio.com", lib="R_libs")' 
Rscript -e 'install.packages("reshape2", repos="https://cran.rstudio.com", lib="R_libs")' 
Rscript -e 'install.packages("plyr", repos="https://cran.rstudio.com", lib="R_libs")' 
Rscript -e 'install.packages("gridExtra", repos="https://cran.rstudio.com", lib="R_libs")' 
Rscript -e 'install.packages("RJSONIO", repos="https://cran.rstudio.com", lib="R_libs")' 
Rscript -e 'install.packages("digest", repos="https://cran.rstudio.com", lib="R_libs")' 
Rscript -e 'install.packages("optparse", repos="https://cran.rstudio.com", lib="R_libs")' 
Rscript -e 'install.packages("Matrix", repos="https://cran.rstudio.com", lib="R_libs")' 
Rscript -e 'install.packages("labeling", repos="https://cran.rstudio.com", lib="R_libs")' 
Rscript -e 'install.packages("ggplot2", repos="https://cran.rstudio.com", lib="R_libs")'


# FAPROTAX Python libraries
conda create -n microbetag
conda activate microbetag
python -m pip install "numpy<1.24" pytest-shutil biom-format pandas


