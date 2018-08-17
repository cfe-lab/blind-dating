# R base
Bootstrap: docker
From: r-base

%labels
AUTHOR Bradley R. Jones, BC CfE in HIV/AIDS

%files
 src/*.R /opt/blind-dating/
 src/blind-dating /opt/blind-dating/
 src/blind-dating2 /opt/blind-dating/

%environment
ENV BDSRC=/opt/blind-dating

%post

# prerequistes
 apt-get update --fix-missing && apt-get install -y \
  zlib1g \
  unzip \
  wget \
  make \
  git \
  libssl-dev \
  libcurl4-openssl-dev \
  libxml2-dev \
  && rm -rf /var/lib/apt/lists/*

# RAxML
 wget -q -O raxml.zip https://github.com/stamatak/standard-RAxML/archive/v8.2.12.zip && \
  unzip raxml.zip -d /tmp && \
  rm raxml.zip
 cd /tmp/standard-RAxML-8.2.12
 make -f Makefile.PTHREADS.gcc && \
  rm *.o && \
  ln -s /tmp/standard-RAxML-8.2.12/raxmlHPC-PTHREADS /opt/raxml

# CRAN R packages and ggtree
 R --vanilla --slave -e 'local({r <- getOption("repos"); r["CRAN"] <- "http://cran.stat.sfu.ca"; options(repos=r)}); install.packages(c("ape", "chemCal", "ggplot2", "optparse", "phylobase", "seqinr", "weights")); update.packages(ask=FALSE); source("https://bioconductor.org/biocLite.R"); biocLite("ggtree")'

# node.dating
 cd /tmp
 git clone https://github.com/brj1/node.dating.git
 cd /tmp/node.dating
 git checkout random && \
  ln -s /tmp/node.dating/R/node.dating.R /opt/node.dating.R 

# scripts
 chmod +x /opt/blind-dating/blind-dating
 chmod +x /opt/blind-dating/blind-dating2

