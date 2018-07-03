# R base
FROM r-base

MAINTAINER Bradley R. Jones, BC CfE in HIV/AIDS

# Prerequistes
RUN apt-get update --qq --fix-missing && apt-get install -qq -y \
  unzip \
  wget \
  make \
  git \
  && rm -rf /var/lib/apt/lists/*

# RAxML
RUN wget -q -O raxml.zip https://github.com/stamatak/standard-RAxML/archive/v8.2.12.zip && \
  unzip raxml.zip -d /tmp && \
  rm raxml.zip
WORKDIR /tmp/standard-RAxML-8.2.12
RUN make -f Makefile.PTHREADS.gcc && \
  rm *.o &&
  ln -s raxmlHPC-PTHREADS /opt/raxml

# CRAN R packages
RUN R --silent --slave -c 'local({r <- getOption("repos") r["CRAN"] <- "http://cran.stat.sfu.ca" options(repos=r)})
  install.packages(c(
    "ape",
    "chemCal",
    "ggplot2",
    "optparse",
    "phylobase",
    "seqinr",
    "TreeSim"
  ));'

# NELSI
WORKDIR /tmp
RUN git clone https://github.com/sebastianduchene/NELSI.git
RUN R CMD INSTALL NELSI

# node.dating
WORKDIR /tmp
RUN git clone https://github.com/brj1/node.dating.git
WORKDIR /tmp/node.dating
RUN git branch random && \
 ln -s R/node.dating.R /opt/node.dating.R

# ggtree
RUN R --silent --slave -c 'source("https://bioconductor.org/biocLite.R")
  biocLite()
  biocLite("ggtree")'

# scripts
COPY *.R /opt/blind-dating/
