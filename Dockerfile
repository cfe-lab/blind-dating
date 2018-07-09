# R base
FROM r-base

MAINTAINER Bradley R. Jones, BC CfE in HIV/AIDS

# Prerequistes
RUN apt-get update --fix-missing && apt-get install -y \
  unzip \
  wget \
  make \
  git \
  libssl-dev \
  libcurl4-openssl-dev \
  libxml2-dev \
  && rm -rf /var/lib/apt/lists/*

# RAxML
RUN wget -q -O raxml.zip https://github.com/stamatak/standard-RAxML/archive/v8.2.12.zip && \
  unzip raxml.zip -d /tmp && \
  rm raxml.zip
WORKDIR /tmp/standard-RAxML-8.2.12
RUN make -f Makefile.PTHREADS.gcc && \
  rm *.o && \
  ln -s raxmlHPC-PTHREADS /opt/raxml

# CRAN R packages and ggtree
RUN R --slave -e 'local({r <- getOption("repos"); r["CRAN"] <- "http://cran.stat.sfu.ca"; options(repos=r)}); install.packages(c("ape", "chemCal", "ggplot2", "optparse", "phylobase", "seqinr", "TreeSim", "parallel", "phylobase")); source("https://bioconductor.org/biocLite.R"); biocLite(); biocLite("ggtree")'

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

# scripts
COPY src/*.R /opt/blind-dating/
COPY src/blind-dating /opt/blind-dating/
RUN chmod 777 /opt/blind-dating/blind-dating

# command
CMD ["/opt/blind-dating/blind-dating", "info_csv", "align_fasta", "settings_txt"]
