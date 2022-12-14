
FROM rocker/verse:4.2.1
## tag 1.1 

ARG R_VERSION
ARG BUILD_DATE
ARG CRAN
ENV BUILD_DATE ${BUILD_DATE:-2022-10-26}
ENV R_VERSION=${R_VERSION:-4.2.1} \
    CRAN=${CRAN:-https://cran.rstudio.com} \ 
    TERM=xterm

# Set the locale
#RUN sed -i -e 's/# en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen && \
#    locale-gen
ENV LANG en_US.UTF-8  
ENV LANGUAGE en_US:en  
ENV LC_ALL en_US.UTF-8   


RUN useradd docker \
	&& mkdir /home/docker \
	&& chown docker:docker /home/docker \
	&& addgroup docker staff

RUN  apt-get update \
	&& DEBIAN_FRONTEND="noninteractive" apt-get install -y --no-install-recommends \
  apt-utils \
  gpg gpg-agent \
  python3-pip \
  tmux \
  screen
  
RUN install2.r --error --deps TRUE \
    rstan \
    && rm -rf /tmp/downloaded_packages/ /tmp/*.rds
    
RUN mkdir -p $HOME/.R/ \
    && echo "CXXFLAGS=-O3 -mtune=native -march=native -Wno-unused-variable -Wno-unused-function -flto -ffat-lto-objects  -Wno-unused-local-typedefs \n" >> $HOME/.R/Makevars

 RUN install2.r --error --deps TRUE \
    rstanarm \
    ggmcmc \
    foreach \
    doParallel \
    && rm -rf /tmp/downloaded_packages/ /tmp/*.rds   
    
# fixes a wierd error
RUN install2.r --error Rcpp
RUN R -e 'options(repos="https://cran.rstudio.com/" ) ; install.packages("posterior"); install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")) ); library(cmdstanr) ; install_cmdstan(cores = 4)'




RUN chmod -R 777 /root ; chmod 777 /usr/local/lib/R/site-library /usr/local/lib/R/library
