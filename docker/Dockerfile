FROM r-base:latest

MAINTAINER VP Nagraj "vpnagraj@virginia.edu"

# Install dependencies and Download and install shiny server
RUN apt-get update && apt-get install -y -t unstable \
    sudo \
    gdebi-core \
    pandoc \
    pandoc-citeproc \
    libsodium-dev \
    libssl-dev \
    libv8-dev \
    libxml2-dev \
    r-cran-rmysql \
    libcurl4-gnutls-dev \
    libcairo2-dev/unstable \
    git \
    libxt-dev && \
    wget --no-verbose https://s3.amazonaws.com/rstudio-shiny-server-os-build/ubuntu-12.04/x86_64/VERSION -O "version.txt" && \
    VERSION=$(cat version.txt)  && \
    wget --no-verbose "https://s3.amazonaws.com/rstudio-shiny-server-os-build/ubuntu-12.04/x86_64/shiny-server-$VERSION-amd64.deb" -O ss-latest.deb && \
    gdebi -n ss-latest.deb && \
    rm -f version.txt ss-latest.deb && \
    R -e "install.packages(c('rmarkdown','devtools','binom','dfcrm','nnet','Iso','V8','pocrm'), repos='https://cran.rstudio.com/')" && \
    rm -rf /var/lib/apt/lists/*

# install explicit versions of certain pkgs
RUN R -e "devtools::install_version('shiny', version = '1.3.2', upgrade = FALSE)"
RUN R -e "devtools::install_version('shinyjs', version = '1.0', upgrade = FALSE)"
RUN R -e "devtools::install_version('shinyjqui', version = '0.3.2', upgrade = FALSE)"
RUN R -e "devtools::install_version('shinyBS', version = '0.61', upgrade = FALSE)"
RUN R -e "devtools::install_version('rhandsontable', version = '0.3.7', upgrade = FALSE)"


# move conf and setup files and clear out shiny server demo apps
COPY shiny-server.conf  /etc/shiny-server/shiny-server.conf
COPY shiny-server.sh /usr/bin/shiny-server.sh
RUN rm -rf /srv/shiny-server/*

# get the app code
RUN git clone https://github.com/uvasomrc/uvatrapps.git
RUN cp -rf uvatrapps/apps/*  /srv/shiny-server/.

# get index.html setup
COPY index.R /srv/shiny-server/index.R
RUN Rscript /srv/shiny-server/index.R

# run the server setup script
CMD ["/usr/bin/shiny-server.sh"]
