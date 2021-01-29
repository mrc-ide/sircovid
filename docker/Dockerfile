FROM rocker/r-ver:4.0.3

RUN apt-get update &&  apt-get install -y --no-install-recommends \
        gnupg2 \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        valgrind \
        wget \
        zlib1g-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Without this, we are unable to pick up more recent packages
COPY docker/Rprofile.site /usr/local/lib/R/etc/Rprofile.site

RUN install2.r --error \
        cpp11 \
        remotes \
        socialmixr

RUN Rscript -e 'remotes::install_github(c("mrc-ide/dust", "mrc-ide/mcstate"))'

RUN Rscript -e 'remotes::install_github("mrc-ide/dust")'

COPY . /src
RUN R CMD INSTALL --install-tests /src && rm -rf /src
