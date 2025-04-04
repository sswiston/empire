# BASE IMAGE
FROM alpine:3.19.1
ENV HOME="/"

# BASE SOFTWARE
WORKDIR "/"
RUN apk update && \
    apk upgrade && \
    apk add bash git unzip vim gcc g++ cmake protoc make automake autoconf curl curl-dev gfortran build-base

# R
WORKDIR "/"
RUN apk add --no-cache bzip2-dev ca-certificates pcre2-dev readline-dev xz-dev tzdata
RUN wget https://cran.r-project.org/src/base/R-4/R-4.4.0.tar.gz && \
    tar xfz R-4.4.0.tar.gz; rm R-4.4.0.tar.gz
WORKDIR "/R-4.4.0"
RUN ./configure --with-x=no --enable-R-shlib --enable-memory-profiling && \
    make -s ; make install
ENV PATH /R-4.4.0/bin:$PATH
ENV R_HOME /R-4.4.0
ENV R_LIBS_USER /R-4.4.0/library
ENV LD_LIBRARY_PATH /R-4.4.0/lib:$LD_LIBRARY_PATH

# R PACKAGES
# "Rcpp" "devtools" "R6" "ape" "phytools" "distributions3" "sf" "hms" "condMVNorm" "mcmcse" "rase" "spatstat"
WORKDIR "/"
RUN apk add linux-headers libgit2-dev freetype-dev libxml2-dev fontconfig-dev harfbuzz-dev fribidi-dev libpng-dev tiff-dev jpeg-dev imagemagick-dev udunits-dev gdal-dev font-terminus font-noto font-noto-extra font-arabic-misc gmp-dev mpfr-dev gsl-dev proj-dev geos-dev texlive libpq-dev openssl-dev unixodbc-dev fftw-dev

RUN mkdir /.R
RUN touch /.R/Makevars
RUN echo "PKG_CXXFLAGS += -D__MUSL__" >> /.R/Makevars
RUN echo "CXXFLAGS += -D__MUSL__" >> /.R/Makevars
RUN echo "CXX1XFLAGS += -D__MUSL__" >> /.R/Makevars
RUN cat /.R/Makevars
RUN R -e "install.packages('Rcpp',dependencies=TRUE,repos='http://cran.rstudio.com/')" && \
    R -e "install.packages('devtools',dependencies=TRUE,repos='http://cran.rstudio.com/')" && \
    R -e "install.packages('R6',dependencies=TRUE,repos='http://cran.rstudio.com/')" && \
    R -e "install.packages('ape',dependencies=TRUE,repos='http://cran.rstudio.com/')" && \
    R -e "install.packages('phytools',dependencies=TRUE,repos='http://cran.rstudio.com/')" && \
    R -e "install.packages('distributions3',dependencies=TRUE,repos='http://cran.rstudio.com/')" && \
    R -e "install.packages('sf',dependencies=TRUE,repos='http://cran.rstudio.com/')" && \
    R -e "install.packages('hms',dependencies=TRUE,repos='http://cran.rstudio.com/')" && \
    R -e "install.packages('condMVNorm',dependencies=TRUE,repos='http://cran.rstudio.com/')" && \
    R -e "install.packages('mcmcse',dependencies=TRUE,repos='http://cran.rstudio.com/')" && \
    R -e "library(devtools);install_github(\"ignacioq/rase\")" && \
    R -e "library(devtools);install_github(\"spatstat/spatstat\")"

# FINISHING
WORKDIR "/"
COPY dockerfile_ellipses.txt /dockerfile_ellipses.txt
COPY test.sh /test.sh
CMD sleep infinity

# Command for building Docker image:
# cp Dockerfile dockerfile_ellipses.txt; docker build --no-cache -t sswiston/ellipses:2 . --progress=plain > ellipses.log 2>&1