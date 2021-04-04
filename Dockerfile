FROM rocker/r-base:4.0.3
RUN install2.r Rcpp igraph parallel getopt
WORKDIR /scripts/sim
COPY ./scr/*{R,cpp} /scripts/sim/
WORKDIR /scripts/plots/
COPY Figures/*R /scripts/plot/
