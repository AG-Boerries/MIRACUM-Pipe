#!/usr/bin/env bash

cd ${ana}

${Rscript} ${ana}/Main.R ${case} ${num} ${file_a} ${file_b} ${mtb} ${RscriptPath} ${DatabasePath}

${Rscript} -e "library(knitr); knit('Report.Rnw')"
pdflatex -interaction=nonstopmode Report.tex
pdflatex -interaction=nonstopmode Report.tex