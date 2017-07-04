ifeq (,$(shell sh -c 'cygpath --version 2> /dev/null'))
  # Unix
  pwd := $$(pwd)
  translate = $1
else
  # Windows mit MSys2/Cygwin
  pwd := $$(cygpath -m "$$(pwd)")
  translate = $(shell echo '$1' | sed 's/:/;/g')
endif

all: build/thesis.pdf

# hier Python-Skripte:

build/Hubb_Ham_Zeit_gzplot.pdf: Hubb_Ham_Zeit_gzplot.py matplotlibrc header-matplotlib.tex Stationäre_Systeme/Hubb_Ham/Hubb_Zust.txt Stationäre_Systeme/Hubb_Ham/Hubb_Ham_j.txt Stationäre_Systeme/Hubb_Ham/Hubb_Ham_d.txt| build
	TEXINPUTS="$(call translate,$(pwd):)" python Hubb_Ham_Zeit_gzplot.py

build/Hubb_gzv_plot.pdf: Hubb_gz.py matplotlibrc header-matplotlib.tex Stationäre_Systeme/Hubb_Eig_Ergebn/Hubb_gzv.txt Stationäre_Systeme/Hubb_Eig_Ergebn/Hubb_gze.txt | build
	TEXINPUTS="$(call translate,$(pwd):)" python Hubb_gz.py

build/Hubb_gze_plot.pdf: Hubb_gz.py matplotlibrc header-matplotlib.tex Stationäre_Systeme/Hubb_Eig_Ergebn/Hubb_gzv.txt Stationäre_Systeme/Hubb_Eig_Ergebn/Hubb_gze.txt | build
	TEXINPUTS="$(call translate,$(pwd):)" python Hubb_gz.py

build/Hubb_Grenz_Plot.pdf: Hubb_Grenz.py matplotlibrc header-matplotlib.tex Stationäre_Systeme/Hubb_Grenz/Hubb_Grenz.txt Stationäre_Systeme/Hubb_Grenz/Hubb_Grenz_lin.txt| build
	TEXINPUTS="$(call translate,$(pwd):)" python Hubb_Grenz.py

build/Hubb_diff.pdf: Hubb_diff.py matplotlibrc header-matplotlib.tex Stationäre_Systeme/Hubb_Eig_Ergebn/Hubb_anr_diff.txt | build
	TEXINPUTS="$(call translate,$(pwd):)" python Hubb_diff.py

# hier weitere Abhängigkeiten für build/thesis.pdf deklarieren:
build/thesis.pdf: build/Hubb_Ham_Zeit_gzplot.pdf build/Hubb_gzv_plot.pdf build/Hubb_gze_plot.pdf build/Hubb_Grenz_Plot.pdf build/Hubb_diff.pdf

build/thesis.pdf: FORCE | build
	  TEXINPUTS="$(call translate,build:)" \
	  BIBINPUTS=build: \
	  max_print_line=1048576 \
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=nonstopmode \
	  --halt-on-error \
	thesis.tex

build:
	mkdir -p build

clean:
	rm -rf build

FORCE:

.PHONY: all clean
