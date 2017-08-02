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

# Ein_Elektron_System:
#build/Ham_Strom.pdf: Ham_Zeit.py matplotlibrc header-matplotlib.tex | build
#  TEXINPUTS="$(call translate,$(pwd):)" python Ham_Zeit.py

#build/RungeKuttaVgl.pdf: RungeKutta.py matplotlibrc header-matplotlib.tex Hubb_Ham_Zeit_txt/Runge_Kutta_vgl.txt | build
#	TEXINPUTS="$(call translate,$(pwd):)" python RungeKutta.py

#build/Hubb_Ham_Zeit_Besetzung.pdf: Hubb_Ham_Zeit_Besetzung.py matplotlibrc header-matplotlib.tex Stationäre_Systeme/Hubb_Ham/Hubb_Zust.txt Stationäre_Systeme/Hubb_Ham/Hubb_Besetzung.txt Hubb_Ham_Zeit_txt/Hubb_Ham_Zeit_Lsg.txt Hubb_Ham_Zeit_txt/Hubb_Ham_Zeit_Linspace.txt | build
#	TEXINPUTS="$(call translate,$(pwd):)" python Hubb_Ham_Zeit_Besetzung.py

#build/Hubb_Strom.pdf: Hubb_Strom.py matplotlibrc header-matplotlib.tex Stationäre_Systeme/Hubb_Ham/Hubb_Zust.txt Hubb_Ham_Zeit_txt/Hubb_Ham_Zeit_Lsg.txt Hubb_Ham_Zeit_txt/Hubb_Ham_Zeit_Linspace.txt | build
#	TEXINPUTS="$(call translate,$(pwd):)" python Hubb_Strom.py

build/Ham2e_Strom.pdf: Ham2e_Strom.py matplotlibrc header-matplotlib.tex | build
	TEXINPUTS="$(call translate,$(pwd):)" python Ham2e_Strom.py

#build/Hubb_Ham_Zeit_Mittelung.pdf: Hubb_Ham_Zeit_Mittelung.py matplotlibrc header-matplotlib.tex Stationäre_Systeme/Hubb_Ham/Hubb_Zust.txt Stationäre_Systeme/Hubb_Ham/Hubb_Ham_j.txt Stationäre_Systeme/Hubb_Ham/Hubb_Ham_d.txt Stationäre_Systeme/Hubb_Ham/Hubb_Besetzung.txt| build
#	 TEXINPUTS="$(call translate,$(pwd):)" python Hubb_Ham_Zeit_Mittelung.py

build/Hubb_eplot.pdf: Hubb_eplot.py matplotlibrc header-matplotlib.tex Stationäre_Systeme/Hubb_Eig_Ergebn/Hubb_eplot.txt | build
	TEXINPUTS="$(call translate,$(pwd):)" python Hubb_eplot.py

build/Hubb_ediff_plot.pdf: Hubb_ediff.py matplotlibrc header-matplotlib.tex Stationäre_Systeme/Hubb_Eig_Ergebn/Hubb_ediff.txt| build
	TEXINPUTS="$(call translate,$(pwd):)" python Hubb_ediff.py

build/Hubb_gzv_plot.pdf: Hubb_gz.py matplotlibrc header-matplotlib.tex Stationäre_Systeme/Hubb_Eig_Ergebn/Hubb_gzv.txt Stationäre_Systeme/Hubb_Eig_Ergebn/Hubb_eplot.txt | build
	TEXINPUTS="$(call translate,$(pwd):)" python Hubb_gz.py

build/Hubb_gze_plot.pdf: Hubb_gz.py matplotlibrc header-matplotlib.tex Stationäre_Systeme/Hubb_Eig_Ergebn/Hubb_gzv.txt Stationäre_Systeme/Hubb_Eig_Ergebn/Hubb_eplot.txt | build
	TEXINPUTS="$(call translate,$(pwd):)" python Hubb_gz.py

build/Hubb_Grenz_Plot.pdf: Hubb_Grenz.py matplotlibrc header-matplotlib.tex Stationäre_Systeme/Hubb_Grenz/Hubb_Grenz.txt Stationäre_Systeme/Hubb_Grenz/Hubb_Grenz_lin.txt| build
	TEXINPUTS="$(call translate,$(pwd):)" python Hubb_Grenz.py

build/Hubb_diff.pdf: Hubb_diff.py matplotlibrc header-matplotlib.tex Stationäre_Systeme/Hubb_Eig_Ergebn/Hubb_anr_diff.txt | build
	TEXINPUTS="$(call translate,$(pwd):)" python Hubb_diff.py

# hier weitere Abhängigkeiten für build/thesis.pdf deklarieren:
build/thesis.pdf:  build/Ham2e_Strom.pdf build/Hubb_eplot.pdf build/Hubb_ediff_plot.pdf build/Hubb_gzv_plot.pdf build/Hubb_gze_plot.pdf build/Hubb_Grenz_Plot.pdf build/Hubb_diff.pdf
#build/Hubb_Ham_Zeit_Besetzung.pdf build/Hubb_Strom.pdf
#build/RungeKuttaVgl.pdf build/Ham_Strom.pdf build/Hubb_Ham_Zeit_Mittelung.pdf
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
