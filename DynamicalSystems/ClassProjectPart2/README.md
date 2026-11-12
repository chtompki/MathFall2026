# MATH 535: Class Project, Part 2

### Build requirements
The following must be in your path. We rely upon `python` as a symlink to `python3` (native python install) or
`/opt/homebrew/bin/python3.14`. Furthermore, you need `pip` simlinked to the associated `pip3` or `pip3.14`.
Next we rely upon `pdflatex` (use whatever the latest version is for compilation). The libraries needed in your path
follow as: 

* Python - 3.14
* LaTeX - latest version available.

The report reproduces Barnsley and Vince's Example 1 and Bak's Figure 5(c).
The original project overview is included as the abstract.

From this directory, install NumPy and Matplotlib in your Python environment, then run:

```sh
python -m venv .venv
source .venv/bin/activete
python -m pip install -r requirements.txt
python chaos_game.py
python devil_staircase.py
rm -rf .venv
pdflatex classProjectPart2.tex
pdflatex classProjectPart2.tex
```

or run `./build.sh` as it contains this series of commands.

The scripts write CSV tables to `data/` and Matplotlib previews to `plots/`.
With only NumPy installed, pass `--data-only` to regenerate the CSV tables.
The report draws the Python-based figures directly from those tables with TikZ/PGFPlots;
it does not embed the Matplotlib previews or require Python during compilation.
The CSV files are included so the report can be compiled immediately with a
TeX distribution that provides PGFPlots 1.18.

### MATLAB figures

Run the following from the MATLAB Command Window with the project directory as the current folder:

~~~matlab
circle_map_dimension
~~~

The script saves three 300 dpi PNG files in the project's own plots folder:
matlab_circle_map.png, matlab_staircase.png, and matlab_dimension_fit.png.
Recompile classProjectPart2.tex twice to include them in the appendix.
The paper also compiles before any MATLAB PNGs have been generated.
