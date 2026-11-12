#!/bin/bash

function build_python_venv() {
    echo ""
    echo "INFO: Build Python using a virtual environment"
    echo ""
    rm -rf .venv
    python -m venv .venv
    source .venv/bin/activate
    pip install -r requirements.txt
}

function run_python() {
    echo ""
    echo "INFO: Run Python programs chaos_game.py and devil_staircase.py"
    echo ""
    python chaos_game.py
    python devil_staircase.py
}

function clean_python() {
    echo ""
    echo "INFO: Clean Up Python virtual environment"
    echo ""
    rm -rf .venv
}

function build_latex() {
    echo ""
    echo "INFO: pdflatex build classProject2.tex twice so figure references and bibliography references are built"
    echo ""
    pdflatex classProjectPart2.tex
    pdflatex classProjectPart2.tex
}

function clean() {
    echo ""
    echo "INFO: deleting extra LaTeX output (.aux, .log, .out) and .venv if present"
    rm -f classProjectPart2.aux classProjectPart2.log classProjectPart2.out texput.log
    rm -rf .venv
}

# Mode 1: default build with no command line args, build everything except for clean as declared
if [ $# -eq 0 ]; then
    build_python_venv
    run_python
    clean_python
    build_latex
fi
# Mode 2: Run task directly from command line (e.g., ./build.sh compile)
for arg in "$@"; do
    case "$arg" in
        build_python_venv)
            build_python_venv
            ;;
        run_python)
            run_python
            ;;
        clean_python)
            clean_python
            ;;
        build_latex)
            build_latex
            ;;
        clean)
            clean
            ;;
        *)
            echo "Error: Unknown task '$arg'"
            exit 1
            ;;
    esac
done