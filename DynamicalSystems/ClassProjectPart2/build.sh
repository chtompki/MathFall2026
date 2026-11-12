#!/bin/bash

echo ""
echo "INFO: Build Python using a virtual environment"
echo ""
rm -rf .venv
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

echo ""
echo "INFO: Run Python programs chaos_game.py and devil_staircase.py"
echo ""
python chaos_game.py
python devil_staircase.py

echo ""
echo "INFO: Clean Up Python virtual environment"
echo ""
rm -rf .venv


echo ""
echo "INFO: pdflatex build classProject2.tex twice so figure references and bibliography references are built"
echo ""
pdflatex classProjectPart2.tex
pdflatex classProjectPart2.tex