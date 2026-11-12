#!/usr/bin/env bash
set -e
echo "======================================================================"
echo "  Sour Rice Microbiome-Nutrition Reproducibility Pipeline (Linux/Mac)"
echo "======================================================================"

if ! command -v python3 &> /dev/null; then
    echo "[ERROR] Python 3 could not be found. Please install Python 3.9+."
    exit 1
fi

echo "[1/3] Setting up virtual environment..."
if [ ! -d "venv" ]; then
    python3 -m venv venv
fi
source venv/bin/activate

echo "[2/3] Installing dependencies..."
pip install -q -r requirements.txt

echo "[3/3] Executing analysis pipeline..."
python3 run_analysis.py

echo "======================================================================"
echo "  Pipeline finished! Outputs saved in ./results"
echo "======================================================================"
