@echo off
echo ======================================================================
echo   Sour Rice Microbiome-Nutrition Reproducibility Pipeline (Windows)
echo ======================================================================
echo.
python --version >nul 2>&1
if %errorlevel% neq 0 (
    echo [ERROR] Python is not installed or not in your PATH.
    echo Please install Python 3.9+ from https://www.python.org/
    pause
    exit /b 1
)

echo [1/3] Setting up Python virtual environment...
if not exist "venv" (
    python -m venv venv
)
call venv\Scripts\activate

echo [2/3] Installing verified dependencies...
pip install -q -r requirements.txt

echo [3/3] Executing complete analysis pipeline...
python run_analysis.py

echo.
echo ======================================================================
echo   Pipeline finished! Figures and tables are in the ./results directory.
echo ======================================================================
pause
