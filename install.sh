#!/bin/bash
# GATG Project Installation Script

echo "Installing GATG Project Dependencies..."
echo "======================================"

# Function to check if SageMath is installed
check_sage() {
    if command -v sage &> /dev/null; then
        echo "✓ SageMath found at: $(which sage)"
        return 0
    else
        return 1
    fi
}

# Function to install packages into SageMath environment
install_with_sage() {
    echo ""
    echo "Installing Python packages into SageMath environment..."
    echo "--------------------------------------------------------"

    # Use sage -pip to install into SageMath's Python environment
    sage -pip install --upgrade pip
    sage -pip install -r requirements.txt

    if [ $? -eq 0 ]; then
        echo "✓ Successfully installed all dependencies into SageMath environment"
    else
        echo "⚠ Some packages may have failed to install. Please check the output above."
    fi
}

# Main installation logic
if check_sage; then
    # SageMath is installed, use it
    install_with_sage

elif command -v conda &> /dev/null; then
    # Conda is available but SageMath not installed
    echo "✓ Conda found but SageMath not installed"
    echo ""
    echo "Installing SageMath via conda..."
    conda install -c conda-forge sage -y

    if check_sage; then
        install_with_sage
    else
        echo "✗ Failed to install SageMath via conda"
        exit 1
    fi

else
    # Neither SageMath nor conda found
    echo "✗ SageMath not found. Please install it first:"
    echo ""
    echo "Option 1: Install via Conda (recommended)"
    echo "  1. Install Miniconda from: https://docs.conda.io/en/latest/miniconda.html"
    echo "  2. Run: conda install -c conda-forge sage"
    echo "  3. Run this script again"
    echo ""
    echo "Option 2: Install via system package manager"
    echo "  macOS:  brew install sage"
    echo "  Ubuntu: sudo apt install sagemath"
    echo "  Fedora: sudo dnf install sagemath"
    echo ""
    echo "Option 3: Download binary from https://www.sagemath.org/download.html"
    echo ""
    exit 1
fi

echo ""
echo "======================================"
echo "Installation complete!"
echo ""
echo "Usage Instructions:"
echo "  • Run verification scripts:  sage -python verification.py"
echo "  • Run tests:                 sage -python -m pytest"
echo "  • Interactive SageMath:      sage"
echo ""
echo "Available verification modules:"
echo "  schwarzschild/           - Schwarzschild solution equivalence"
echo "  kerr/                    - Kerr rotating black hole equivalence"
echo "  linearized_gravity/      - Linearized GR wave equation equivalence"
echo "  cosmology/               - FLRW cosmological equivalence"
echo "  gravitoelectromagnetism/ - GEM weak field equivalence"
echo "  experimental_predictions/ - Observable effects equivalence"
echo "  coordinate_transforms/   - Coordinate transformation equivalence"
echo "  flux_law/                - Conservation law equivalence"
echo ""
echo "Example commands:"
echo "  cd schwarzschild && sage -python verification.py"
echo "  cd kerr && sage -python verification.py"
echo "  cd flux_law && sage -python verification.py"
echo "  cd cosmology && sage -python verification.py"
echo "  sage -python -m pytest core/tests/ -v"
echo ""
echo "Expected output for all verifications:"
echo "  ✓ EQUIVALENCE VERIFIED: Standard GR ≡ Lapse-First GR"
echo ""
echo "Note: Always use 'sage -python' instead of plain 'python' to ensure"
echo "      the SageMath environment is properly loaded for symbolic computation."