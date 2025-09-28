# Claude Code Instructions for GATG Project

## Primary Purpose

This project performs direct symbolic computation with SageMath
to prove that **Gravity as Temporal Geometry (GATG)** (lapse-first General Relativity)
is mathematically equivalent to standard GR.

**Mathematical Goal:** Prove that two different formulations of General Relativity
yield identical physical results:

1. **Standard GR**: Start with 4D spacetime metric perturbations
2. **Lapse-First GR (GATG)**: Start with temporal potential Φ (lapse function)

**Key Insight:** The Schwarzschild solution, which requires solving 10 complex
coupled Einstein equations in standard GR, reduces to a single simple ODE
when approached through the lapse-first GATG formulation: `r A'(r) + A(r) - 1 = 0`

**Verification Target:** Prove mathematically that both formulations are equivalent.

**Research Target:** Explore the implications of time first gravity

## Python/Sage Execution Instructions

**ALWAYS use `sage -python` to run Python scripts in this project**

- Do NOT use plain `python` command - it will fail
- The project uses SageMath's Python environment
- All scripts should be run as: `sage -python script.py`

## Critical Scientific Computing Rules

### **RULE 1: NEVER USE FAKE DATA IN ANY MODULE OF THIS PROGRAM**

This is a rigorous mathematical verification program.
Every computation **MUST** use SageMath symbolic computation
Every value, data point, every number, every plot must use ONLY
actual computed results from the mathematical computation.

**Absolutely prohibited in ALL modules:**

- Random number generation for "demonstration"
- Synthetic data for "visualization purposes"
- Made-up values for "testing"
- Fabricated examples for "illustration"

**Required for ALL modules:**

- Only real computed results from actual calculations
- Actual symbolic computation outputs
- True mathematical verification values
- Genuine scientific data

### **RULE 2: ALL FUNCTIONS MUST PERFORM MATHEMATICAL COMPUTATION**

**ABSOLUTELY PROHIBITED:**

- "Textbook mode" functions that explain concepts without calculating
- Hardcoding result values

**REQUIRED FOR ALL FUNCTIONS:**

- Take mathematical inputs (metrics, coordinates, fields, etc.)
- Perform actual symbolic or numerical computation
- Return computed mathematical results (expressions, matrices, tensors, etc.)
- Enable verification between Standard GR and Lapse-First GR formulations

**EXAMPLES:**

- ❌ BAD: `function() { return "Einstein tensor describes curvature" }`
- ✅ GOOD: `function(metric, coordinates) { return computed_einstein_tensor }`

**VERIFICATION PURPOSE:** Every function must contribute to the rigorous mathematical proof
that Standard GR and Lapse-First GR are equivalent.
Documentation belongs in comments and docstrings, NOT in function return values.

The integrity of the scientific computation is paramount across the entire GATG project.

## Debug Script Management Rule

**ALWAYS CLEAN UP DEBUG SCRIPTS: DEBUG → EXTRACT → INTEGRATE → DELETE**

When creating debugging/testing scripts during development,
you MUST follow this process:

1. **DEBUG**: Create temporary scripts to investigate issues or test improvements
2. **EXTRACT**: All computational improvements, precision enhancements, and fixes
3. **INTEGRATE**: Apply ALL improvements to the main codebase
4. **DELETE**: Remove the debug script to prevent codebase clutter

**The main codebase should contain ONLY:**

- Core functionality modules
- Main verification scripts
- Essential physics modules

**NO debug/test/investigate scripts should remain**
after improvements are extracted and integrated.

This systematic approach prevents accumulation of
debugging scripts and ensures all computational improvements
are properly consolidated into the core scientific codebase.
