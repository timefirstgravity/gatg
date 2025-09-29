# Claude Code Instructions for GATG Project

## Primary Purpose

This project has two complementary scientific goals:

### Goal 1: Mathematical Equivalence Verification
Perform direct symbolic computation with SageMath to prove that
**Gravity as Temporal Geometry (GATG)** (lapse-first General Relativity)
is mathematically equivalent to standard GR.

**Mathematical Target:** Prove that two different formulations of General Relativity
yield identical physical results:

1. **Standard GR**: Start with 4D spacetime metric perturbations
2. **Lapse-First GR (GATG)**: Start with temporal potential Φ (lapse function)

**Key Insight:** The Schwarzschild solution, which requires solving 10 complex
coupled Einstein equations in standard GR, reduces to a single simple ODE
when approached through the lapse-first GATG formulation: `r A'(r) + A(r) - 1 = 0`

### Goal 2: Theoretical Exploration and Experimental Support
Explore the theoretical implications of the lapse-first reformulation and
develop experimental protocols to test predictions unique to this perspective.

**Research Targets:**
- Quantum origin of classical spacetime from fixed-point emergence
- Gravitational dephasing signatures in precision atomic clocks
- Quantum commutator witness protocols for temporal geometry
- Novel physical insights from the temporal potential perspective
- Observable experimental predictions for laboratory and astronomical tests

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

### **RULE 3: ABSOLUTELY NO SILENT FALLBACKS OR TRY/EXCEPT HIDING**

**CRITICAL RULE:** Computations must work correctly or fail loudly.

**ABSOLUTELY PROHIBITED:**

- Silent fallbacks (try/except with alternative computation)
- Returning placeholder strings when computation fails
- Catching exceptions and continuing with different logic
- Any form of "graceful degradation" that hides errors

**REQUIRED:**

- All computations must succeed or raise explicit errors
- Failed computations must halt execution with clear error messages
- If syntax is wrong, fix the syntax - don't work around it
- If mathematics fails, identify and fix the mathematical error

**EXAMPLES:**

- ❌ BAD: `try: result = complex_calculation() except: result = "simplified_result"`
- ❌ BAD: `if computation_fails: return "fallback_value"`
- ✅ GOOD: `result = complex_calculation()  # Let it fail if wrong`
- ✅ GOOD: `if precondition_not_met: raise ValueError("Clear error message")`

**RATIONALE:** Silent fallbacks hide bugs, produce incorrect scientific results,
and compromise the mathematical rigor of the verification. Every computation
must be mathematically exact or the program must stop with a clear error.

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
