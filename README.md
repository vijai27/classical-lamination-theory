# Classical Lamination Theory (CLT) Analysis

MATLAB implementation of Classical Lamination Theory for IM7/HexPly 8552 carbon fibre composite laminates.

## What it does

- Builds the **ABD stiffness matrix** from ply angles and material properties
- Computes **flexural modulus** and **tensile modulus**
- Performs **first ply failure** analysis using the Tsai-Hill criterion
- Runs **progressive ply failure** analysis at 1% applied strain
- Estimates **prepreg material area** required for fabrication

## Material System

| Property | Value |
|----------|-------|
| E₁ | 162.0 GPa |
| E₂ | 8.96 GPa |
| G₁₂ | 4.69 GPa |
| ν₁₂ | 0.316 |
| Xt / Xc | 2.558 / 1.689 GPa |
| Yt / Yc | 0.064 / 0.286 GPa |
| S | 0.0911 GPa |

## Usage

1. Open `Project.m` in MATLAB
2. Edit the `layup` array to define your ply angles (in degrees):
   ```matlab
   layup = [27.5 0 0 0 0 0 0 0 0 0 0 27.5];
   ```
3. Run the script — results are printed to the Command Window

## Output

```
=== CLT ANALYSIS FOR LAYUP ===
Flexural Modulus (Ef): XX.XX GPa  (error vs target: X.X%)
Tensile Modulus (Ex):  XX.XX GPa

Ultimate Tensile Strength: XXX.X MPa

PROGRESSIVE FAILURE ANALYSIS (1% strain):
Iteration 1: X° ply failed
...
Remaining plies: X

Prepreg Area Required: XXXX.X cm²

=== SUMMARY ===
```

## Requirements

- MATLAB R2019b or later (uses `linspace`, `cosd`, `sind`, local functions)
