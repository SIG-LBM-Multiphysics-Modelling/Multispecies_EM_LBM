# Multispecies EM-LBM Solver

This repository contains the simulation code used in the accompanying paper.

## Paper

Status: Submitted, under review.

---

## Repository Structure

The repository is organised into two main groups of simulation cases corresponding to the paper sections.

### `multispecies/` (Section 3 of the paper)
Contains multispecies fluid simulations with:
- drag forcing
- BVE coupling (where applicable)

### `multispecies_em/` (Section 4 of the paper)
Contains multispecies simulations with:
- drag forcing
- electromagnetic coupling via Poisson solves

Each folder contains:
- `main.cpp` (simulation driver)
- `VTIWriter.h` (output utility specific to the case)

---

## Model Configuration

This repository provides the base implementation used in the paper.

All results reported in the paper are obtained from the same underlying codebase by:
- enabling or disabling specific coupling terms (e.g. drag, BVE, electromagnetic coupling)
- selecting boundary condition options
- modifying parameters within the existing solver structure
- calling or omitting specific model components already implemented in the code

No structural changes to the solver are required between cases; variations correspond to controlled changes in model configuration within each case implementation.

---

## Running the Cases

Each case is self-contained.

To run a simulation:

### Compile

```bash
g++ main.cpp -O2 -std=c++17 -o sim
```

### Run

Run the executable depending on your terminal:

- Linux / WSL / MSYS2 (Git Bash): `./sim`
- Windows PowerShell / VS Code terminal: `.\sim.exe`
- Windows Command Prompt: `sim.exe`
