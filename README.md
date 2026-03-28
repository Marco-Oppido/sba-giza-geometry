# SBA-Giza Geometric Scaling Proofbench

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/)

This repository contains the numerical proofbench and the LaTeX manuscript for the paper **"Giza Pyramids & SBA"**. It serves as a falsifiable geometric case study testing the macro-scale limits of the **Specular Bit Architecture (SBA)** framework against the architectural dimensions of the Giza plateau.

⚠️ **Note on Scientific Context:** This repository is an applied geometric test. For the foundational theorems, deterministic update semantics, and cosmological implications of the SBA framework, please refer to the primary theoretical repository: `[https://github.com/Marco-Oppido/Cosmological-SBA-Framework]`

## 📌 Overview

The SBA model predicts that the continuum limit of its discrete, Fibonacci-weighted causal expansion will exhibit specific scaling constants. If the Giza pyramids (Menkaure, Khafre, Khufu) represent a macroscopic architectural reflection of this informational growth, their dimensions should align with these constants.

This repository tests two primary observables using Menkaure as the geometric seed (unit 1):
1. **Base Ratio ($R_B$):** The scaling from Menkaure's base to Khufu's base tested against the target **$\sqrt{5}$**.
2. **Volume Ratio ($R_V$):** The scaling from Menkaure's volume to Khufu's volume tested against the target **$\varphi^5$** (where $\varphi$ is the Golden Ratio).

Current results using standard archaeological working values show residuals of **$\approx 1.54\%$** for the base ratio and **$\approx 1.38\%$** for the volume ratio. 

## 📂 Repository Structure

* `main.tex` - The IEEE-style LaTeX source code of the paper.
* `[sba_pyramid_proofbench.py](https://github.com/Marco-Oppido/sba-giza-geometry/blob/main/sba_pyramid_proofbench.py)` - The Python script containing the dimensional data, target constant calculations, and Monte Carlo uncertainty analysis.
* `figures/` - Directory containing the generated output plots (Endpoint Ratios and Monte Carlo Intervals) required for the LaTeX compilation.

## 🚀 How to Run

### 1. Run the Python Proofbench
Ensure you have Python 3.8+ and the required scientific libraries installed:
```bash
pip install numpy matplotlib
```
Execute the analysis script:
```bash
python sba_pyramid_proofbench.py
```
*This will print the endpoint calculations, run the Monte Carlo uncertainty simulations, and automatically save the required `.png` plots into the `figures/` folder.*

### 2. Compile the Paper
Once the figures are generated, compile the LaTeX document:
```bash
pdflatex main.tex
pdflatex main.tex
```

## 🔬 Falsifiability & Methodology

Following strict scientific methodology, this project enforces a one-way evidentiary flow:
`Formal Layer (Proved SBA Algebra) ➔ Model Layer (Derived Scaling) ➔ Hypothesis Layer`

The Python script explicitly separates data inputs from target constants to prevent target-locking and includes uncertainty propagation to account for measurement variations in historical reconstructions. 

## ✍️ Author
**Marco Oppido**  
Independent Researcher  
*Contact: scaccomarco@gmail.com*

## 📜 License
This project is licensed under the MIT License - see the LICENSE file for details. If you use this code or the associated geometric scaling hypothesis in your own research or media, please cite this repository and the author.
