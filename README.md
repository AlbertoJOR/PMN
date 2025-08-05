# Pre-Metastatic Niche Model
An **agent-based simulation** of the **pre-metastatic niche (PMN)** formation process using [Mesa](https://mesa.readthedocs.io/en/stable/), a Python framework for agent-based modeling. This project models the interactions between immune cells, exosomes, and signaling molecules involved in preparing the microenvironment for metastasis.

---

## Purpose

To explore and visualize how different cellular and molecular components contribute to the formation of the pre-metastatic niche, a key step in cancer progression.

---

## Modeled Entities

- **Immune Cells**: CD8, CD4, NK, Treg, M1, M2, N1, N2, BMDC, etc.
- **Tumor Cells**: CTC (Circulating Tumor Cells)
- **Microenvironment Cells**: Stromal cells, fibroblasts
- **Exosomes**: PD-L1, MDE, MPCDE, ODEV, MCDEV
- **Molecules**: IL-12, VEGF, IL-6, S100, TGFB, fibronectin

---

## Technologies

- Python 3
- [Mesa](https://github.com/projectmesa/mesa)
- Interactive visualization (`CanvasGrid`, `ChartModule`)
- CSV data export (cell and molecule abundance per step)

---

## Installation

```bash
git clone https://github.com/AlbertoJOR/PMN.git
cd PMN
pip install -r requirements.txt  # Make sure to include this or install Mesa manually
