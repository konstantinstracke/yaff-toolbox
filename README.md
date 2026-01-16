# Toolbox for YAFF

This repository contains a small collection of utility scripts for file handling, preprocessing, and analysis in workflows based on YAFF.

The focus is mainly on conversion scripts that help bridge different file formats commonly used in atomistic simulations (e.g. MFPX, XYZ, CHK), along with lightweight analysis helpers.

---

## Repository Structure

.
├── 1-conversion/
│   └── Scripts for converting and preprocessing input/output files
│
├── 2-analysis/
│   └── Simple analysis tools for YAFF-based simulations
│
└── README.md

---

## Features

- Conversion of structure files to YAFF-compatible formats
- Preprocessing helpers (e.g. cleaning .mfpx files)
- Lightweight analysis scripts for simulation output
- Designed to be simple, modular, and easy to adapt

---

## Citation

If you use YAFF in connection with this toolbox, please cite:

Verstraelen, T.; Vanduyfhuys, L.; Vandenbrande, S.; Rogge, S.  
YAFF, Yet Another Force Field.  
Available at: http://molmod.ugent.be/software  
2022. Accessed on 25-11-2024.  
(No corresponding bibliographic record is available.)

---

## Notes

- These scripts are not a polished package, but a practical toolbox developed for day-to-day research use.
- File formats and conventions are assumed to match the author’s simulation setup.
- Use and modify at your own discretion.
