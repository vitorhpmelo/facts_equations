# FACTS Equations

This repository contains symbolic derivations of current, power, and measurement equations for FACTS devices using SymPy.

The equations are used during the development of my PhD thesis on power system state estimation and are intended to generate analytical expressions and Julia code.

## Devices

- SVC
- TCSC
- UPFC

## Structure

- `src/`: reusable symbolic functions
- `*_hvector.py`: symbolic derivations for each device

## Requirements

```bash
pip install -r requirements.txt
```
