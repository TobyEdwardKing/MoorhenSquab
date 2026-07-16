# MoorhenSquab

A WebGPU-accelerated Small-Angle X-ray Scattering (SAXS) profile calculator for macromolecular structures.

## Overview

MoorhenSquab is a browser-based plugin for Moorhen designed to compute scattering profiles directly from atomic coordinates using the Debye equation.

The central aim is to reproduce to an appropriate standard of accuracy a computed SAXS curve that may be produced by PEPSI-SAXS or CRYSOL.

---

## Features

Current implementation includes:

* PDB parsing using Gemmi
* Atomic coordinate extraction
* Cromer–Mann atomic form factors
* Excluded solvent correction
* Hydration shell generation
* WebAssembly backend (C++)
* WebGPU computation of Debye equation
* Excluded volume and hydration shell consideration
* React + TypeScript frontend



---

## Repository Structure

```text
.
├── include/           # Public C++ headers and static resources
│   ├── DebyeCalculator.hpp
│   ├── DebyePipeline.hpp
│   ├── atomExtractor.hpp
│   ├── formFactor.hpp
│   ├── hydrationShell.hpp
│   ├── hydrationShell.cpp
│   ├── scatteringModel.hpp
│   ├── elementsInfo.txt
│   ├── shaders/
│   └── gemmi/
│
├── src/               # C++ implementation
│   ├── main.cpp
│   ├── ApplicationState.cpp
│   ├── GPUContext.cpp
│   ├── DebyeCalculator.cpp
│   ├── DebyePipeline.cpp
│   ├── atomExtractor.cpp
│   ├── formFactor.cpp
│   └── scatteringModel.cpp
│
├── frontend/          # React + TypeScript frontend
│   ├── src/
│   ├── public/
│   ├── package.json
│   └── vite.config.ts
│
├── checkout/          # Third-party libraries
│   └── gemmi/
│
├── CMakeLists.txt
├── Doxyfile
└── README.md
```

## Calculation Pipeline

The SAXS calculation currently follows the workflow:

```text
PDB/CIF
   │
   ▼
AtomExtractor
   │
   ▼
FormFactorTable
   │
   ▼
ScatteringModel
   │
   ├── effective atomic amplitudes
   └── hydration shell
            │
            ▼
DebyeCalculator
   │
   ├── CPU Debye equation
   └── GPU buffers
            │
            ▼
DebyePipeline
            │
            ▼
WebGPU shader
            │
            ▼
Intensity I(q)

I(q)
=
I_atom-atom
+
I_atom-shell
+
I_shell-shell
```

---
## Component summary
| Component         | Purpose                                                 |
| ----------------- | ------------------------------------------------------- |
| `DebyeCalculator` | CPU Debye scattering calculation                        |
| `DebyePipeline`   | WebGPU pipeline interface                               |
| `GPUContext`      | WebGPU device/buffer management                         |
| `scatteringModel` | Atomic amplitudes, excluded solvent and hydration shell |
| `hydrationShell`  | Generates hydration shell scatterers                    |
| `atomExtractor`   | Reads coordinates from Gemmi structures                 |
| `formFactor`      | Atomic form factors and excluded-volume amplitudes      |

## Technologies

* C++20
* WebAssembly (Emscripten)
* WebGPU
* WGSL
* React
* TypeScript
* Vite

---

## Building

### Backend

```mkdir build-wasm
cd build-wasm

emcmake cmake ..
cmake --build . -j

cd ..

cp build-wasm/saxs.wasm frontend/public/
cp build-wasm/saxs.js frontend/src/wasm/
```

### Frontend

```bash
npm install
npm run dev
```

---

## Documentation

API documentation is generated using Doxygen.

```bash
doxygen Doxyfile
```

Generated documentation is written to

```text
docs/html/index.html
```

---

## Current Status

The CPU implementation currently supports

* atomic scattering
* excluded solvent correction
* hydration shell generation
* shell–shell scattering
* atom–shell scattering

The WebGPU implementation currently computes atomic Debye scattering and is under active development to support the complete SAXS model.

---


## Author

Toby E King, PhD
 
---

## Third-party libraries

This project vendors the Gemmi library for crystallographic and structural data handling.

Gemmi is developed by Global Phasing Ltd. and contributors:
https://github.com/project-gemmi/gemmi

See `checkout/gemmi/LICENSE.txt` for its license.


## License

This project is released under the Apache 2.0 license.
