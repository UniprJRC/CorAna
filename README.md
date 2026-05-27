# Calibrated Robust Correspondence Analysis with Applications

**Marco Riani**<sup>1</sup>, **Anthony C. Atkinson**<sup>2</sup>, **Francesca Torti**<sup>3</sup>, **Aldo Corbellini**<sup>1</sup>, **Gianluca Morelli**<sup>1</sup>

<sup>1</sup> Dipartimento di Scienze Economiche e Aziendale and Interdepartmental Centre for Robust Statistics, Universita di Parma, 43100 Parma, Italy

<sup>2</sup> Department of Statistics, The London School of Economics, London WC2A 2AE, UK

<sup>3</sup> European Commission, Joint Research Centre (JRC), 21027 Ispra, Italy

## Abstract

Correspondence analysis is a method for the visual display of information from two-way contingency tables. Often, the rows are subjects (in our major examples countries) and the columns are response categories. Our main result is a robustly calculated plot showing the structure of the data, including tests for outlying rows. We use simulation to calibrate the statistical properties of the procedure. The aim is both to detect outlying rows, if any, and to explore the homogeneity of the non-outlying rows.

## Keywords:

Forward search; inertia; informative plotting; outlier detection; simulation envelopes; trade restrictions.

---

This repository contains the MATLAB code associated with the paper. In the table below you can find the original MATLAB scripts (`.m`), the corresponding MATLAB live scripts (`.mlx`), and exported Jupyter notebook (`.ipynb`) files.

**MATLAB live script files**

The `.mlx` files contain both the code and the output produced by the code. To run a live script in MATLAB Online, click the "Open in MATLAB Online" button. The repository will be cloned automatically.

**Jupyter notebook files**

The `.ipynb` files contain both the code and the output produced by the code, including the generated figures and textual output from the live-script execution.

**Requirements**

To run the files below, install the [FSDA toolbox](https://www.mathworks.com/matlabcentral/fileexchange/72999-fsda-flexible-statistics-data-analysis-toolbox). The scripts use FSDA correspondence-analysis routines and the FSDA data sets `ExportShifts.mat`, `clothes.mat`, and `ChristmasSales.mat`. MATLAB with Optimization Toolbox is required, and Parallel Computing Toolbox is used by the figure-generation scripts.

| FileName | View :eyes: | Run | Jupyter notebook | m format |
| -------- | ----------- | --- | ---------------- | -------- |
| `createFigures1_4.mlx`: Code to reproduce Figures 1 to 4 by setting `FigNumber`. | [![File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://viewer.mathworks.com/addons/183971/1.0.0/files/createFigures1_4.mlx) | [![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=UniprJRC/CorAna&file=createFigures1_4.mlx) | [createFigures1_4.ipynb](https://github.com/UniprJRC/CorAna/blob/main/createFigures1_4.ipynb) | [createFigures1_4.m](https://github.com/UniprJRC/CorAna/blob/main/createFigures1_4.m) |
| `createFig5.mlx`: Code to reproduce Figure 5. | [![File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://viewer.mathworks.com/addons/183971/1.0.0/files/createFig5.mlx) | [![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=UniprJRC/CorAna&file=createFig5.mlx) | [createFig5.ipynb](https://github.com/UniprJRC/CorAna/blob/main/createFig5.ipynb) | [createFig5.m](https://github.com/UniprJRC/CorAna/blob/main/createFig5.m) |
| `Export_Shifts.mlx`: Code to reproduce Figures 6 to 10 using the Export Shifts data. | [![File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://viewer.mathworks.com/addons/183971/1.0.0/files/Export_Shifts.mlx) | [![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=UniprJRC/CorAna&file=Export_Shifts.mlx) | [Export_Shifts.ipynb](https://github.com/UniprJRC/CorAna/blob/main/Export_Shifts.ipynb) | [Export_Shifts.m](https://github.com/UniprJRC/CorAna/blob/main/Export_Shifts.m) |
| `Clothes.mlx`: Code to reproduce Figures 11 to 14 using the clothes data. | [![File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://viewer.mathworks.com/addons/183971/1.0.0/files/Clothes.mlx) | [![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=UniprJRC/CorAna&file=Clothes.mlx) | [Clothes.ipynb](https://github.com/UniprJRC/CorAna/blob/main/Clothes.ipynb) | [Clothes.m](https://github.com/UniprJRC/CorAna/blob/main/Clothes.m) |
| `ChristmasSales.mlx`: Code to reproduce Figure 15 using the Christmas Sales data. | [![File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://viewer.mathworks.com/addons/183971/1.0.0/files/ChristmasSales.mlx) | [![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=UniprJRC/CorAna&file=ChristmasSales.mlx) | [ChristmasSales.ipynb](https://github.com/UniprJRC/CorAna/blob/main/ChristmasSales.ipynb) | [ChristmasSales.m](https://github.com/UniprJRC/CorAna/blob/main/ChristmasSales.m) |

![GitHub top language](https://img.shields.io/github/languages/top/UniprJRC/CorAna)
![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/UniprJRC/CorAna)

[![GitHub contributors](https://img.shields.io/github/contributors/UniprJRC/CorAna)](https://github.com/UniprJRC/CorAna/graphs/contributors)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://github.com/UniprJRC/CorAna/graphs/commit-activity)
![main](https://img.shields.io/github/last-commit/UniprJRC/CorAna/main)
