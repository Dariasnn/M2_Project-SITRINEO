# SiTrInEO Development & Uncertainty Analysis

This repository is composed of **C programs** designed to be compiled and executed using **CERN's ROOT framework**. These programs are essential for the **development of the SITRINEO project** and for improving the **understanding of uncertainties** related to the experiment.

Additionally, the repository includes a **Graphical User Interface (GUI)** to assist in experiment visualization and data analysis.

## ⚙️ Compilation & Execution

### **🔧 Prerequisites**
- **CERN ROOT framework** ([Installation Guide](https://root.cern/install/))

### **🛠 Compiling the Programs**
To compile the C files with ROOT, navigate to the repository and run:

```bash
# Load ROOT environment (modify path if needed)
source /path/to/root/bin/thisroot.sh

# Start ROOT
root

# Compile the simulation of trajectory program or the Monte Carlo simulation
.L trajectory_final.C
.L ToyMC.C
.L fitting.C

# if error message on virtual machine like WSL, close root with .q and start again
````

### **Running the Programs**
```bash
# Run the needed function among the following
ToyMC()
trajectory()
trajectory3B()
trajectory_comparison()
trajectory_scattered()
simple_vs_scattered()
fitting()
trajectory_error()
trajectory_error3B()


````
## Purpose of the Programs

**Simulation program** (trajectory_final.C)

Models particle trajectories within the SITRINEO setup.
Implements multiple scattering effects and different B field conditions.
Outputs simulation data in ROOT format, and uncertainty results. 
```bash
trajectory() : trajectory through a uniforme B field
trajectory3B() : trajectory through a 3 region B field 
trajectory_comparison() : comparison of the two trajectories above
trajectory_scattered() : trajectory through a uniform B field but with multiple scattering at each plane
simple_vs_scattered() : comparison of the uncertainty on position with multiple scattering
fitting() : histograms and gaussian fits on simple and scattered positions
trajectory_error() : uncertainty on position with only a uniform B field
trajectory_error3B() : uncertainty on position using a 3 region B field

````
**Uncertainty on position** (fitting.C)
**Get global uncertainty** (ToyMC.C)
