# SiTrInEO Development & Uncertainty Analysis

This repository is composed of **C programs** designed to be compiled and executed using **CERN's ROOT framework**. These programs are essential for the **development of the SITRINEO project** and for improving the **understanding of uncertainties** related to the experiment.

Additionally, the repository includes a **Graphical User Interface (GUI)** to assist in experiment visualization and data analysis.

## ⚙️ Compilation & Execution

### **🔧 Prerequisites**
- **CERN ROOT framework** ([Installation Guide](https://root.cern/install/))
- **GCC Compiler** (`sudo apt install gcc` for Linux)

### **🛠 Compiling the Programs**
To compile the C files with ROOT, navigate to the repository and run:

```bash
# Load ROOT environment (modify path if needed)
source /path/to/root/bin/thisroot.sh

# Start ROOT
root

# Compile the simulation of trajectory program for example
.L trajectory_final.C
````

### **Running the Programs**
```bash
# Run the needed function
trajectory()

````
## Purpose of the Programs

**Simulation program** (trajectory_final.C)

Models particle trajectories within the SITRINEO setup.
Implements multiple scattering effects and different B field conditions.
Outputs simulation data in ROOT format, and uncertainty results. 

