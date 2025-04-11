
**"Fast Stable Parameter Estimation for Linear Dynamical Systems"**  
*Michelle Carey & James O. Ramsay*  
Published in *Computational Statistics & Data Analysis*, Volume 156, 2021.  
DOI: [10.1016/j.csda.2020.107124](https://doi.org/10.1016/j.csda.2020.107124)

---

# Data2LD-Matlab: Fast and Stable Estimation for Linear Dynamical Systems

This repository provides MATLAB code implementing the Data2LD (Data to Linear Dynamics) framework introduced in the aforementioned paper. The methodology offers a fast and stable approach for estimating parameters in linear ordinary differential equations (ODEs) from noisy and incomplete data, leveraging the linearity of the system to simplify computations and improve estimation accuracy.

##  Overview

The Data2LD framework builds upon the parameter cascading approach, where a linear combination of basis functions approximates the solution of the dynamical system. The system's parameters are then estimated to ensure that this approximating solution adheres to the observed data. Key features include:

- Efficient estimation of parameters in linear ODEs.
- Incorporation of B-spline basis functions for flexible function approximation.
- Implementation of an iterative scheme for fast and stable computation.
- Support for real-world applications, such as biomechanics data analysis.

## 📁 Repository Structure

The repository contains the following key files and directories:

- `Data2LD.m`: Main function implementing the Data2LD estimation procedure.
- `Data2LD_Opt.m`: Optimization routines for parameter estimation.
- `Data2LD_ISE.m`: Computes the Integrated Squared Error for model evaluation.
- `create_modelCell.m`: Constructs the model cell structure defining the ODE system.
- `inprod_basis_Data2LD.m`: Computes inner products of basis functions specific to Data2LD.
- `Setup:Install.m`: Script to set up the environment and install necessary components.
- `Examples/`: Directory containing example scripts and data demonstrating the usage of Data2LD.
- `fdaM.zip`: Archive containing the Functional Data Analysis MATLAB toolbox required for basis function operations.

## 🛠️ Requirements

- MATLAB (version R2014b or later recommended)
- Functional Data Analysis (FDA) MATLAB toolbox
- Basic knowledge of ODEs and numerical methods

## 🚀 Getting Started

1. **Clone the repository:**

   ```bash
   git clone https://github.com/mcareyucd/Data2LD-Matlab.git
   cd Data2LD-Matlab
   ```

2. **Install the FDA toolbox:**

   Unzip `fdaM.zip` and add the extracted directory to your MATLAB path.

3. **Run setup script:**

   In MATLAB, execute:

   ```matlab
   Setup:Install
   ```

4. **Run example scripts:**

   Navigate to the `Examples/` directory and run the desired example script to see Data2LD in action.

## 📈 Applications

The Data2LD framework is applicable in various domains where systems are modeled by linear ODEs, including:

- Biomechanics (e.g., modeling human movement)
- Pharmacokinetics/pharmacodynamics modeling
- Environmental modeling
- Engineering systems analysis

## 📖 Citation

If you utilize this code in your research, please cite the original paper:

Carey, M., & Ramsay, J. O. (2021). Fast stable parameter estimation for linear dynamical systems. *Computational Statistics & Data Analysis*, 156, 107124. [https://doi.org/10.1016/j.csda.2020.107124](https://doi.org/10.1016/j.csda.2020.107124)

## 📝 License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
