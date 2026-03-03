# Electromagnetic Potential Solver (Finite Element Method)

[cite_start]This project was developed as part of the **Differential and Difference Equations** course[cite: 1]. [cite_start]It implements the **Finite Element Method (FEM)** to solve a 1D boundary value problem (BVP) representing electromagnetic potential[cite: 3, 66].

## 1. Problem Definition
[cite_start]The project solves the following second-order differential equation[cite: 42, 71]:
$$\frac{d^2\phi}{dx^2} = -\frac{\rho}{\epsilon_r}$$

### Boundary Conditions:
* [cite_start]**Robin Condition** (at $x=0$): $\phi'(0) - \phi(0) = 2$ [cite: 71, 72]
* [cite_start]**Dirichlet Condition** (at $x=4$): $\phi(4) = 1$ [cite: 43, 71]

### Parameters:
* [cite_start]Charge density ($\rho$): $2$ [cite: 43]
* [cite_start]Relative permittivity ($\epsilon_r$): A piecewise constant function defined over the intervals $[0,1]$, $[1,2]$, and $(2,4]$[cite: 44, 73].

## 2. Methodology
[cite_start]The solution follows a standard FEM pipeline[cite: 5]:
1.  [cite_start]**Variational Formulation**: Deriving the weak form $b(w, v) = l(v)$ using integration by parts[cite: 7, 76, 88].
2.  [cite_start]**Discretization**: Using $n$ elements with linear basis functions (hat functions)[cite: 10, 103].
3.  **Numerical Integration**: 2-point **Gauss-Legendre quadrature** is used to calculate the stiffness matrix and load vector.
4.  **System Solver**: The resulting tri-diagonal system of linear equations is solved efficiently using the **Thomas Algorithm**.

## 3. Implementation
The solver is written in **C++**. 
* [cite_start]**Input**: The user provides the number of elements ($n$) as a runtime parameter[cite: 10].
* **Output**: The program generates a `.csv` file containing the calculated potential values $\Phi(x)$ for visualization.

## 4. Results
Visual representation of the electromagnetic potential for $n = 50$:

![Potential Plot](results/visualization.jpg)

## 5. How to Run
1.  Compile the source code:
    ```bash
    g++ src/main.cpp -o fem_solver
    ```
2.  Run the executable:
    ```bash
    ./fem_solver
    ```
3.  Enter the number of elements when prompted.
