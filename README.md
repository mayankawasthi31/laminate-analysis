# ES - 671: Mechanics of Composite Materials - Course Project

This repository contains the code for the course project, developed in phases. The code performs micromechanics, macromechanics, and failure analysis of composite laminas and laminates.

---

## Phase 1 (Due: Apr 12, 2025)

### Objectives
Develop a code to perform macromechanical analysis of a single lamina.

### Tasks
* Take as inputs the **fibre and matrix properties** and the **angle of the lamina**.
* Compute the **engineering constants** for the lamina.
* Compute the **compliance ($[S]$) and stiffness ($[C]$) matrices**.
* Obtain the **transformed compliance ($[\bar{S}]$) and stiffness ($[\bar{Q}]$) matrices**.
* Given the global strains, compute the **local strains, global stresses, and local stresses** in the lamina.

---

## Extension of Phase 1 (Due: Apr 12, 2025)

### Objectives
Implement various failure criteria for the lamina.

### Tasks
* Considering the global and local stresses (computed in Phase 1):
    * Implement **Maximum Stress** failure criterion.
    * Implement **Maximum Strain** failure criterion.
    * Implement **Tsai-Hill** failure criterion.
    * Implement **Tsai-Wu** failure criterion.
* Identify the failure mode based on the criteria.
* If the lamina is safe, compute the **strength ratio**.

---

## Final Phase of the Project (Due: April 25, 2025)

### Objectives
Extend the code to analyze a multi-layered laminate using Classical Laminate Theory (CLT).

### Tasks
* Using the code developed so far:
    * Take as inputs the **lamina properties**, **number of layers**, **stacking sequence**, and **thickness** of each layer.
* Compute the **$[A]$, $[B]$, and $[D]$ matrices** (extensional, coupling, and bending stiffness matrices).
* Given the resultant forces $[N]$ and moments $[M]$, find the **midplane strains $[\epsilon^0]$ and curvatures $[\kappa]$**.
* Find the **global strains distribution** in the laminate based on midplane strains and curvatures.
* Find the **global stresses distribution** in the laminate based on the global strains.
