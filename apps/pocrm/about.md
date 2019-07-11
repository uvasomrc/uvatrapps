
<br>

### Table of Contents

- [Getting Started](#getting-started)
- [Simulation](#simulation)
- [Implementation](#implementation)
- [How to cite POCRM](#how-to-cite)
- [References](#references)

---

<a name="getting-started"></a>
### Getting Started

Begin by specifying the number of dose levels of each agent, A and B, under investigation to form a matrix of K drug combinations labeled 1 through K. 

---

<a name="simulation"></a>

### Simulation

Input the hypothesized true DLT probabilities for each combination in the grid provided on the far left of the application. Specify the possible orderings for the DLT probabilities for your trial according to the labels given on the far left side of the application. It is imperative in using pocrm correctly that the user maintains the monotonicity assumption across rows and up columns of the matrix of combinations. This can be done by adhering to the across row, up columns, up/down diagonal specification technique outlined by Wages and Conaway (2013). On the far right side of the application, complete the design specifications based on the details of your trial. Click the **Simulate box** to generate operating characteristics under the scenario provided. Click on the question mark next to the possible orderings or simulation specifications for more information about each input.

---

<a name="implementation"></a>
### Implementation

Specify the possible orderings for the DLT probabilities for your trial according to the labels given on the far left side of the application. It is imperative in using pocrm correctly that the user maintains the monotonicity assumption across rows and up columns of the matrix of combinations. This can be done by adhering to the across row, up columns, up/down diagonal specification technique outlined by Wages and Conaway (2013). On the far right side of the application, complete the design specifications based on the details of your trial, and upload a .csv file with your current accumulated trial data. Click the **Get updated recommendation** box to generate an updated combination recommendation for the next participant based on the data provided. Click on the question mark next to the combination labels, possible orderings, or simulation specifications for more information about each input.

---

<a name="how-to-cite"></a>
### How to Cite

> Wages NA, Conaway MR, O’Quigley J. Dose-finding design for multi-drug combinations. Clin Trials 2011; 8: 380-9. PMCID: PMC3485079.

---
<a name="references"></a>
### Additional References

> Wages NA, Conaway MR, O’Quigley J. Continual reassessment method for partial ordering. Biometrics 2011; 67: 1555-63. PMCID: PMC3141101.

> Wages NA, Conaway MR. Specifications of a continual reassessment method design for phase I trials of combined drugs. Pharm Stat 2013; 12: 217-24. PMCID: PMC3771354.

> Wages NA, Varhegyi N. pocrm: an R package for phase I trials of combinations of agents. Comput Methods Programs Biomed 2013; 112: 211-8. PMCID: PMC3775989.


