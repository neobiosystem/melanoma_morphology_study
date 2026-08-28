# Part 3 — Quantitative Alignment of Simulation to Real Images (Core Result)

## 1. Positioning and Role in the Closed Loop

This part is the core-result stage of the closed loop and carries the cross-scale mapping between image features and biological mechanism. Its input is the fractal-dimension ranges of real benign and malignant lesions obtained in Part 1, together with the parameter-driven virtual tumour morphologies produced in Part 2; its output is a phase diagram showing how parameters such as contact inhibition, migration, and proliferation determine the fractal dimension of the simulated boundary, and whether that trend matches the observed change from benign to malignant in real images. This part measures simulated boundaries with exactly the same fractal algorithm as Part 1, thereby removing the systematic bias that could arise from algorithmic differences and giving the alignment an interpretable causal meaning.

## 2. Method and Technical Route

The key to alignment is consistency of measurement. For every virtual tumour produced in Part 2, the fractal dimension is computed on the simulated contour with the same box-counting method used in Part 1, yielding the simulated fractal dimension. The parameter space is then swept using grid sampling or Latin hypercube sampling to cover the plausible ranges of proliferation rate, migration probability, and contact-inhibition factor, running many simulations and recording the simulated fractal dimension for each parameter combination.

A phase diagram is constructed from the sweep, with the contact-inhibition factor and migration probability as axes and the simulated fractal dimension as the response, to see whether there exists a boundary in the phase diagram separating a low-fractal-dimension regular-circular region from a high-fractal-dimension finger-like infiltration region. The measured benign and malignant fractal-dimension ranges from Part 1 are overlaid on the phase diagram to test whether the simulation can reproduce the benign-to-malignant transition along a parameter path. To disentangle the relative contributions of the parameters, a variance-based global sensitivity analysis is applied, computing first-order and total-effect indices to identify the parameters and interactions that most influence the fractal dimension. Alignment quality is measured as the error between simulated fractal dimensions and the observed ranges.

| Parameter / step | Meaning | Note |
|---|---|---|
| Measurement consistency | Same fractal algorithm for simulation and real images | Box-counting, removing algorithmic bias |
| Parameter sweep | Covering the parameter space | Grid or Latin hypercube sampling |
| Phase-diagram axes | Visualising the response surface | β, P_move as axes, D_f as response |
| Sensitivity analysis | Relative parameter contributions | Sobol'/Saltelli first-order and total-effect indices |
| Alignment error | Distance between simulation and observation | RMSE or range coverage |

## 3. Key Algorithm and Pseudocode

Pseudocode for the parameter sweep, phase-diagram construction, and sensitivity analysis follows.

```
Input:  parameter space Ω = {θ_1 .. θ_M} (grid or Latin hypercube), observed D_f ranges
for m = 1 .. M:
    run the Part-2 simulation with parameters θ_m, obtaining a virtual tumour
    compute D_f_sim(m) with the Part-1 box-counting method
build phase diagram: (β, P_move) as coordinates, D_f_sim as the response
overlay: superpose the observed benign/malignant D_f ranges from Part 1
sensitivity: compute Sobol' first-order S_i and total-effect S_Ti indices
alignment error: RMSE = sqrt(mean((D_f_sim - D_f_obs)^2))
output: phase diagram + parameter importance ranking + alignment error
```

## 4. Expected Output and Evaluation Metrics

The expected output is a phase diagram and a set of sensitivity indices. The phase diagram should show that the simulated fractal dimension rises monotonically as contact inhibition weakens or migration strengthens, and that its range covers and matches the observed benign-to-malignant range from Part 1, thereby explaining mechanistically why the boundary fractal dimension increases. The sensitivity analysis is expected to show that contact inhibition and migration have the most prominent total effects on the fractal dimension, with the effect of proliferation appearing partly indirectly through spatial crowding.

| Metric | Definition | Interpretation |
|---|---|---|
| Simulated fractal dimension D_f_sim | Box-counting on simulated contour | Degree of overlap with observed ranges |
| Alignment error RMSE | Distance between simulated and observed D_f | Smaller is better |
| Phase-diagram transition | Boundary between circular and fingering regions | Whether it matches the benign/malignant parameter threshold |
| Sobol' total effect S_Ti | Total variance contribution of a parameter | Ranking to identify dominant parameters |

## 5. Feasibility, Risks, and Alternatives

The feasibility of the method hinges on the balance between the cost of a single simulation and the size of the sweep. If the parameter space is high-dimensional, a full grid sweep becomes intractable; in that case Latin hypercube sampling can control the scale, or a coarse grid can first locate the transition region and then be refined locally. Global sensitivity analysis requires a large number of Monte Carlo samples and is computationally demanding; the alternative is to run a one-at-a-time screening first and then compute Sobol' indices only for the few dominant parameters. Another risk is that the observed fractal-dimension ranges themselves carry uncertainty from segmentation quality, so the phase-diagram overlay should present confidence ranges rather than single values, to avoid mistaking noise for mechanistic signal.

## 6. Data and Ethical Compliance

This part reuses the public image descriptors from Part 1 and the simulation results from Part 2, adds no new human-subject data, and carries no ethical risk. The phase diagram and sensitivity results must be reported faithfully, including goodness of fit and range uncertainty, and must not present only the subset of parameters consistent with the hypothesis.

## 7. References

[1] A. Brú, J. M. Pastor, I. Fernaud, I. Brú, and S. Melle, "Super-rough dynamics on tumor growth," Phys. Rev. Lett., vol. 81, no. 18, pp. 4008–4011, 1998, doi: 10.1103/PhysRevLett.81.4008.

[2] T. Homma and A. Saltelli, "Importance measures in global sensitivity analysis of nonlinear models," Reliab. Eng. Syst. Saf., vol. 52, no. 1, pp. 1–17, 1996, doi: 10.1016/0951-8320(96)00002-6.
