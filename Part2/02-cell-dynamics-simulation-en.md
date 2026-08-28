# Part 2 — Two-Dimensional Cell-Dynamics Simulation (Simulation Core)

## 1. Positioning and Role in the Closed Loop

This part is the simulation core of the closed loop. It translates cell-scale behavioural rules such as proliferation, migration, and contact inhibition into an observable macroscopic tumour morphology. Its input is a set of interpretable cell-dynamics parameters, and its output is the two-dimensional morphology of a virtual tumour and its boundary, which Part 3 measures with the same fractal algorithm that is applied to real images. The scientific task of this part is not to reproduce the lesion of a specific patient but to build a mechanistic model with tunable parameters and variable morphology, so that hypotheses such as loss of contact inhibition and enhanced migration can be tested through quantifiable changes in shape.

## 2. Method and Technical Route

The model adopts a hybrid cellular-automaton framework, representing individual cancer cells discretely on an N×N two-dimensional lattice and describing the nutrient or oxygen concentration with a continuous reaction-diffusion field. This framework inherits the hybrid mathematical model proposed by Anderson, whose core idea is that the discrete states of cells on the lattice and the continuous microenvironment field are coupled through consumption and diffusion. Each cell carries attributes such as proliferation rate, migration probability, contact-inhibition strength, and an adhesion threshold, which jointly determine at each time step whether the cell divides, whether it migrates, and how it responds to the local nutrient.

The proliferation rule is gated by maturity and available space: a cell becomes capable of division only after reaching the cell-cycle duration threshold, and division is allowed only when an empty neighbouring site exists, otherwise division is delayed and the cell becomes quiescent. The contact-inhibition factor determines whether mechanical crowding blocks division. Benign nevi retain this inhibition, so cells stop proliferating when neighbours are occupied, whereas malignant phenotypes weaken or lose the inhibition and divide even when surrounded by cells of the same type, pushing outward; this is the mechanistic source of the malignant finger-like infiltrative morphology. The migration rule is a random walk on the lattice, modified by a chemotactic or haptotactic bias term that steers the direction of movement along nutrient or matrix gradients; higher migration probability and lower adhesion threshold make cells more likely to detach from the aggregate and spread outward. The nutrient field is updated at every time step by solving a diffusion-minus-consumption reaction-diffusion equation, and cells become quiescent or apoptotic when the local nutrient falls below a threshold.

| Parameter | Meaning | Benign vs malignant tendency |
|---|---|---|
| Proliferation rate α | Probability or rate of division | Higher in malignant cells |
| Migration probability P_move | Tendency to move per time step | Higher in malignant cells |
| Contact-inhibition factor β | Strength of mechanical crowding blocking division | Retained in benign, lost in malignant |
| Adhesion threshold A_i | Neighbour constraint before migration | High in benign, low in malignant |
| Cell-cycle duration T_cell | Minimum interval between divisions | Shorter in malignant cells |
| Apoptosis rate δ | Probability of cell death | Can be lower in malignant cells |
| Nutrient diffusion/consumption | Physical parameters of the field C(x,y) | Determines necrotic core, not overall size |

## 3. Key Algorithm and Pseudocode

Pseudocode for a single time-step dynamics loop follows.

```
Input:  cell lattice L, nutrient field C, parameters θ = {α, P_move, β, T_cell, δ}
Output: updated L and C

for each cell in random order:
    consume local nutrient; if C < C_min mark quiescent or apoptotic (prob. δ)
    if cell age >= T_cell and an empty neighbour exists:
        apply contact inhibition: if β high and neighbours occupied, delay division
        otherwise divide with probability α, one daughter occupies a random empty neighbour
for each cell in random order:
    migrate with probability P_move to a random empty neighbour
    (or choose direction weighted by chemotactic/haptotactic bias)
numerically solve the nutrient field: dC/dt = diffusion - consumption
```

## 4. Expected Output and Evaluation Metrics

The output is the two-dimensional morphology of the virtual tumour and its boundary. Under a parameter combination of high proliferation, weak migration, and strong contact inhibition, the tumour is expected to grow as a regular circle; under faster proliferation, weakened contact inhibition, and enhanced migration, the boundary is expected to transition from smooth to finger-like infiltration. Quantitative evaluation uses the same metrics as Part 1: the fractal dimension of the simulated boundary, the compactness index, and the radius growth rate, which serve as the outputs for the parameter sweep and alignment in Part 3.

| Output | Definition | Purpose |
|---|---|---|
| Simulated boundary fractal dimension D_f | Box-counting on the simulated contour | Alignment with real-image descriptors |
| Tumour radius and cell count | Macroscopic-scale evolution | Calibrating growth dynamics |
| Boundary morphology | Circular vs finger-like | Qualitative reading of mechanism effects |
| Nutrient-field snapshot | Spatial distribution of C(x,y) | Explaining necrotic-core formation |

## 5. Feasibility, Risks, and Alternatives

The cellular automaton is simple to implement and computationally efficient, and even the hybrid framework with a continuous field can complete a single simulation within minutes on an ordinary workstation, so feasibility is high. The main risk is that lattice discreteness introduces artefacts tied to the grid orientation, giving the circular morphology a spurious anisotropy. The alternative is an off-lattice single-cell biophysical model, which treats each cell as an individual object in continuous space with mechanical interactions; this approach represents contact inhibition more faithfully but raises the implementation and calibration cost substantially, and can serve as a methodological upgrade in a later phase. Another risk is parameter identifiability, since several parameters can produce equivalent effects on the same morphology; this must be resolved in Part 3 through sensitivity analysis.

## 6. Data and Ethical Compliance

This part is purely computational and involves no human subjects or clinical data. Parameter values must be grounded in published cell-biology and melanoma-modelling literature rather than set arbitrarily; any value without a literature basis is marked with a to-be-determined placeholder until it is measured experimentally or confirmed against the literature.

## 7. References

[1] A. R. A. Anderson, "A hybrid mathematical model of solid tumour invasion: the importance of cell adhesion," Math. Med. Biol., vol. 22, no. 2, pp. 163–186, 2005, doi: 10.1093/imammb/dqi005.

[2] D. Drasdo and S. Höhme, "A single-cell-based model of tumor growth in vitro: monolayers and spheroids," Phys. Biol., vol. 2, no. 3, pp. 133–147, 2005, doi: 10.1088/1478-3975/2/3/001.

[3] J. Wang, L. Zhang, C. Jing, G. Ye, and H. Wu, "Multi-scale agent-based modeling on melanoma and its related angiogenesis analysis," Theor. Biol. Med. Model., vol. 10, 41, 2013, doi: 10.1186/1742-4682-10-41.
