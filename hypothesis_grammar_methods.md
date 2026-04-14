# Human Interpretable Grammar for Multicellular ABMs
## Method Implementation: English → Mathematics → Validation

> Johnson, Bergman, Rocha et al. *Cell* 188, 4711–4733 (2025)
> DOI: 10.1016/j.cell.2025.06.048

---

## Table of Contents

1. [Conceptual Architecture](#1-conceptual-architecture)
2. [The Vocabulary: Dictionaries of Signals and Behaviors](#2-the-vocabulary-dictionaries-of-signals-and-behaviors)
3. [The Grammar: From English Sentence to Mathematical Rule](#3-the-grammar-from-english-sentence-to-mathematical-rule)
4. [Single-Rule Mathematics](#4-single-rule-mathematics)
5. [Multi-Rule Mathematics](#5-multi-rule-mathematics)
6. [Response Function Options](#6-response-function-options)
7. [Implementation in PhysiCell](#7-implementation-in-physicell)
8. [Parameter Estimation and Sensitivity Analysis](#8-parameter-estimation-and-sensitivity-analysis)
9. [Data-Driven Model Initialization from Multi-Omics](#9-data-driven-model-initialization-from-multi-omics)
10. [Validation Case Studies](#10-validation-case-studies)
11. [Limitations and Future Directions](#11-limitations-and-future-directions)

---

## 1. Conceptual Architecture

### 1.1 The Core Problem

Agent-based models (ABMs) are powerful but inaccessible. Previously:
- Custom C++ code was required for each new model
- Biological hypotheses were buried deep in source code
- High-dimensional molecular data could not be systematically encoded into equations
- There was no reproducible mapping between human biological knowledge and mathematical rules

### 1.2 The Proposed Solution

A **cell behavior hypothesis grammar**: a controlled natural language that maps one-to-one onto mathematical expressions that can be automatically compiled into executable ABMs.

```
Human hypothesis statement
        ↓
  Grammar parser
        ↓
  Mathematical rule (Hill function)
        ↓
  PhysiCell ABM runtime
        ↓
  In silico simulation
```

### 1.3 Key Design Requirements

| Requirement | Implementation |
|---|---|
| Human-readable | Plain text CSV rules |
| Mathematically unambiguous | Controlled vocabulary + Hill response |
| Modular | Adding rules never modifies existing ones |
| Data-compatible | Cell type labels map directly to omics annotations |
| Reproducible | Rules saved to HTML/text after parsing |

---

## 2. The Vocabulary: Dictionaries of Signals and Behaviors

The grammar is built from two controlled dictionaries. Ambiguity is eliminated by restricting users to predefined terms.

### 2.1 Dictionary of Signals (S)

Signals are inputs—stimuli that a cell can sense. They are broadly categorized as:

| Category | Examples |
|---|---|
| Chemical factors | Oxygen, IL-10, IFN-γ, EGF, cisplatin |
| Mechanical cues | Pressure, volume |
| Cell-cell contact | Contact with dead cells, contact with specific cell types |
| Internal state | Cell volume, cycle phase |
| Accumulated damage | From effector attack |
| Simulation time | Used to trigger developmental timing events |

### 2.2 Dictionary of Behaviors (B)

Behaviors are outputs—quantifiable cell-scale processes parameterized by continuous rates, magnitudes, or probabilities:

| Category | Behaviors |
|---|---|
| Proliferation | Cycle entry rate, exit rate |
| Death | Apoptosis rate, necrosis rate |
| Motility | Migration speed, migration bias, chemotaxis direction |
| Mechanics | Adhesion, repulsion |
| Secretion/Uptake | Secretion rate, uptake rate for each substrate |
| State Transitions | Transition rate to another cell type/state |
| Cell-cell interactions | Phagocytosis rate, effector attack rate, fusion rate |
| Immunogenicity | Antigen presentation level |

### 2.3 Why This Vocabulary Works

The controlled vocabulary creates a **bijection**:

```
Biological noun (signal)   ←→  Mathematical variable s
Biological noun (behavior) ←→  Mathematical parameter b
Biological verb (verb)     ←→  Direction of response (+/-)
```

This means every valid English sentence in the grammar has exactly one mathematical representation, and every parameter in the resulting model has a biological interpretation.

---

## 3. The Grammar: From English Sentence to Mathematical Rule

### 3.1 Canonical Sentence Structure

```
In [cell type T], [signal S] increases/decreases [behavior B]
             [with half-max H] [and Hill power n]
```

**Examples:**

```
In MCF-7 breast cancer cells, cisplatin increases apoptosis.

In naïve T cells, IL-10 decreases transition to CD8+ T cells.

In tumor cells, oxygen increases cycle entry.

In tumor cells, pressure decreases cycle entry.

In M0-like macrophages, contact with dead cells increases 
    transition to M1-like macrophage.
```

### 3.2 CSV Encoding

Each sentence becomes one row in a CSV file with a fixed column schema:

```
cell_type ; signal ; response ; behavior ; max_response ; half_max ; hill_power ; applies_to_dead
```

**Concrete example:**

```csv
tumor ; oxygen ; increases ; cycle entry ; 0.0005 ; 5.0 ; 4 ; 0
tumor ; pressure ; decreases ; cycle entry ; 0.0 ; 0.25 ; 4 ; 0
MCF-7 ; cisplatin ; increases ; apoptosis ; 0.005 ; 0.5 ; 4 ; 0
naive_T_cell ; IL-10 ; decreases ; transition to CD8_T_cell ; 0.001 ; 0.25 ; 2 ; 0
```

**Column semantics:**

| Column | Mathematical Role | Notes |
|---|---|---|
| `cell_type` | Specifies which agents this rule applies to | Must match dictionary |
| `signal` | The variable s in the rule | Must match signal dictionary |
| `response` | Direction: increases → up-rule, decreases → down-rule | |
| `behavior` | The behavioral parameter b being modified | Must match behavior dictionary |
| `max_response` | The saturation value bM (for up) or bm (for down) | Units match behavior |
| `half_max` | s* — the signal level at half-maximal response | Units match signal |
| `hill_power` | n — controls the steepness/cooperativity | Typically 1–8 |
| `applies_to_dead` | Whether this rule applies after cell death | Boolean 0/1 |

---

## 4. Single-Rule Mathematics

### 4.1 The General Form

For a single rule relating behavior B to signal S:

$$b(s) = b_0 + (b_M - b_0) \cdot R(s)$$

Where:
- $b_0$ = base value of behavioral parameter in the absence of signal
- $b_M$ = maximally upregulated value (for an increasing rule)  
- $b_m$ = minimally downregulated value (for a decreasing rule)
- $R(s) \in [0, 1]$ = response function

### 4.2 The Default: Hill (Sigmoidal) Response Function

$$R(s) = H(s;\, s^*, n) = \frac{s^n}{(s^*)^n + s^n}$$

Where:
- $s^*$ = half-max value (signal level that produces 50% of maximal response)
- $n$ = Hill power (cooperativity/steepness parameter)

**Properties:**
- $R(0) = 0$ (no response in absence of signal)
- $R(s^*) = 0.5$ (half-max by definition)
- $R(\infty) \to 1$ (saturates at full response)
- Higher $n$ → sharper, more switch-like response

**Biological rationale:** Hill functions arise naturally from cooperative receptor binding and signaling network motifs. They are the standard dose-response curve in pharmacodynamics. A Hill power of 1 gives a simple Michaelis-Menten response; higher powers model ultrasensitivity.

### 4.3 Worked Example: Hypoxia-Induced Necrosis

English rule:
```
In tumor cells, oxygen decreases necrosis.
```

Mathematical instantiation:
$$\text{necrosis\_rate}(O_2) = b_M + (b_0 - b_M) \cdot H(O_2;\, O_2^*, n)$$

Where:
- $b_0$ = baseline necrosis rate (high, in absence of oxygen)
- $b_M$ = minimum necrosis rate (low, in excess oxygen)
- $O_2^*$ = oxygen level at half-max suppression of necrosis
- $n$ = Hill power

For the MDAMB-231 breast tumor calibration, $O_2^*$ and $b_0$ were derived from prior experimental validation (Rocha et al. 2021, iScience).

### 4.4 Worked Example: EGF-Induced Motility (Breast Cancer)

English rule:
```
In malignant epithelial cells, EGF increases migration speed.
```

Mathematical instantiation:
$$v_\text{mig}(\text{EGF}) = v_0 + (v_M - v_0) \cdot H(\text{EGF};\, \text{EGF}^*, n)$$

This is the "go hypothesis" — tested against the "grow hypothesis" (EGF increases proliferation instead) to determine which better reproduced macrophage-induced invasive tumor structures from DeNardo et al. 2009.

---

## 5. Multi-Rule Mathematics

### 5.1 The Extensibility Problem

When multiple signals regulate the same behavior, naïve approaches (e.g., summing Hill functions) can violate physical bounds or require modifying earlier rules when new rules are added.

The grammar solves this with a **bilinear interpolation** formulation that is:
- Bounded to $[b_m, b_M]$ by construction
- Modular: new rules can be appended without altering existing ones
- Reducible to single-rule form when only one signal is present

### 5.2 Pooling Up-Regulators

Given $m$ up-regulating signals $u_1, u_2, \ldots, u_m$ with respective half-maxes $u_1^*, \ldots, u_m^*$ and Hill powers $p_1, \ldots, p_m$, define the **total up-response**:

$$U = H_M(\mathbf{u};\, \mathbf{u}^*, \mathbf{p}) = \frac{\displaystyle\sum_{i=1}^{m} \left(\frac{u_i}{u_i^*}\right)^{p_i}}{1 + \displaystyle\sum_{i=1}^{m} \left(\frac{u_i}{u_i^*}\right)^{p_i}}$$

Note this is a **multivariate Hill function** — it reduces to the standard Hill function when $m=1$, and saturates at $U \to 1$ as any signal becomes dominant.

### 5.3 Pooling Down-Regulators

Given $n$ down-regulating signals $d_1, d_2, \ldots, d_n$:

$$D = H_M(\mathbf{d};\, \mathbf{d}^*, \mathbf{q}) = \frac{\displaystyle\sum_{j=1}^{n} \left(\frac{d_j}{d_j^*}\right)^{q_j}}{1 + \displaystyle\sum_{j=1}^{n} \left(\frac{d_j}{d_j^*}\right)^{q_j}}$$

### 5.4 Combined Response via Bilinear Interpolation

$$b(\mathbf{u}, \mathbf{d}) = (1 - D) \cdot \left[(1 - U) \cdot b_0 + U \cdot b_M\right] + D \cdot b_m$$

**Interpretation:**
1. The up-signals first set a **target** between $b_0$ and $b_M$: target $= (1-U) b_0 + U b_M$
2. The down-signals then **suppress** that target toward $b_m$
3. When $D=0$ (no inhibition), the expression reduces to the single up-regulator case
4. When $U=0$ (no activation), the expression reduces to the single down-regulator case

**Limiting cases:**

| U | D | b(u,d) |
|---|---|---|
| 0 | 0 | $b_0$ (baseline) |
| 1 | 0 | $b_M$ (fully activated) |
| 0 | 1 | $b_m$ (fully inhibited) |
| 1 | 1 | $b_m$ (inhibition dominates) |

### 5.5 Worked Example: Oxygen + Pressure Regulating Cycle Entry

Two rules applied simultaneously to the same behavior:

```
In tumor cells, oxygen increases cycle entry.
In tumor cells, pressure decreases cycle entry.
```

Results in:

$$b_\text{cycle}(O_2, P) = (1 - D_P) \cdot \left[(1-U_{O_2}) \cdot b_0 + U_{O_2} \cdot b_M\right] + D_P \cdot b_m$$

Where:
$$U_{O_2} = \frac{(O_2 / O_2^*)^{p_1}}{1 + (O_2 / O_2^*)^{p_1}}, \quad D_P = \frac{(P / P^*)^{q_1}}{1 + (P / P^*)^{q_1}}$$

This is the 3D surface shown in Figure S3 of the paper. The surface is bounded, smooth, and biologically interpretable at every point.

### 5.6 Non-Monotonic Response via Opposing Rules

For ECM-density-dependent PDAC cell motility (biphasic: motility first increases then decreases with collagen density), two opposing rules for the same signal are used:

```
In mesenchymal_tumor cells, ECM increases migration speed.  [low ECM]
In mesenchymal_tumor cells, ECM decreases migration speed.  [high ECM]
```

These are encoded as separate Hill functions with different half-max values. The U and D terms together create a non-monotonic (peaked) combined response that was fit to the experimental monoculture migration data from Panc10.05 cells across ECM densities of 2–7 mg/mL.

---

## 6. Response Function Options

The default Hill function can be replaced with linear or step responses:

### 6.1 Capped Linear Response

$$R_\text{linear}(s) = \begin{cases} 0 & s < s_\text{min} \\ \dfrac{s - s_\text{min}}{s_\text{max} - s_\text{min}} & s_\text{min} \leq s \leq s_\text{max} \\ 1 & s > s_\text{max} \end{cases}$$

Where $s_\text{min}$ and $s_\text{max}$ are threshold parameters replacing $s^*$ and $n$.

### 6.2 Step Function

$$R_\text{step}(s) = \begin{cases} 0 & s < s^* \\ 1 & s \geq s^* \end{cases}$$

Useful for discrete switch-like transitions (e.g., triggering a developmental event at a precise time).

### 6.3 Equivalence

Hill functions with high $n$ approximate step functions; linear functions can approximate Hill functions with $n \approx 1$ near the half-max. This interconvertibility (shown in Figure S2) means the grammar is flexible without sacrificing interpretability.

---

## 7. Implementation in PhysiCell

### 7.1 PhysiCell Agent Architecture

Each cell agent stores a **phenotype object**:

```
Phenotype
├── Cycle
│   ├── entry_rate
│   └── exit_rates
├── Death
│   ├── apoptosis_rate
│   └── necrosis_rate  
├── Volume
├── Mechanics
│   ├── adhesion
│   └── repulsion
├── Motility
│   ├── migration_speed
│   ├── migration_bias
│   └── chemotaxis_direction
├── Secretion
│   └── {substrate: secretion_rate, uptake_rate}
└── Cell_Interactions
    ├── phagocytosis_rate
    ├── attack_rate
    └── transition_rates
```

At each simulation timestep, each agent:
1. Samples signals from its local microenvironment (chemical concentrations, contact information, pressure)
2. Evaluates all applicable rules via the combined response function
3. Updates its phenotype parameters accordingly
4. Executes stochastic decisions based on updated rates

### 7.2 Rule Parsing Pipeline

```
CSV rule file
     ↓
  Parser reads each row
     ↓
  Validates: cell_type ∈ dict_cells, signal ∈ dict_signals, behavior ∈ dict_behaviors
     ↓
  For each agent of matching type, at each timestep:
      1. Sample signal value s from microenvironment
      2. Compute R(s) using stored {half_max, hill_power}
      3. Update behavioral parameter b using b₀, bM/bm, R(s)
     ↓
  Reproducibility: rules exported to HTML + text at parse time
```

### 7.3 Reaction-Diffusion Substrate Solver (BioFVM)

Chemical signals diffuse through the extracellular space. Each substrate obeys:

$$\frac{\partial c}{\partial t} = D \nabla^2 c - \lambda c + \sum_\text{cells} \left(S_i - U_i \cdot c\right) \delta(\mathbf{x} - \mathbf{x}_i)$$

Where:
- $D$ = diffusion coefficient
- $\lambda$ = decay rate
- $S_i$ = secretion rate of cell $i$
- $U_i$ = uptake rate of cell $i$

Cell agents contribute point sources/sinks at their positions. The secretion and uptake rates $S_i$, $U_i$ are themselves behavioral parameters modifiable by grammar rules.

### 7.4 Cell Mechanics

Cell position evolves under overdamped mechanics:

$$\frac{d\mathbf{x}_i}{dt} = \frac{1}{\xi} \left(\mathbf{F}_\text{adhesion} + \mathbf{F}_\text{repulsion} + \mathbf{F}_\text{motility}\right)$$

Where $\xi$ is a drag coefficient and motility force is biased random walk:

$$\mathbf{F}_\text{motility} = v_\text{mig} \left[\phi \hat{\mathbf{d}} + (1-\phi)\hat{\mathbf{r}}\right]$$

With $\phi$ = migration bias (0 = random walk, 1 = directed), $\hat{\mathbf{d}}$ = preferred direction (e.g., chemotaxis gradient), $\hat{\mathbf{r}}$ = random unit vector. Both $v_\text{mig}$ and $\phi$ are behavioral parameters controllable by grammar rules.

---

## 8. Parameter Estimation and Sensitivity Analysis

### 8.1 Sources of Parameter Values

The paper uses a hierarchy of evidence for parameter estimation:

| Priority | Source | Example |
|---|---|---|
| 1 (best) | Direct experimental measurement from the modeled system | PDAC motility assays vs. ECM density |
| 2 | Published experimental data from same cell type | MDA-MB-231 orthotopic tumor parameters (Rocha 2021) |
| 3 | Literature-derived values for related systems | T cell killing rates, macrophage polarization kinetics |
| 4 | Order-of-magnitude estimates with sensitivity analysis | Hypoxia model initial exploration |

### 8.2 Local Sensitivity Analysis

The paper quantifies parameter sensitivity via **multiplicative perturbation**:

$$p_j' = (1 + \epsilon) \cdot p_j, \quad \epsilon \in \{0.01, 0.05, 0.10, 0.20\}$$

Each perturbed parameter set is used to run replicate simulations. **Quantities of Interest (QoIs)** are computed from each run:

- Area under the curve (AUC) of live/dead cell populations over time
- Wasserstein distance between radial distributions of cell subtypes
- Spatial distribution metrics at endpoint

**Sensitivity Index (SI)** for parameter $j$:

$$\text{SI}_j = \frac{\overline{|\Delta \text{QoI}|}}{|\overline{\text{QoI}}|}$$

The mean and standard deviation of $\text{SI}_j$ across all QoIs is plotted to rank parameters by influence (Figure 2E, Figure S11E).

### 8.3 Key Sensitivity Finding: Half-Max Values Dominate

Across both the hypoxia model and the tumor-immune model, **half-max values** (the $s^*$ parameters in Hill functions) had the largest sensitivity indices, especially:

- Oxygen half-max for necrosis onset (`rule3_hfm`, `rule7_hfm`)
- Oxygen half-max for motile/non-motile phenotype transitions (`rule4_hfm`, `rule8_hfm`)
- Oxygen half-max for macrophage polarization (tumor-immune model)

This is biologically meaningful: the oxygen level at which cells cross phenotypic thresholds determines where in the tumor these behaviors manifest spatially.

**Base values** (the $b_0$ parameters) had consistently low sensitivity indices, indicating order-of-magnitude estimates are sufficient for qualitative behavior reproduction.

### 8.4 Motility Parameter Calibration from Tracking Data

For the PDAC model, motility parameters were fit directly from 3D single-cell tracking data:

1. Panc10.05 cells and HT-231 CAFs were embedded in collagen-I at 7 densities (2–7 mg/mL)
2. Cells tracked at 5-minute intervals over 8 hours
3. Trajectories fit to an **anisotropic persistent random walk model** (Wu et al. 2014, 2015)
4. The average speed parameter $v_\text{mig}$ extracted as a function of ECM density for each condition

**Fitting the biphasic response:**

The non-monotonic speed-vs-density curve was fit using MATLAB's `fmincon` (minimizing residual sum of squares) to find the parameters of two combined Hill functions:

$$v_\text{mig}(c) = v_0 + (v_M - v_0) \cdot H_\text{up}(c;\, c_1^*, n_1) - (v_\text{mig}(c) - v_m) \cdot H_\text{down}(c;\, c_2^*, n_2)$$

Where $c_1^* < c_2^*$: the first Hill function captures the motility increase at intermediate ECM, the second captures the suppression at high ECM density. The fitted curves are shown in Figure S9 for all four conditions (CAF monoculture, CAF coculture, PANC monoculture, PANC coculture).

---

## 9. Data-Driven Model Initialization from Multi-Omics

### 9.1 Cell Type Label Mapping

Because grammar rules use human-readable cell type names, these names can be **directly mapped to cell type annotations from omics data**:

```
scRNA-seq cluster annotation "epithelial_tumor"  
                    ↕  (same string)
Grammar cell type   "epithelial_tumor"
```

This is the key enabling insight for omics integration. No intermediate translation layer is needed.

### 9.2 Spatial Transcriptomics → ABM Initialization (PDAC)

For two human PDAC Visium samples (PDAC01, PDAC02 from GEO: GSE254829):

**Step 1: Spot annotation via transfer learning**  
ProjectR (transfer learning) applied NMF patterns from scRNA-seq to ST spots to distinguish:
- Proliferative signaling → `epithelial_tumor` phenotype
- EMT + inflammatory co-occurrence → `mesenchymal_tumor` phenotype

**Step 2: Fibroblast scoring**  
Seurat module scores computed from a pan-CAF gene signature per spot.

**Step 3: H&E image segmentation**  
CODA (3D reconstruction tool) used to annotate:
- Acinar cells, islet cells, smooth muscle (inert scaffolding)
- Collagen distribution → initial ECM density heatmap

**Step 4: Coordinate transformation**  
Affine linear transformation maps ST spatial coordinates (μm scale, tissue geometry) → ABM simulation domain, preserving aspect ratio.

**Step 5: ABM initialization**  
Each spot's cell type proportions seed agents at the corresponding domain location. Initial ECM density is set from collagen annotations.

### 9.3 scRNA-seq → Virtual Clinical Trial (PDAC)

For the immunotherapy simulation (Figure 6), immune cell compositions were derived from Steele et al. 2020 (GEO: GSE155698):

**Cell subtyping:**
- T cells classified hi/lo by median-normalized expression of `PDCD1` (PD-1) and `TNFRSF9` (CD137)
- Tumor cells classified by `CD274` (PD-L1) expression
- Resulting 4 CD8+ T cell subtypes (PD-1hi/lo × CD137hi/lo) × 2 tumor subtypes × macrophages

**Virtual patient initialization:**
- 1,000 tumor cells per simulation
- PD-L1hi:PD-L1lo ratio set from per-patient scRNA-seq measurements
- Immune cell counts set directly from per-patient cell type proportions

**Therapy simulation:**
- GVAX → double all T cell populations
- Nivolumab (ICI) → convert PD-1hi agents to PD-1lo; convert PD-L1hi tumor to PD-L1lo
- Urelumab (URU) → convert CD137lo agents to CD137hi

### 9.4 Allen Brain Atlas → Cortical Development Calibration

For the neurodevelopment model, the Allen Brain Atlas was used differently — **not for initialization but for parameter fitting**:

1. Single z-slices extracted from somatosensory (SOM) and auditory (AUD) cortex
2. Rectilinear subregions isolated; layer annotations used to count cells per layer
3. Target: layer cell counts at simulation endpoint (fully formed cortex)
4. **Fitting objective:** minimize residual sum of squares between simulated final layer thicknesses and atlas-derived layer thicknesses

$$\hat{\theta} = \arg\min_\theta \sum_\text{layers} \left(T_\text{layer}^\text{simulated}(\theta) - T_\text{layer}^\text{atlas}\right)^2$$

---

## 10. Validation Case Studies

### 10.1 Hypoxia Model: Qualitative Replication + Sensitivity Robustness

**English hypothesis set:**
```
In tumor cells, oxygen decreases necrosis.
In tumor cells, oxygen decreases transformation to motile tumor cells.
In motile tumor cells, oxygen increases transformation to tumor cells.
In tumor cells, oxygen increases cycle entry.
In motile tumor cells, oxygen decreases migration speed.  [paradoxically: 
                                                           post-hypoxic cells seek
                                                           low-oxygen regions]
```

**Mathematical model:** 8-rule system generating 24-dimensional parameter space.

**Validation approach:**  
Parameters from prior experimental calibration (Rocha et al. 2021) using MDA-MB-231 orthotopic murine breast tumor data. Simulation run for 5 days from 2,000 cells in 38 mmHg oxygenation.

**Predicted behaviors (qualitatively validated):**
- Oxygen-poor necrotic core forms at tumor center ✓
- Hypoxic cells concentrate near peri-necrotic boundary ✓
- Post-hypoxic cells form invasive "plumes" in normoxic regions ✓
- Plumes fail to exit tumor when hypoxic memory is short ✓ (consistent with Godet et al. 2019)

**Quantitative robustness:**  
QoI medians (AUC of populations, Wasserstein distances) remained stable under 1–20% parameter perturbations. Variability increased with perturbation magnitude but central tendency was preserved, confirming that qualitative behaviors are robust to order-of-magnitude parameter uncertainty.

### 10.2 PDAC CAF Model: Experimental Validation of Invasion Hypothesis

**English hypothesis set:**
```
In fibroblasts, ECM increases ECM secretion.
In epithelial_tumor cells, fibroblast_factor increases 
    transformation to mesenchymal_tumor cells.
In mesenchymal_tumor cells, ECM increases migration speed.
In mesenchymal_tumor cells, ECM decreases migration speed.  [biphasic]
In epithelial_tumor cells, inflammatory_factor increases cycle entry.
```

**Experimental validation — Step 1: Motility assay**  
PANC10.05 cells co-cultured with HT-231 CAFs in collagen-I at 7 ECM densities. Co-culture increased motility vs. monoculture at all but the highest ECM density (p<0.05 for 14/15 density comparisons, 30-frame trajectories). This supported the hypothesis that CAF-secreted factors increase tumor cell motility.

**Experimental validation — Step 2: Organoid invasion assay**  
15 patient-derived PDAC organoids (PDOs) treated with iCAF-conditioned and myCAF-conditioned media vs. control (collagen alone). Result: significantly increased invasion in CAF-conditioned media (Wilcoxon p = 6.1×10⁻⁵). Both iCAF and myCAF induced similar invasion levels (p = 1 between subtypes), validating the single-CAF-subtype abstraction.

**In silico prediction validated:**  
Simulations predicted that higher fibroblast density initially promotes invasion but at very high ratios (1:10 PANC:CAF) the dense ECM physically blocks invasion — matching the experimental observation that co-culture at highest ECM density showed decreased motility vs. lower densities.

### 10.3 EGF/EGFR Breast Cancer: Hypothesis Discrimination

**Two competing hypotheses tested in silico:**

| Hypothesis | English Rule | Predicted Behavior |
|---|---|---|
| "Grow" | EGF increases cycle entry in epithelial cells | Expanded tumor volume, no invasion |
| "Go" | EGF increases migration speed in epithelial cells | Invasive protrusions matching DeNardo 2009 |
| "Go+Grow" | Both simultaneously | Invasive protrusions + expansion |

**Computational conclusion:** Only "Go" and "Go+Grow" reproduced the macrophage-induced invasive structures.

**Experimental validation:**
1. MMTV-PyMT organoids treated with gefitinib (EGFR inhibitor): dose-dependent reduction in invasion (Figure 5F), but no effect on colony formation (IC50 = 2.48 μM), confirming motility as primary EGFR effect
2. MCF10A cells + EGF: both increased motility AND proliferation confirmed (Figures 5G, S14), supporting the "Go+Grow" hypothesis

**Quantitative validation:**  
Gefitinib colony formation assay showed no colony reduction until extremely high doses (IC50 = 2.48 μM), while invasion was inhibited at low doses — confirming the computational prediction that EGFR signaling acts primarily through motility (go) rather than proliferation (grow) in this invasion context.

### 10.4 Virtual Clinical Trial: Macrophage Hypothesis

**Novel hypothesis generated by simulation:**  
From the PDAC virtual clinical trial (Figure 6G), simulations predicted that **macrophage abundance** (not CD8+ T cell abundance) correlates with response to triple therapy (GVAX + ICI + URU). Specifically:
- Macrophage percentage significantly higher in simulated responders (p = 0.0056)
- CD8+ T cell percentage showed no significant difference (p = 0.5509)

**Biological hypothesis generated:** Macrophage clearing of dead tumor cells is essential for lymphocyte trafficking and effective cytotoxic killing.

**Consistency with clinical data:**  
This computationally generated hypothesis is consistent with published clinical observations of increased TREM2+ macrophage signaling in the triple combination arm of the Heumann et al. 2023 neoadjuvant PDAC trial.

### 10.5 Neurodevelopment: Quantitative Layer Count Matching

**English hypothesis set (key rules):**
```
In stem cells, time increases asymmetric division rate.
In stem cells, time decreases symmetric division rate.
In progenitor cells, time increases transition to layer_N cells.
In layer_N cells, pial_contact decreases migration speed.
```

**Calibration approach:** Residual sum of squares minimization against Allen Brain Atlas layer counts for SOM and AUD cortex regions.

**Result:** Simulated layer cell counts at endpoint matched atlas-derived counts for both SOM and AUD regions (Figure 7D), demonstrating that the same grammar framework that models cancer can capture the spatial and temporal dynamics of cortical laminar formation when parameterized from static atlas data.

---

## 11. Limitations and Future Directions

### 11.1 Current Grammar Limitations

| Limitation | Biological Impact | Proposed Extension |
|---|---|---|
| Statements treated as independent (inclusive OR) | Cannot encode AND-gated or REQUIRES logic | Additional operators: AND, REQUIRES |
| No hysteresis | Cannot model persistent epigenetic states | Hysteresis extensions to response curves |
| No access to contacting cell's properties as signals | Cannot model delta-Notch signaling | Allow `signal = property_of_contacting_cell` |
| No gene regulatory networks | Cannot model transcription factor cascades | Integration with PhysiBoSS Boolean networks |
| No pharmacokinetics | Virtual clinical trials lack drug distribution | Extension to quantitative systems pharmacology models |
| All rules currently active simultaneously | No temporal sequencing of rule activation | `time` signal can partially address this |

### 11.2 Data Integration Gaps

- Current data connections use macroscale summaries (cell type counts, layer thicknesses, tumor size)
- Future need: spatially, temporally, and phenotypically deeper connections to data
- Parameter fitting methods needed: Bayesian inference, data assimilation for omics-informed ABMs
- Sensitivity analysis of cellular phenotype granularity (how many subtypes are necessary?)

### 11.3 LLM Integration Opportunity

The paper explicitly notes that LLMs (e.g., ChatGPT) could facilitate **translation** from familiar biological language (e.g., "fibrosis") into the constrained grammar vocabulary — a natural agentic application given the structured, bounded output space required.

---

## Appendix: Complete Mathematical Reference

### A.1 Single Up-Regulating Rule

Given rule: *"In T, S increases B"* with parameters $(b_0, b_M, s^*, n)$:

$$b(s) = b_0 + (b_M - b_0) \cdot \frac{s^n}{(s^*)^n + s^n}$$

### A.2 Single Down-Regulating Rule

Given rule: *"In T, S decreases B"* with parameters $(b_0, b_m, s^*, n)$:

$$b(s) = b_0 - (b_0 - b_m) \cdot \frac{s^n}{(s^*)^n + s^n}$$

Equivalently: $b(s) = b_m + (b_0 - b_m) \cdot \left[1 - \frac{s^n}{(s^*)^n + s^n}\right]$

### A.3 General Multi-Rule Combined Response

$$\boxed{b(\mathbf{u}, \mathbf{d}) = (1 - D) \cdot \left[(1 - U) \cdot b_0 + U \cdot b_M\right] + D \cdot b_m}$$

With:

$$U = \frac{\displaystyle\sum_{i=1}^{m} \left(\frac{u_i}{u_i^*}\right)^{p_i}}{1 + \displaystyle\sum_{i=1}^{m} \left(\frac{u_i}{u_i^*}\right)^{p_i}}, \qquad D = \frac{\displaystyle\sum_{j=1}^{n} \left(\frac{d_j}{d_j^*}\right)^{q_j}}{1 + \displaystyle\sum_{j=1}^{n} \left(\frac{d_j}{d_j^*}\right)^{q_j}}$$

### A.4 CSV Row Schema

```
cell_type ; signal ; increases|decreases ; behavior ; b_sat ; s_half ; hill_n ; dead_flag
```

### A.5 Parameter Roles Summary

| Symbol | Name | Grammar Column | Role |
|---|---|---|---|
| $b_0$ | Base value | (inferred from PhysiCell defaults) | Behavior in absence of all signals |
| $b_M$ | Saturation value | `max_response` | Maximum behavior under full up-regulation |
| $b_m$ | Minimum value | `max_response` for down-rules | Minimum behavior under full down-regulation |
| $s^*$ | Half-max | `half_max` | Signal level producing 50% of max response |
| $n$ | Hill power | `hill_power` | Cooperativity / steepness |
| $U$ | Up-response | (computed) | Aggregate activation level ∈ [0,1] |
| $D$ | Down-response | (computed) | Aggregate inhibition level ∈ [0,1] |

---

*Document compiled from Johnson et al. (2025), Cell 188:4711–4733. Code available at https://github.com/physicell-models/grammar_samples*
