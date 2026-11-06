# TMImix

TMImix is a MATLAB extension of the [SMI signal model](https://github.com/NYU-DiffusionMRI/SMI). Choose the compartments in a fascicle, convolve their response with a sampled fibre orientation distribution (fODF), and fit the normalized rotational invariants **S₀ and S₂** across b-values and diffusion times. TMImix runs alongside SMI and does not require the SMI toolbox.

Start MATLAB in this directory and run `addpath(pwd,'examples')`. The forward model and regression use base MATLAB; NLLS uses `lsqnonlin` from Optimization Toolbox. Tested with MATLAB R2026a.

## Run a noise experiment

Choose the active compartments and supply one **[b, Delta, Ndirs]** row per shell. This example uses all six compartments, eight diffusion-weighted shells and six b=0 repetitions.

```matlab
addpath(pwd,'examples')
parts = ["stick","ecs","permeableStick","sphere","dot","water"];
[b,Delta] = ndgrid([0.6 1.5 3 6],[18 40]);
Ndirs = repmat([12; 20; 32; 48],2,1);
protocol = [0 0 6; b(:) Delta(:) Ndirs];
c = noise_config(parts,protocol);
c.delta = 5;                 % PGSE pulse duration in ms
c.snr = [30 60];             % per-measurement SNR at b=0
c.mode = "free";
nlls = noise_nlls(c);
poly = noise_polynomial(c);
```

Every active parameter has an editable range in `c.prior`. With six compartments there are **14 independent fit parameters**, including five independent fractions because the six fractions sum to one. With one compartment its fraction is one. The default selection in `noise_config()` is permeable sticks, ECS and spheres, with all their parameters free. Set `c.prior.Dsoma=2`, for example, to keep that diffusivity fixed throughout a study. `c.mode="free"` uses the fit bounds directly and ignores `c.fixed`.

### Specify the test and training populations

`c.testPrior` controls held-out tissues for **both** methods. `c.trainPrior` controls the **polynomial training population**. Set only the fields you want to change. An unspecified field uses its fit bounds from `c.prior`. Test tissue may lie outside the fit bounds to study model mismatch. Training labels for estimated fields must lie within their fit bounds.

```matlab
c.testPrior.radius = [6 10];                 % uniform in µm
c.trainPrior.radius = @() 6 + 4*rand()^2;    % same support, more small radii
c.testPrior.tau = @() 10*exp(rand()*log(6)); % log-uniform from 10 to 60 ms
c.trainPrior.tau = [5 100];
c.testPrior.fractionAlpha = [2 4 2 2 1 1];  % joint Dirichlet, in parts order
```

Each field accepts a scalar for a point mass, a two-element interval, or a zero-argument function that returns one draw in physical units. A function can implement a discrete, bounded, or empirical distribution, for example `@() samples(randi(numel(samples)))` for a vector of measured radii.

Fractions share a unit total. Individual fraction callbacks draw first and preserve their values. For the remaining fractions you can specify a joint Dirichlet law with positive concentrations in `fractionAlpha`. It allocates the residual mass among free pools, conditional on fixed fractions and truncated to any specified fraction bounds. Without `fractionAlpha`, numeric intervals draw uniformly within the feasible remaining interval at each step. Explicit fraction fields are allocated first in declaration order, then other pools in compartment order. The last free fraction receives the remainder. To give one fraction a prescribed marginal, specify it individually and leave a pool to absorb the remaining mass. Incompatible draws raise an error.

### Compare radius constraints

Use the same `c.testPrior.radius=[6 10]` and seed for the following experiments. The mode applies to each field in `c.fixed`.

| Mode | MATLAB settings | Radius in test tissue | Radius held fixed in the fit |
|---|---|---|---|
| Matched | `c.mode="matched"; c.fixed=struct('radius',6);` | 6 µm | 6 µm |
| Oracle | `c.mode="oracle"; c.fixed="radius";` | Uniform from 6 to 10 µm | That tissue's true radius |
| Mismatched | `c.mode="mismatched"; c.fixed=struct('radius',7);` | Uniform from 6 to 10 µm | 7 µm |

Run `noise_nlls(c)` and `noise_polynomial(c)` after setting each mode. The two estimators see identical test measurements. Oracle and mismatched also share the same tissues and noise draws. If several parameters are known per tissue, set `c.fixed=["radius","Dperm"]` in oracle mode. A struct of field names and values also works in oracle mode, where only the names are used.

In matched and mismatched modes, polynomial training fixes the selected parameters at the estimator's prescribed values. In oracle mode, known parameters vary according to their training distributions and enter the regression as known covariates. Give those training distributions support over the values you intend to supply at prediction.

### Acquisition, noise and outputs

- Supply `c.protocol` as a numeric matrix or a table with columns named `b`, `Delta`, and `Ndirs`. Units are ms/µm² for b and ms for Δ. `c.delta` sets the rectangular PGSE pulse width as one value or one per row. Replace `c.protocol` to change the acquisition after creating the config. Each nonzero-b row uses its requested number of directions and requires at least six for S₂. The SH projection uses the highest supported even order up to `c.L` and retains S₀ and S₂ for fitting.
- Set `c.snr` to one or more b=0 SNRs. The unweighted signal is one, so each directional measurement gets independent noise with SD `1/SNR`. Set `c.noise="gaussian"` for real Gaussian noise or `"rician"` for complex Gaussian noise followed by magnitude. With b=0 rows, `Ndirs` gives the repeat count and the pooled noisy baseline mean normalizes each realization. With diffusion-weighted rows only, normalization uses the known baseline of one.
- Each simulated fODF has an antipodal fibre of mass `p2` and a uniform background of mass `1-p2`. The same fibre orientation applies to every shell of a tissue. NLLS starts at the true parameters, made feasible under the imposed constraints. Underdetermined protocols can be run to examine the local noise response around that start. Polynomial fitting uses an independent training set and a separate ridge regression at each SNR, with degree `c.degree` from 1 to 3.
- Set `c.nTissue`, `c.nRep`, and `c.nTrain` to control the test population, noise draws per tissue and training population. Results include the fitted `protocol`, held-out `truth`, `trainingTruth` for regression, `observed` invariants, `estimate`, and per-trial `b0Scale`. The `metrics` table gives bias, error SD, within-tissue SD and RMSE by parameter and SNR. NLLS also returns starts, exit flags and residual sums of squares. Saved configurations containing the original `model` and `nDirs` fields continue to use that model and its known baseline.

## Build a kernel and supply an fODF

The low-level forward constructor takes **[b, Δ, δ]** in ms/µm², ms and ms. Δ is the leading-edge separation of rectangular PGSE pulses; δ is their duration. Divide b-values in s/mm² by 1000. The model supports **LTE only**, a shared fODF across fascicles/pools, a common echo time and signals normalized to b=0. Different rows can have identical b and different timing.

```matlab
m = TMImix([0 20 5; 2 20 5; 2 50 5], ...
    ["stick","ecs","permeableStick","sphere","dot","water"]);
p = m.defaults();  % equal fractions; inspect p to see all active fields
p.fStick=.23; p.fECS=.35; p.fPerm=.20; p.fSphere=.15; p.fDot=.02; p.fWater=.05;
axes = [1 0 0; cosd(65) sind(65) 0]; mass = [.65; .35];
Y = m.signal(p,eye(3),axes,mass);  % three gradients x three shells
p.p2 = TMImix.odfP2(axes,mass);
S = m.invariants(p);                % three shells x [S0 S2]
K = m.kernel(p,linspace(0,1,101));  % full angular response per shell
```

| Compartment key | Fields | Response |
|---|---|---|
| `stick` | `fStick, Dstick` | Impermeable zero-radius stick |
| `ecs` | `fECS, DePar, DePerp` | Axially symmetric Gaussian; equal diffusivities give isotropic ECS |
| `permeableStick` | `fPerm, Dperm, tau` | Stick exchanging with the selected ECS; requires `ecs` |
| `sphere` | `fSphere, radius, Dsoma` | Reflecting sphere, Gaussian phase approximation with finite gradient pulses |
| `dot` | `fDot` | Immobile water; normalized signal is one at all b |
| `water` | `fWater, Dwater` | Isotropic free diffusion (default 3 µm²/ms) |

Diffusivities are in µm²/ms; radius is in µm; exchange time is in ms. Fractions are normalized b=0 signal fractions. For angular cosine μ, the stick diffusivity is D∥μ² and the ECS diffusivity is D⊥ + (D∥ − D⊥)μ². ECS contributes once: either as a Gaussian or as part of the exchanging pair. Other pools do not exchange. Choose `["stick","ecs","water"]` for the LTE SMI model, `["stick","ecs","sphere"]` with isotropic ECS for SANDI-like simulations, or `["permeableStick","ecs"]` for exchange alone.

**Exchange timing.** With η = fPerm/(fPerm + fECS), the stick-to-ECS rate is (1−η)/τ and the reverse rate is η/τ. Thus τ is the inverse **sum** of rates, not the stick residence time. `tau=Inf` disables exchange. The default `Exchange="finite"` propagates the two-site Kärger equations along the actual rectangular PGSE q(t)² history. The plateau is exact; each ramp uses 16 subintervals with exact mean q² and a matrix exponential. Increase `Steps` to check convergence for a demanding protocol. `Exchange="narrow"` uses a constant q² and effective time Δ−δ/3; δ=0 gives the strict narrow-pulse limit. The sphere uses the timing provided, independently of the exchange option, with 64 reflecting-sphere modes and an explicit δ=0 limit.

**Angular conventions.** The full forward signal is the probability-weighted sum of K(g·n) over fODF nodes n. Supply any list of unit-vector rows and nonnegative masses; for sampled fODF densities, multiply by the quadrature areas first. The code normalizes the masses and does not truncate the forward fODF or kernel at order 2. The `p2` field is ignored by `signal`: directions and weights determine it. For the invariant path, Kℓ = ∫₀¹ K(μ)Pℓ(μ)dμ; S₀ = K₀, S₂ = p₂|K₂|, in SMI's normalization. `Order` sets the number of Gauss–Legendre nodes (48 by default). Only p₂ of the fODF can enter this fit; S₀/S₂ do not determine the full fODF.

## Fit your own measurements

```matlab
% Y: baseline-normalized directional signals [directions x nonzero-b shells].
% dirs: matching unit direction rows. Choose an SH order supported by the sampling.
m = TMImix([.8 20 5; 2 20 5; 4 40 5],["stick","ecs"]); % match Y's shells
W = tmix.projector(dirs,4);
obs = tmix.measure(Y,W);   % one [S0,S2] row per shell in acquisition order
s = m.defaults();         % every scalar is fixed; open only desired parameters
s.fStick=[.05 .65]; s.fECS=[.05 .8]; s.Dstick=[.5 3]; s.p2=[0 1];
start = m.defaults();     % initial physical parameters
[estimate,info] = m.fitNLLS(obs,s,start);  % start: a full physical parameter struct
```

Supply one invariant row per protocol row. For shell-specific directions/counts, construct a projector and call `measure` for each shell separately. For b=0 use `[1 0]` after normalization; those rows carry no fit information. Nuisance even SH orders through `L` can reduce projection leakage; only the degree-0 and degree-2 norms enter the fit. Insufficient angular rank raises an error. NLLS defaults to equal weights; pass `Weights=w` (shells × 2, inverse variances) for weighted least squares.

A fit specification has exactly the fields returned by `m.defaults()`. A scalar fixes that physical value; an interval frees it within bounds. Select at least two free fractions to vary their distribution; with one free fraction its mass is determined by the others. Internally the fit allocates remaining mass successively in compartment order within feasible bounds, including bounds on subsequent pools. Uniform sampling of the independent coordinates induces an order-dependent fraction prior, not a uniform distribution over all feasible fractions. Keep the compartment order and record the prior when comparing estimators.

For custom supervised training, `R=tmix.PolyRegressor(m,s,S_train,p_train,Degree=2)` takes an array shaped shells × 2 × samples and matching truth structs. `R.predict(S_test)` returns one physical struct per observation. For per-observation constraints, pass `Known=["radius",...]` with scalar placeholders in `s`, vary those fields over the training population, and use `R.predict(S_test,knownStructs)`. Train on the same ordered protocol, noise model and preprocessing as the test data.

An arbitrary mixture can be simulated; that does not mean each of its parameters can be identified from S₀/S₂. For example, p₂ is unobservable in a purely isotropic signal, and τ is unobservable if either exchanging fraction vanishes. Start with the smallest plausible mixture, fix independently known quantities, and inspect noise propagation before opening additional parameters.

## Checks and references

Run `assertSuccess(runtests('tests'))` from the repository root. Tests cover finite-pulse exchange against an independent ODE, narrow-pulse exchange against a matrix exponential, SANDIX sphere fixtures, the SMI subset, full-fODF convolution vs SH extraction, noiseless recovery and all three constraint modes for both estimators.

The model follows [Coelho et al., SMI](https://doi.org/10.1016/j.neuroimage.2022.119290), [Jelescu et al., NEXI](https://doi.org/10.1016/j.neuroimage.2022.119277), [Olesen et al., finite-pulse exchange](https://doi.org/10.1016/j.neuroimage.2022.118976), and [Palombo et al., SANDI](https://doi.org/10.1016/j.neuroimage.2020.116835). The [SANDIX implementation](https://github.com/Neurophysics-CFIN/SANDIX) and [C2 protocol-design library](https://github.com/leehhtw/C2_protocoldesign) provided numerical/convention references.
