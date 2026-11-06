# TMImix

TMImix is a MATLAB extension of the [SMI signal model](https://github.com/NYU-DiffusionMRI/SMI). Choose the compartments in a fascicle, convolve their response with a sampled fibre orientation distribution (fODF), and fit the normalized rotational invariants **S₀ and S₂** across b-values and diffusion times. TMImix runs alongside SMI and does not require the SMI toolbox.

Start MATLAB in this directory and run `addpath(pwd,'examples')`. The forward model and regression use base MATLAB; NLLS uses `lsqnonlin` from Optimization Toolbox. Tested with MATLAB R2026a.

## Run a noise experiment

```matlab
addpath(pwd,'examples')
c = noise_config();
c.fixed = struct('radius',7);   % µm; use struct() for no overrides
c.testPrior.radius = [4 10];    % ground-truth range when allowed to vary
c.prior.radius = [4 10];        % cover that range in oracle regression training
c.mode = "oracle";              % "matched", "oracle", or "mismatched"
nlls = noise_nlls(c);
poly = noise_polynomial(c);
nlls.metrics(nlls.metrics.parameter=="tau",:)
```

| Mode | Radius in test tissue | Value supplied to the estimator |
|---|---|---|
| `matched` | Always 7 µm | Always 7 µm |
| `oracle` | Drawn from 4–10 µm for each tissue | That tissue's true radius, held fixed across noise repetitions |
| `mismatched` | Drawn from 4–10 µm for each tissue | Always 7 µm, even when it is wrong |

The flag applies to all fields of `c.fixed`. For example, `struct('radius',7,'Dperm',2,'fDot',0.04)` fixes three parameters in a model containing all three. Any other parameter can remain fixed across all modes by making it a scalar in both `c.prior` and `c.testPrior`. To fit radius freely, remove it from `c.fixed` and set a range in `c.prior`. A single scalar in `c.fixed` is a placeholder in oracle mode: set a range in both priors for a known parameter that should vary. Fraction bounds must still admit total mass one.

The default experiment uses an exchanging stick/ECS pair and impermeable spheres. It measures b = 0.8, 2.5, 6 ms/µm² at Δ = 16, 32, 56 ms (δ = 5 ms), with 64 LTE directions per shell. Eight test tissues have five noise repetitions at each SNR. Increase `nTissue` for population coverage and `nRep` for precision of within-tissue noise estimates. To test a different mixture, change `c.model` and rebuild both priors from that model's `defaults()`.

- **Noise generation.** Each synthetic fODF is a randomly oriented antipodal fibre of mass `p2` plus a uniform background of mass `1-p2`. Independent Gaussian noise is added to each directional measurement before projecting onto SHs. Set `c.noise="rician"` to take the magnitude after adding complex Gaussian noise. SNR is relative to a known, unit baseline; the experiment does not simulate noise from estimating b=0. Rician noise and the nonnegative S₂ norm can bias uncorrected NLLS.
- **NLLS.** Every repetition starts at its actual tissue parameters. Under a deliberately wrong constraint, the start keeps the imposed values and reconciles the other parameters with their bounds and the fraction sum. The actual start is returned. This measures local noise sensitivity; it does not demonstrate recovery from an unknown starting point.
- **Polynomial regression.** For each SNR a separate ridge regression learns from independent synthetic tissues, using S₀ and S₂ with all monomials and interactions through `c.degree` (1, 2 or 3). The default is quadratic. In matched/mismatched mode the training tissue obeys the assumed constraint. In oracle mode only the explicitly known parameters enter as extra covariates, so one conditional model spans the prior; unknown ground truth is never passed to prediction. Predictions obey the same fixed values, bounds and fraction sum as NLLS. Results depend on the protocol, prior, training size, polynomial degree and noise level.
- **Outputs.** `metrics` gives bias (estimate − truth), SD of errors across trials, pooled within-tissue SD, RMSE and nonconvergence counts. `truth`, `estimate`, `observed` and `trials` preserve every realization, including NLLS exit flags and squared residuals. `regressors` contains the trained predictors. Use `writetable(poly.metrics,'poly_metrics.csv')` to export. The two methods use identical test measurements for a given seed/configuration; oracle and mismatched modes also share identical test data.

As a reproducibility check, with the default seed and **SNR 60**, the two variable-radius modes gave the following results on the *same* 8 tissues × 5 noise draws. A wrong radius affected estimates of both sphere fraction and exchange time:

| Estimator | Sphere fraction RMSE: true radius supplied | Sphere fraction RMSE: assume 7 µm | τ RMSE: true radius supplied | τ RMSE: assume 7 µm |
|---|---:|---:|---:|---:|
| Ground-truth-start NLLS | 0.0055 | 0.0371 | 8.73 ms | 14.62 ms |
| Quadratic regression | 0.0066 | 0.0536 | 8.58 ms | 11.55 ms |

These are results for this synthetic protocol and prior, with other diffusivities fixed. To reproduce, run each estimator twice with `c.mode="oracle"` and `c.mode="mismatched"`; leave the other defaults unchanged.

## Build a kernel and supply an fODF

Each row of the protocol is **[b, Δ, δ]** in ms/µm², ms and ms. Δ is the leading-edge separation of rectangular PGSE pulses; δ is their duration. Divide b-values in s/mm² by 1000. The model supports **LTE only**, a shared fODF across fascicles/pools, a common echo time and signals normalized to b=0. Different rows can have identical b and different timing.

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
