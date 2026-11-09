# TMImix

TMImix is a MATLAB extension of the [SMI biophysical model](https://github.com/NYU-DiffusionMRI/SMI). Choose the compartments in a fascicle, convolve their response with a sampled fibre orientation distribution (fODF), and fit the normalized rotational invariants **S₀ and S₂** across b-values and diffusion times. TMImix runs alongside SMI and does not require the SMI toolbox.

Start MATLAB in this directory and run `addpath(pwd,'examples')`. The forward model and regression use base MATLAB; NLLS uses `lsqnonlin` from Optimization Toolbox. Tested with MATLAB R2026a.

## Fit exchanging tissue with impermeable sticks

Run `run('examples/demo_exchange_mismatch.m')` from the repository root, or pass the full script path to `run` from any MATLAB folder. The [population example](examples/demo_exchange_mismatch.m) samples permeable sticks, ECS and impermeable spheres with the same exchange time (τ = 10 ms) in every test tissue. NLLS and polynomial regression fit the **same noisy S₀ and S₂** using impermeable sticks, ECS and spheres. The stick fraction and diffusivity in each fit correspond to the forward tissue's `fPerm` and `Dperm`. Each NLLS fit starts at its tissue's shared ground-truth values. A separate quadratic regressor trains at each SNR on an independent population of impermeable tissues with the same acquisition and noise process.

The default protocol crosses b = [0, 1, 2, 5, 10, 15] ms/µm² with Δ = 10:10:80 ms using δ = 6 ms throughout. At each Δ it acquires one b₀ and 64 directions per nonzero b. The noise standard deviation per direction is 1 / b₀ SNR before normalization by the pooled noisy b₀ mean. The example uses 300 shared test tissues, one realization per tissue per SNR, b₀ SNRs of 50, 100, 250 and 500, and 6,000 training tissues. Edit `countsPerB`, `c.snr`, `c.nTissue`, `c.nRep`, `c.nTrain`, `c.testPrior` and `c.trainPrior` in the [shared settings](examples/exchange_population_config.m). Set `c.noise="rician"` for magnitude data. A scalar in `c.prior` fixes that parameter in the estimators and their polynomial training labels.

The example writes `protocol.png`, the paired results in `population.mat`, and four scatter figures under `results/exchange_mismatch/`. Each estimator has a figure for fractions and p₂ and a figure for diffusivities and soma radius. Rows show SNRs and columns show parameters. The dashed line is identity. The brown line joins the mean fits in ground-truth quantile bins. The orange band shows each mean ± one within-bin standard deviation of **fit minus truth**. Bias and SD in the corner use all trials in the panel, including flagged NLLS exits marked in red. Axes match across SNRs and estimators for each parameter.

To redraw the scatter plots from the saved fits, run the following from the repository root.

```matlab
addpath(pwd,'examples')
outputDir = fullfile(pwd,'results','exchange_mismatch');
load(fullfile(outputDir,'population.mat'),'population')
tmixplots.plot_exchange_population(population,outputDir);
```

For the original noiseless, single-tissue signal comparison, run [demo_exchange_single_tissue.m](examples/demo_exchange_single_tissue.m). Set `truth.tau = Inf` there for a no-exchange control.

### Matched high-SNR sanity check

Run `run('examples/demo_exchange_high_snr.m')` for the [b₀ SNR = 10⁶ sanity check](examples/demo_exchange_high_snr.m). Signal generation, NLLS and polynomial training all use permeable sticks exchanging with ECS plus impermeable spheres. The default forward τ is 10 ms for every test tissue. All ten fit parameters are free, including τ. The three fractions sum to one. The training population spans the fit range for τ (5 to 100 ms by default).

NLLS starts from each tissue's full ground truth. For every candidate parameter vector, it predicts directional signals using the simulated fibre axis, then applies the same finite-direction SH projection used to form the observations. Only the resulting S₀ and S₂ enter the residual. This synthetic check treats the fibre axis as known acquisition geometry relative to the tissue. The saved `nlls.start`, `nlls.trials.startSSE` and `nlls.trials.sse` allow comparison of the actual optimizer starts and the initial and final residuals. The four scatter figures and `population.mat` go to `results/exchange_high_snr_matched/`.

### Compare stick permeability and fixing both ECS diffusivities

Run `run('examples/demo_compare_exchange_assumptions.m')` or pass the full [script path](examples/demo_compare_exchange_assumptions.m) to MATLAB's `run`. Edit **`targetSNR`** at the top (default 100). The generator varies every tissue parameter, including both ECS diffusivities and exchange time. Set **`testTau`** to the desired uniform range in ms (default `[5 40]`). The default polynomial training law for τ uses the same range. Every estimator and assumption receives the same tissues, orientations, normalized S₀/S₂ measurements and noise draws. Both stick models use the same finite-direction projection in NLLS. NLLS starts at the feasible truth for each assumption. Polynomial training uses that assumption's kernel and constraints at the requested SNR.

| Stick model in fit | ECS diffusivities in fit |
|---|---|
| Permeable, with free exchange time | Both free |
| Impermeable | Both free |
| Permeable, with free exchange time | Both fixed at their population means |
| Impermeable | Both fixed at their population means |

The two **`fixedECS`** values are the means of the generating priors, 1.6 µm²/ms axial and 0.525 µm²/ms radial for the default uniform distributions. For custom sampling functions, enter their known means in this setting. The defaults use 300 tissues, four independent noise draws per tissue and 6,000 independent training examples per regressor. Keep at least two repeats to measure precision. The script can use four local workers for NLLS with Parallel Computing Toolbox. Set `c.nllsWorkers=0` for serial execution.

The output directory is `results/exchange_assumptions/snr_<SNR>/`. Open **`all_parameters.png`** for every ground-truth-versus-fit scatter panel. Its four rows group the fitting assumptions, with NLLS on the left and polynomial regression on the right. Each assumption subdirectory also contains `nlls_all_parameters.png` and `polynomial_all_parameters.png`, with two rows of parameters per figure. To redraw all of these from saved estimates, run the following from the repository root.

```matlab
addpath(pwd,'examples')
outputDir = fullfile(pwd,'results','exchange_assumptions','snr_100');
load(fullfile(outputDir,'comparison.mat'),'fits')
tmixplots.plot_exchange_assumptions(fits,outputDir,1);
```

The default scatter plots show one paired noise draw per tissue. Pass `[]` in place of `1` to show every draw. The ground truths vary along the horizontal axis even when the estimator fixes that parameter. This produces horizontal bands at the prescribed ECS values. The permeable fits include the varying τ truth and its fitted value. The bias and SD printed inside these scatter panels describe the displayed errors across tissues. **`bias_precision.png` and `bias_precision.csv` use all repeats**, including flagged NLLS fits. Bars show ±1 repeat SD, computed as the square root of the mean of the per-tissue sample variances. A signed bias of one repeat SD means that the average error is as large as a typical noise fluctuation.

The CSV gives signed bias, across-population error SD, within-tissue repeat SD, RMSE and nonconvergence counts. `biasOverRepeatSD` uses each fit's own precision. `biasShiftOverMatchedSD` measures the change in bias from the permeable, free-ECS fit in units of that estimator's matched repeat SD. Fixed parameters have zero repeat SD and an undefined own-SD ratio. Their RMSE still measures errors across tissues. The impermeable fits have no exchange-time estimate. `comparison.mat` contains per-trial truths, estimates, observations, exit flags and both regressor and fit configurations.

### Compare free, fixed and very slow exchange with isotropic ECS

Run [demo_compare_sandix_exchange.m](examples/demo_compare_sandix_exchange.m) to simulate permeable sticks exchanging with Gaussian isotropic ECS plus impermeable spheres. Each tissue has one ECS diffusivity, `De`. The default acquisition uses b = [0, 1, 2, 5, 10, 15] ms/µm², Δ = 10:10:80 ms, δ = 6 ms, 30 directions per diffusion-weighted shell and one b₀ image per Δ. Per-image b₀ SNR is 50. Independent Gaussian noise is added before pooled b₀ normalization and S₀/S₂ extraction. Set `c.noise="rician"` for magnitude noise.

The same 300 tissues and four noise draws are fitted by both NLLS and quadratic polynomial regression under three assumptions.

| Fit | Exchange time | Independent free parameters |
|---|---|---|
| SANDIX-like | Estimated | 8 |
| SANDIX-fixed | Fixed at 10 ms | 7 |
| SANDI-like | Fixed at 1,000,000 ms | 7 |

The nine displayed parameters are three fractions (which sum to one), stick diffusivity, ECS diffusivity, sphere radius, intrasomal diffusivity, fODF p₂ and exchange time. All ground-truth parameters vary across tissues in every fit. Ground-truth τ is uniform from 5 to 30 ms. The default free-τ fit bounds and training distribution also span 5 to 30 ms. Edit `testTau`, `c.prior.tau`, and `c.trainPrior.tau` independently. The default ECS diffusivity is uniform from 0.5 to 1.5 µm²/ms. Other distribution settings are in `c.testPrior` and `c.trainPrior`.

Set **`fixedTau`** to the list of assumed values to compare. For example, `fixedTau = [5 10 15 20 25 30 1e6]` runs the free-τ baseline and each listed fixed-τ fit on identical noisy observations. Each fixed fit imposes its τ on the polynomial training population while the test population retains its varying τ. NLLS starts at each tissue's ground truth for all free parameters. It evaluates candidate S₀/S₂ through the same finite-direction projection as the simulated data, using the known simulated fibre axis. The fixed value 1e6 ms is evaluated by the finite-pulse exchange kernel and approximates impermeable sticks. τ follows TMImix's inverse-sum-of-rates convention described below.

Results go to `results/sandix_exchange/snr_<SNR>/`. **`all_parameters.png`** contains all GT-versus-fit scatter panels, with one assumption per block row and NLLS and polynomial regression in its two columns. Each assumption subdirectory contains a 3×3 figure per estimator. The fixed-1e6-ms τ scatter uses separate x and y ranges to show both varying truths and the prescribed fitted value. **`bias_precision.png`**, **`rmse.png`**, and **`repeat_precision.png`** compare signed bias, total error and within-tissue noise SD across all repeats. Fixed outputs have zero repeat variance and can have large errors across tissues. The CSV includes `rmseRatio` and `precisionRatio`, each relative to that estimator's free-τ baseline. A ratio below one indicates a reduction. Flagged NLLS exits are included in the metrics and marked in red on the scatter plots.

```matlab
% Redraw saved estimates from the repository root.
addpath(pwd,'examples')
outputDir = fullfile(pwd,'results','sandix_exchange','snr_50');
load(fullfile(outputDir,'comparison.mat'),'fits','metrics')
tmixplots.plot_exchange_assumptions(fits,outputDir,1);
tmixplots.plot_bias_precision(metrics,outputDir);
tmixplots.plot_bias_precision(metrics,outputDir,"RMSE");
tmixplots.plot_bias_precision(metrics,outputDir,"withinTissueSD");
```

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

### Plot noise propagation across all four modes

Run `run('examples/demo_compare_estimators.m')` from the repository root. The [example script](examples/demo_compare_estimators.m) fits permeable sticks, ECS, and impermeable spheres with both NLLS and polynomial regression. At the top, set **`matchedFixed`** and **`mismatchedFixed`** independently. Each accepts any active parameter names and physical values, including multiple fields such as `struct('radius',6,'Dsoma',2)`. Set `oracleKnown` separately to the names whose true values are supplied for each tissue. The same settings section defines the protocol, b=0 SNR sweep, priors, and sample counts.

The script writes `results/compare_estimators/bias_rmse_by_parameter.pdf` and `metrics.csv`. Each PDF page shows one parameter in its native units, with one column per mode. Each point pools the held-out tissues and noise repetitions at one SNR. Bias is the mean of estimate minus truth. The CSV includes the per-mode constraints, sample counts and NLLS nonconvergence counts. The individual estimates and truths remain available in the MATLAB workspace, for example in `comparison.mismatched.nlls`.

The final script section writes `results/compare_estimators/ground_truth_vs_fit.pdf` directly from `comparison`. Rerun that section to redraw the scatter plots from existing fits. Each parameter has ground truth on x and the fitted value on y, with equal scales and a dashed identity line. Columns show the four modes and rows show b=0 SNR, with up to four SNRs per page. Each marker represents one tissue and one noise realization. Repeated or exactly constrained fits can coincide. Blue crosses identify NLLS estimates with a nonpositive exit flag.

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

For Gaussian isotropic ECS, use `TMImix(protocol,parts,ECS="isotropic")`. This selects a single `De` parameter for the ECS response, including its role in exchange. In the noise examples, `noise_config(parts,protocol,"isotropic")` creates the corresponding fit schema. The default `ECS="axial"` uses separate `DePar` and `DePerp`.

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
