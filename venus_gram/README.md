# Venus-GRAM flyby

This folder adds atmospheric propagation and a post-flyby Earth correction
to the Earth–Venus–Earth project. The runner calls
`AA279lambert_curtis`, `AA279j2000_planetary_elements`, and `rot`
directly from the adjacent `Pearl's Code` folder.

## Run the example

From the repository root in MATLAB:

```matlab
addpath('venus_gram');
mission = run_venus_gram_flyby(Beta=100);
plot_venus_gram_flyby(mission);
assertSuccess(runtests('venus_gram/test_venus_gram_flyby.m'));
```

The example launches on 24 December 2024, passes Venus on 25 April 2025,
and returns to Earth on 6 January 2026. Set `Beta` to the spacecraft's
`mass / (drag coefficient * reference area)` in kg/m².

## Use a candidate from Pearl's search

```matlab
addpath('venus_gram');
addpath("Pearl's Code");
AA279B_HM_claude_lowfly

reported = cand(idx_report,:);
insideBand = reported.low == 1 & reported.h_p > 80 & reported.h_p < 300;
lowCandidates = reported(insideBand,:);
disp(lowCandidates(:,{'dep_MJD','fb_MJD','ret_MJD','h_p'}));
```

Choose a row from `lowCandidates` and pass it to the runner:

```matlab
mission = run_venus_gram_flyby(PearlCandidate=lowCandidates(1,:), Beta=100);
plot_venus_gram_flyby(mission);
```

The runner converts the row's Modified Julian Dates and reconstructs the
inbound aim from the two Lambert legs. It checks the recomputed vacuum
periapsis against the altitude floor and interface. You can also supply
`DatesJD=[departureJD encounterJD returnJD]` directly. The full Pearl
search uses `wrapTo360` from Mapping Toolbox in addition to the
functions shipped in this repository.

## Set up the local Windows Venus-GRAM runtime

The runner looks for `venus_gram/runtime` by default. To point it to your
own GRAM Suite installation, set the root once in MATLAB:

```matlab
setenv('VENUS_GRAM_ROOT','C:/path/to/GRAM');
```

You can also pass `GramRoot="C:/path/to/GRAM"` on an individual run.
The root is the folder containing these files:

```text
Windows/VenusGRAM.exe
SPICE/lsk/naif0012.tls
SPICE/pck/pck00011.tpc
SPICE/spk/planets/de440_GRAM.bsp
```

The Git ignore rule keeps the local runtime outside version control.
Each collaborator can use their own NASA GRAM distribution.

## Read the results

- `DragBallistic` contains the unpowered atmospheric pass and its closest
  approach, peak dynamic pressure, and exit status. The model uses a mean
  Venus-GRAM density and wind profile at the nominal periapsis location.
- `DragPatch` contains the correction vector, the burn epoch, and the
  propagated Earth return. The default burn location is the outbound
  crossing of 616,000 km from Venus's center.
- `VacuumBallistic` and `VacuumPatch` provide gravity-only comparisons from
  the same inbound state. When `DragPatch.Status` is `"targeted"`, the
  modeled propulsive budget is
  `mission.EarthInjection_kms + norm(mission.DragPatch.Burn_kms)`.
- Pearl's `dV_flyby` is part of the vacuum candidate ranking. The GRAM
  mission uses the computed post-flyby correction for its Venus maneuver.
