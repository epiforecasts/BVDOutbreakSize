```@meta
EditURL = "../examples/province.jl"
```

# Province estimates

The outbreak by province.
The model runs one renewal equation per province and fits the national
streams against the summed provinces, so the national estimates on the
[analysis](@ref "Results") page are the sum of the provinces here.

The per-province case-fatality ratio sits with the national one in the
[confirmed case-fatality ratio](@ref "Confirmed case-fatality ratio").
The per-province forecast is in the
[one-week-ahead forecast results](@ref "One-week-ahead forecast results")
and its scoring in the [forecast by province](@ref "Forecast by province").
Each needs a quantity the page it sits on already computes.

```@raw html
<details><summary>Load packages, data and fitted chains</summary>
```

````julia
# Shared setup: packages, observations, the fit registry and every model fit
# (loaded from the content-addressed cache). See `docs/examples/_setup.jl`.
using BVDOutbreakSize
include(joinpath(pkgdir(BVDOutbreakSize), "docs", "examples", "_setup.jl"))
````

````
_release_data (generic function with 1 method)
````

```@raw html
</details>
```

## Size and infections

The national outbreak size on the [analysis](@ref "Joint model estimates") page is the sum of the four patches' renewal equations.
The reproduction number and the relative case ascertainment are identified only as a product, and the per-province deaths break the tie.

```@raw html
<details><summary>Cross-province overview table</summary>
```

````julia
province_overview_table = patch_overview_table(chn_joint, N_PATCHES);
````

```@raw html
</details>
```


```@raw html
<div><div style = "float: left;"><span>4×5 DataFrame</span></div><div style = "clear: both;"></div></div><div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">Province</th><th style = "text-align: left;">Reproduction number</th><th style = "text-align: left;">Cumulative infections</th><th style = "text-align: left;">Share of infections (%)</th><th style = "text-align: left;">Relative ascertainment</th></tr><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;"></th><th title = "String" style = "text-align: left;">String</th><th title = "String" style = "text-align: left;">String</th><th title = "String" style = "text-align: left;">String</th><th title = "String" style = "text-align: left;">String</th><th title = "String" style = "text-align: left;">String</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">Ituri</td><td style = "text-align: left;">1.02 (0.86–1.29)</td><td style = "text-align: left;">9483 (6860–13960)</td><td style = "text-align: left;">69.9 (60.8–76.9)</td><td style = "text-align: left;">1.19 (0.81–1.94)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">Nord-Kivu</td><td style = "text-align: left;">1.1 (0.92–1.4)</td><td style = "text-align: left;">3325 (2250–5438)</td><td style = "text-align: left;">24.6 (18.7–32.2)</td><td style = "text-align: left;">1.0 (0.72–1.33)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">3</td><td style = "text-align: left;">Haut-Uele</td><td style = "text-align: left;">0.97 (0.71–1.3)</td><td style = "text-align: left;">575 (351–1034)</td><td style = "text-align: left;">4.2 (2.8–6.7)</td><td style = "text-align: left;">1.01 (0.81–1.3)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">4</td><td style = "text-align: left;">Other provinces</td><td style = "text-align: left;">1.08 (0.84–1.43)</td><td style = "text-align: left;">123 (52–327)</td><td style = "text-align: left;">0.9 (0.4–2.3)</td><td style = "text-align: left;">0.85 (0.43–1.43)</td></tr></tbody></table></div>
```

The figure below shows the modelled infections behind those totals, daily on the top row and cumulative on the bottom.

```@raw html
<details><summary>Modelled infections by province</summary>
```

````julia
province_infections_fig = plot_infections_patches(
    chn_joint;
    n = obs.n, seeding = obs.seeding, n_patches = N_PATCHES
);
````

```@raw html
</details>
```

![](province-16.png)

The provinces are coupled by a gravity kernel weighted by destination population, described in the [mixing and importation](@ref "Mixing and importation") Methods section, with its intensity estimated.
Every arrival is debited from its origin the same day, so the figure reads as where infection occurred rather than as extra infection.
The distances between the patch capitals are 379 km from Bunia to Goma, 322 km from Bunia to Isiro and 206 km from Goma to the pooled patch's centre, so most of what leaves Nord-Kivu lands in the pooled patch.

```@raw html
<details><summary>Importation intensity and imports by province</summary>
```

````julia
importation_table = summary_table(
    chn_joint, [:importation_epsilon];
    digits = 4,
    labels = Dict(:importation_epsilon => "Importation intensity")
);

province_imports_fig = plot_imports_patches(
    chn_joint;
    n = obs.n, seeding = obs.seeding, n_patches = N_PATCHES
);
````

```@raw html
</details>
```

````julia

````
![](province-21.png)

## Reproduction number by province

The national reproduction number is on the [analysis](@ref "Reproduction number over time") page.
Below it is split by province, one panel per province with the national trajectory in grey behind it.
The deviations sum to zero, so the grey band is the incidence-weighted middle of the panels rather than any one province.
A panel tracking grey says that province moves with the national trajectory.
The pooled patch holds almost no confirmed cases, so its panel is carried by the deviation prior and its width is not a measurement.

```@raw html
<details><summary>Reproduction number by province</summary>
```

````julia
province_rt_fig = plot_rt_patches(
    chn_joint;
    n = obs.n, breakpoint = _BREAKPOINT,
    n_patches = N_PATCHES,
    rt_start = _rt_start_plot,
    rt_walk_start = clamp(_BREAKPOINT - RT_WALK_LEAD, _rt_start_plot, obs.n),
    display_start = _rt_start_plot,
    as_of_date = string(obs.cutoff), seeding = obs.seeding,
    ramp = RT_INTERVENTION_RAMP
);
````

```@raw html
</details>
```

![](province-27.png)

The figure below gives each province's estimates, including its log-Rt deviation from the trend, the deviation's walk scale and the contrast against Ituri.

```@raw html
<details><summary>Per-province summary figure</summary>
```

````julia
province_detail_fig = plot_patch_summary(chn_joint, N_PATCHES);
````

```@raw html
</details>
```

![](province-32.png)

```@raw html
<details><summary>Per-province summary table</summary>
```

````julia
province_detail_table = patch_summary_table(chn_joint, N_PATCHES);

````

```@raw html
<div><div style = "float: left;"><span>28×8 DataFrame</span></div><div style = "clear: both;"></div></div><div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">Patch</th><th style = "text-align: left;">Quantity</th><th style = "text-align: left;">Lower 90%</th><th style = "text-align: left;">Lower 60%</th><th style = "text-align: left;">Lower 30%</th><th style = "text-align: left;">Upper 30%</th><th style = "text-align: left;">Upper 60%</th><th style = "text-align: left;">Upper 90%</th></tr><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;"></th><th title = "String" style = "text-align: left;">String</th><th title = "String" style = "text-align: left;">String</th><th title = "Float64" style = "text-align: left;">Float64</th><th title = "Float64" style = "text-align: left;">Float64</th><th title = "Float64" style = "text-align: left;">Float64</th><th title = "Float64" style = "text-align: left;">Float64</th><th title = "Float64" style = "text-align: left;">Float64</th><th title = "Float64" style = "text-align: left;">Float64</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">Ituri</td><td style = "text-align: left;">Cumulative infections</td><td style = "text-align: right;">6860.0</td><td style = "text-align: right;">7830.0</td><td style = "text-align: right;">8630.0</td><td style = "text-align: right;">10408.0</td><td style = "text-align: right;">11467.0</td><td style = "text-align: right;">13960.0</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">Ituri</td><td style = "text-align: left;">Reproduction number</td><td style = "text-align: right;">0.86</td><td style = "text-align: right;">0.94</td><td style = "text-align: right;">0.99</td><td style = "text-align: right;">1.07</td><td style = "text-align: right;">1.13</td><td style = "text-align: right;">1.29</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">3</td><td style = "text-align: left;">Ituri</td><td style = "text-align: left;">Daily infections at cut-off</td><td style = "text-align: right;">34.0</td><td style = "text-align: right;">51.0</td><td style = "text-align: right;">60.0</td><td style = "text-align: right;">86.0</td><td style = "text-align: right;">108.0</td><td style = "text-align: right;">164.0</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">4</td><td style = "text-align: left;">Ituri</td><td style = "text-align: left;">log-Rt deviation from trend</td><td style = "text-align: right;">-0.13</td><td style = "text-align: right;">-0.07</td><td style = "text-align: right;">-0.04</td><td style = "text-align: right;">0.01</td><td style = "text-align: right;">0.05</td><td style = "text-align: right;">0.11</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">5</td><td style = "text-align: left;">Ituri</td><td style = "text-align: left;">log-Rt vs primary patch</td><td style = "text-align: right;">0.0</td><td style = "text-align: right;">0.0</td><td style = "text-align: right;">0.0</td><td style = "text-align: right;">0.0</td><td style = "text-align: right;">0.0</td><td style = "text-align: right;">0.0</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">6</td><td style = "text-align: left;">Ituri</td><td style = "text-align: left;">Rt deviation drift</td><td style = "text-align: right;">0.004</td><td style = "text-align: right;">0.019</td><td style = "text-align: right;">0.031</td><td style = "text-align: right;">0.05</td><td style = "text-align: right;">0.062</td><td style = "text-align: right;">0.089</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">7</td><td style = "text-align: left;">Ituri</td><td style = "text-align: left;">Relative case ascertainment</td><td style = "text-align: right;">0.81</td><td style = "text-align: right;">0.98</td><td style = "text-align: right;">1.09</td><td style = "text-align: right;">1.3</td><td style = "text-align: right;">1.48</td><td style = "text-align: right;">1.94</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">8</td><td style = "text-align: left;">Nord-Kivu</td><td style = "text-align: left;">Cumulative infections</td><td style = "text-align: right;">2250.0</td><td style = "text-align: right;">2722.0</td><td style = "text-align: right;">3021.0</td><td style = "text-align: right;">3737.0</td><td style = "text-align: right;">4198.0</td><td style = "text-align: right;">5438.0</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">9</td><td style = "text-align: left;">Nord-Kivu</td><td style = "text-align: left;">Reproduction number</td><td style = "text-align: right;">0.92</td><td style = "text-align: right;">1.01</td><td style = "text-align: right;">1.05</td><td style = "text-align: right;">1.15</td><td style = "text-align: right;">1.23</td><td style = "text-align: right;">1.4</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">10</td><td style = "text-align: left;">Nord-Kivu</td><td style = "text-align: left;">Daily infections at cut-off</td><td style = "text-align: right;">40.0</td><td style = "text-align: right;">56.0</td><td style = "text-align: right;">69.0</td><td style = "text-align: right;">95.0</td><td style = "text-align: right;">118.0</td><td style = "text-align: right;">174.0</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">11</td><td style = "text-align: left;">Nord-Kivu</td><td style = "text-align: left;">log-Rt deviation from trend</td><td style = "text-align: right;">-0.07</td><td style = "text-align: right;">-0.0</td><td style = "text-align: right;">0.03</td><td style = "text-align: right;">0.08</td><td style = "text-align: right;">0.12</td><td style = "text-align: right;">0.2</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">12</td><td style = "text-align: left;">Nord-Kivu</td><td style = "text-align: left;">log-Rt vs primary patch</td><td style = "text-align: right;">-0.1</td><td style = "text-align: right;">-0.0</td><td style = "text-align: right;">0.03</td><td style = "text-align: right;">0.11</td><td style = "text-align: right;">0.15</td><td style = "text-align: right;">0.25</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">13</td><td style = "text-align: left;">Nord-Kivu</td><td style = "text-align: left;">Rt deviation drift</td><td style = "text-align: right;">0.002</td><td style = "text-align: right;">0.009</td><td style = "text-align: right;">0.017</td><td style = "text-align: right;">0.038</td><td style = "text-align: right;">0.049</td><td style = "text-align: right;">0.078</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">14</td><td style = "text-align: left;">Nord-Kivu</td><td style = "text-align: left;">Relative case ascertainment</td><td style = "text-align: right;">0.72</td><td style = "text-align: right;">0.85</td><td style = "text-align: right;">0.94</td><td style = "text-align: right;">1.06</td><td style = "text-align: right;">1.14</td><td style = "text-align: right;">1.33</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">15</td><td style = "text-align: left;">Haut-Uele</td><td style = "text-align: left;">Cumulative infections</td><td style = "text-align: right;">351.0</td><td style = "text-align: right;">442.0</td><td style = "text-align: right;">513.0</td><td style = "text-align: right;">651.0</td><td style = "text-align: right;">751.0</td><td style = "text-align: right;">1034.0</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">16</td><td style = "text-align: left;">Haut-Uele</td><td style = "text-align: left;">Reproduction number</td><td style = "text-align: right;">0.71</td><td style = "text-align: right;">0.84</td><td style = "text-align: right;">0.9</td><td style = "text-align: right;">1.03</td><td style = "text-align: right;">1.12</td><td style = "text-align: right;">1.3</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">17</td><td style = "text-align: left;">Haut-Uele</td><td style = "text-align: left;">Daily infections at cut-off</td><td style = "text-align: right;">2.0</td><td style = "text-align: right;">3.0</td><td style = "text-align: right;">4.0</td><td style = "text-align: right;">7.0</td><td style = "text-align: right;">10.0</td><td style = "text-align: right;">18.0</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">18</td><td style = "text-align: left;">Haut-Uele</td><td style = "text-align: left;">log-Rt deviation from trend</td><td style = "text-align: right;">-0.32</td><td style = "text-align: right;">-0.19</td><td style = "text-align: right;">-0.13</td><td style = "text-align: right;">-0.02</td><td style = "text-align: right;">0.04</td><td style = "text-align: right;">0.13</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">19</td><td style = "text-align: left;">Haut-Uele</td><td style = "text-align: left;">log-Rt vs primary patch</td><td style = "text-align: right;">-0.39</td><td style = "text-align: right;">-0.21</td><td style = "text-align: right;">-0.12</td><td style = "text-align: right;">0.01</td><td style = "text-align: right;">0.08</td><td style = "text-align: right;">0.21</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">20</td><td style = "text-align: left;">Haut-Uele</td><td style = "text-align: left;">Rt deviation drift</td><td style = "text-align: right;">0.057</td><td style = "text-align: right;">0.078</td><td style = "text-align: right;">0.09</td><td style = "text-align: right;">0.113</td><td style = "text-align: right;">0.127</td><td style = "text-align: right;">0.158</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">21</td><td style = "text-align: left;">Haut-Uele</td><td style = "text-align: left;">Relative case ascertainment</td><td style = "text-align: right;">0.81</td><td style = "text-align: right;">0.93</td><td style = "text-align: right;">0.97</td><td style = "text-align: right;">1.04</td><td style = "text-align: right;">1.12</td><td style = "text-align: right;">1.3</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">22</td><td style = "text-align: left;">Other provinces</td><td style = "text-align: left;">Cumulative infections</td><td style = "text-align: right;">52.0</td><td style = "text-align: right;">80.0</td><td style = "text-align: right;">100.0</td><td style = "text-align: right;">152.0</td><td style = "text-align: right;">195.0</td><td style = "text-align: right;">327.0</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">23</td><td style = "text-align: left;">Other provinces</td><td style = "text-align: left;">Reproduction number</td><td style = "text-align: right;">0.84</td><td style = "text-align: right;">0.96</td><td style = "text-align: right;">1.02</td><td style = "text-align: right;">1.14</td><td style = "text-align: right;">1.22</td><td style = "text-align: right;">1.43</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">24</td><td style = "text-align: left;">Other provinces</td><td style = "text-align: left;">Daily infections at cut-off</td><td style = "text-align: right;">1.0</td><td style = "text-align: right;">2.0</td><td style = "text-align: right;">3.0</td><td style = "text-align: right;">5.0</td><td style = "text-align: right;">7.0</td><td style = "text-align: right;">12.0</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">25</td><td style = "text-align: left;">Other provinces</td><td style = "text-align: left;">log-Rt deviation from trend</td><td style = "text-align: right;">-0.14</td><td style = "text-align: right;">-0.04</td><td style = "text-align: right;">0.0</td><td style = "text-align: right;">0.06</td><td style = "text-align: right;">0.12</td><td style = "text-align: right;">0.21</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">26</td><td style = "text-align: left;">Other provinces</td><td style = "text-align: left;">log-Rt vs primary patch</td><td style = "text-align: right;">-0.2</td><td style = "text-align: right;">-0.06</td><td style = "text-align: right;">0.0</td><td style = "text-align: right;">0.1</td><td style = "text-align: right;">0.16</td><td style = "text-align: right;">0.29</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">27</td><td style = "text-align: left;">Other provinces</td><td style = "text-align: left;">Rt deviation drift</td><td style = "text-align: right;">0.007</td><td style = "text-align: right;">0.02</td><td style = "text-align: right;">0.036</td><td style = "text-align: right;">0.066</td><td style = "text-align: right;">0.085</td><td style = "text-align: right;">0.115</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">28</td><td style = "text-align: left;">Other provinces</td><td style = "text-align: left;">Relative case ascertainment</td><td style = "text-align: right;">0.43</td><td style = "text-align: right;">0.63</td><td style = "text-align: right;">0.74</td><td style = "text-align: right;">0.94</td><td style = "text-align: right;">1.09</td><td style = "text-align: right;">1.43</td></tr></tbody></table></div>
```

```@raw html
</details>
```

The spread of those deviations is the spatial diagnostic.
The prior admits real divergence, with a 31% prior probability that the Ituri to Nord-Kivu ratio moves by more than 25% over the window, so a shrunken posterior is a finding rather than an artefact of the prior.
With four patches and the pooled one carrying almost no signal the cross-province correlation is not identified, and it tracks its prior.

```@raw html
<details><summary>Spatial hyperparameter summary table</summary>
```

````julia
spatial_hyper_table = summary_table(
    chn_joint,
    [
        :region_sd, :region_halflife, :region_corr_primary_secondary,
        :province_ascertainment_sd, :province_testing_coefficient,
    ];
    digits = 3,
    labels = Dict(
        :region_sd => "Rt deviation spread",
        :region_halflife => "Rt deviation half-life (days)",
        :region_corr_primary_secondary => "Ituri-N.Kivu Rt correlation",
        :province_ascertainment_sd => "Ascertainment spread",
        :province_testing_coefficient => "Testing effect on ascertainment"
    )
);
````

```@raw html
</details>
```


```@raw html
<div><div style = "float: left;"><span>5×7 DataFrame</span></div><div style = "clear: both;"></div></div><div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">Quantity</th><th style = "text-align: left;">Lower 90%</th><th style = "text-align: left;">Lower 60%</th><th style = "text-align: left;">Lower 30%</th><th style = "text-align: left;">Upper 30%</th><th style = "text-align: left;">Upper 60%</th><th style = "text-align: left;">Upper 90%</th></tr><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;"></th><th title = "String" style = "text-align: left;">String</th><th title = "Float64" style = "text-align: left;">Float64</th><th title = "Float64" style = "text-align: left;">Float64</th><th title = "Float64" style = "text-align: left;">Float64</th><th title = "Float64" style = "text-align: left;">Float64</th><th title = "Float64" style = "text-align: left;">Float64</th><th title = "Float64" style = "text-align: left;">Float64</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">Rt deviation spread</td><td style = "text-align: right;">0.107</td><td style = "text-align: right;">0.153</td><td style = "text-align: right;">0.186</td><td style = "text-align: right;">0.247</td><td style = "text-align: right;">0.288</td><td style = "text-align: right;">0.373</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">Rt deviation half-life (days)</td><td style = "text-align: right;">19.571</td><td style = "text-align: right;">26.713</td><td style = "text-align: right;">32.392</td><td style = "text-align: right;">44.117</td><td style = "text-align: right;">55.02</td><td style = "text-align: right;">80.409</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">3</td><td style = "text-align: left;">Ituri-N.Kivu Rt correlation</td><td style = "text-align: right;">-0.565</td><td style = "text-align: right;">-0.351</td><td style = "text-align: right;">-0.178</td><td style = "text-align: right;">0.121</td><td style = "text-align: right;">0.294</td><td style = "text-align: right;">0.536</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">4</td><td style = "text-align: left;">Ascertainment spread</td><td style = "text-align: right;">0.015</td><td style = "text-align: right;">0.062</td><td style = "text-align: right;">0.11</td><td style = "text-align: right;">0.222</td><td style = "text-align: right;">0.309</td><td style = "text-align: right;">0.494</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">5</td><td style = "text-align: left;">Testing effect on ascertainment</td><td style = "text-align: right;">-0.134</td><td style = "text-align: right;">-0.027</td><td style = "text-align: right;">0.029</td><td style = "text-align: right;">0.112</td><td style = "text-align: right;">0.178</td><td style = "text-align: right;">0.313</td></tr></tbody></table></div>
```

## [Province compositions](@id province-compositions)

The per-province confirmed cases and deaths are fitted as compositions conditional on the national total, so what the model predicts is each province's share rather than its count.
The panels below show that modelled share at every spatial vintage against the observed one.
Each panel carries two bands.
The grey band is the posterior predictive interval on the observed share, built by pushing every posterior draw's expected shares back through the composition's own overdispersed allocation at that vintage's observed total.
The overdispersion is what absorbs reporting lags between the provincial and national tables and the reassignment of cases between health zones.
The observed points should fall inside it.
The coloured ribbon inside the grey band is the expected share alone, which is the modelled centre the points scatter around.
A point outside the grey band is a vintage the composition does not reproduce, and points consistently to one side of the coloured ribbon are a province the model splits wrongly on average.
Each panel starts at zero and takes its own upper limit, because the shares differ by orders of magnitude.
The vintages stop before the cut-off, so the panels end earlier than the [national posterior predictive checks](@ref "Posterior predictive checks").

```@raw html
<details><summary>Province composition posterior predictive checks</summary>
```

````julia
province_case_ppc_fig = plot_province_composition_ppc(
    chn_joint;
    share_key = :province_shares,
    obs_increments = province_cases.increments,
    days = province_cases.days, seeding = obs.seeding, n_patches = N_PATCHES,
    title = "Confirmed case share by province"
);

province_death_ppc_fig = plot_province_composition_ppc(
    chn_joint;
    share_key = :province_death_shares,
    obs_increments = province_deaths.increments,
    days = province_deaths.days, seeding = obs.seeding, n_patches = N_PATCHES,
    title = "Confirmed death share by province"
);
````

```@raw html
</details>
```

````julia

````
![](province-45.png)

## Saving province assets

The summary dashboard shows the per-province table and the reproduction
number by province, so they are written here rather than on the analysis
page.

```@raw html
<details><summary>Write the province dashboard assets</summary>
```

````julia
dashboard_dir = joinpath(
    pkgdir(BVDOutbreakSize), "docs", "src", "summary_assets"
)
mkpath(dashboard_dir)
CairoMakie.save(joinpath(dashboard_dir, "rt_provinces.png"), province_rt_fig)
open(joinpath(dashboard_dir, "provinces.md"), "w") do io
    print(io, markdown_table(province_overview_table))
end
````

```@raw html
</details>
```

---

The full analysis code, data and model definitions are in the
[epiforecasts/BVDOutbreakSize](https://github.com/epiforecasts/BVDOutbreakSize)
repository.

