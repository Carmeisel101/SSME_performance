# Space Shuttle Main Engine (SSME) perfomance model

This is a simple model of the Space Shuttle Main Engine (SSME) performance. It is based on the data from the Dr. Rodney Bowersox:

<div align="center">


| Combustion Chamber Data | |
|--------------------------| ----|
|Fuel|  Hydrogen |
|Oxidizer| O2 |
|Fuel/Oxidizer Ratio | 0.166 |
| Hydrogen Injection Temperature |  850 K |
| Oxygen Injection Temperature |  530 K |
| Combustion Chamber Pressure |  204 atm |

| Nozzle Geometry | |
|--------------------------| ----|
| $R_{t}$ |  5.15 in |
| $R_{e}$ |  45.35 in |
| $\theta_{p}$ |  32.00 deg |
</div>


## Oxygen First Enthalpy Exchanger
We begin with computing the equivalence ratio $\phi$ for the oxygen first enthalpy exchanger. The equivalence ratio is defined as the ratio of the mass of fuel to the mass of oxidizer. The mass of fuel is the mass of hydrogen, and the mass of oxidizer is the mass of oxygen. The mass of hydrogen is $\phi \frac{mol_{O_{2}}}{mol_{H_{2}}} = \Phi$

From here we are to calculate the standard heat of formation by:

<div align="center">

$\Delta \hat{h_{i}} =  \hat{h_{H_{2}}}(850K) -  \hat{h_{H_{2}}}(298K) + \hat{R}T^2 \int_{p_{f}}^{p_{i}} \frac{1}{p}  \frac{-\hat{b} p}{\hat{R}T^2} dp|_{T_i}$

</div>

Where $\hat{b}$ is read from a table developed in `exchanger.py`. $\Delta \hat{h_{i}}$ is calculated for both hydrogen and oxygen. Once we calculate $\Delta \hat{h_{i}}$ for each, we then calculate the enthalpy of the mixture:

<div align="center">

$\Delta \hat{h_{mix}} = \Phi \Delta \hat{h_{H_{2}}} + \Delta \hat{h_{O_{2}}}$

</div>

Which is then applied to our chemical equation:

<div align="center">

$O_{2} + \Phi H_{2} \rightarrow (2-\Phi)H_{2} + 2H_{2}O$

</div>

Note: notice the $(2-\Phi)H_{2}$ portion of the equation. This is used since we calculated the enthalpy of formation for hydrogen (H) and not hydrogen gas $H_{2}$.

From here we can sum the molar enthalpies of formation to get the first estimate of the first reaction enthalpy:

<div align="center">

$\Delta \hat{h_{r}} = \Phi \Delta \hat{h_{H_{2}}} + 2( \Delta \hat{h_{H_{2}O}} )$

</div>

We are to then interpolate this value from the given turn tables to get the actual value of the combustion temperature for Oxygen First Enthalpy Exchanger. This is done in `exchanger.py`.