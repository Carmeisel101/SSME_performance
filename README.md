# Space Shuttle Main Engine (SSME) perfomance model

This is a simple model of the Space Shuttle Main Engine (SSME) performance. It is based on the data from the Dr. Rodney Bowersox:

<div align="center">


| Combustion Chamber Data | |
|--------------------------| ----|
|Fuel|  Hydrogen |
|Oxidizer| $O_{2}$ |
|Fuel/Oxidizer Ratio $\phi$ | 0.166 |
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

We are to then interpolate this value from the given turn tables to get the actual value of the combustion temperature for Oxygen First Enthalpy Exchanger. This is done in `exchanger.py`. At the end of this process we have the combustion temperature for the Oxygen First Enthalpy Exchanger, $T_{c}$.

## Oxygen - Combustion Chamber

For the Combustion Chamber we wil using the two sub-reactions:

<div align="center">

$H_{2}O \rightarrow OH+H$

$H_{2}\rightarrow 2H$

</div>

From here we are to calculate the Gibbs Free Energy of each species.

For $H_{2}\rightarrow 2H$

<div align="center">

$\Delta G_{H_{2}} = 2(\bar{g_{H}} -\bar{g_{H_{2}}})$

</div>

For $H_{2}O \rightarrow OH+H$

<div align="center">

$\Delta G_{H_{2}O} = \bar{g_{H}} + \bar{g_{OH}} - \bar{g_{H_{2}O}}$

</div>

Where $\bar{g_{fh}}$ is the Gibbs Free Energy of formation of $H_{2}$, $\bar{g_{foh}}$ is the Gibbs Free Energy of formation of OH, and $\bar{g_{fh_{2}o}}_{fh_{2}o}$ is the Gibbs Free Energy of formation of $H_{2}O$. These values are read from the `Iterater.py` file.

### Law of Mass Action
Looking at For $H_{2}O \rightarrow OH+H$:

<div align="center">

| Species | $\nu_{(i)} '$ | $\nu_{(i)} "$ | $\nu_{(i)} " - \nu_{(i)} '$ |
|---------|-------|-------|---------------|
| $H_{2}O$ | 1 | 0 | -1 |
| $OH$ | 0 | 1 | 1 |
| $H$ | 0 | 1 | 1 |

</div>

We know that:

<div align="center">

$K_{n} = \prod_{n=1}^{N_{s}} X_{(n)}^{(\nu_{(i)} " - \nu_{(i)} ')} = (\frac{p}{p^{ * }})^{- \sigma_{v}} K_{p}$

</div>

Where; 

<div align="center">

$\sigma_{v} = \sum_{n=1}^{N_{s}} \nu_{(i)} " - \nu_{(i)} '$

$K_{p} = e^{- \frac{\Delta G}{\hat{R}T}}$

$\Delta G = \sum_{n=1}^{N_{s}} (\nu_{(i)} " -\nu_{(i)}') \bar{g}_{f(n)}$

</div>

We can then solve for $K_{n}$:

<div align="center">

$K_{H_{2}O \rightarrow OH+H} = (\frac{204 atm}{1 atm})e^{- \frac{\Delta G_{H_{2}O}}{\hat{R}T_{c}}}$

$K_{H_{2}\rightarrow 2H} = (\frac{204 atm}{1 atm})e^{- \frac{\Delta G_{H_{2}}}{\hat{R}T_{c}}}$

</div>

Here we know that the fuel is $H_{2}$ therefore the mass of the fuel is 2.02 $\frac{g}{mol}$ and the mass of the oxidizer is $m_{oxidizer} =\frac{2.02}{\phi}$. And $l = \frac{m_{oxidizer}}{2 m_{O_2}}$. Once we have $l$ we re-visit the Chemistry Problem by solving the Atom Balance:

<div align="center">

$H_{2} +l O_{2} \rightarrow N_{H_{2}O}H_{2}O +N_{H_{2}}H_{2} + N_{H}H +N_{OH}OH$

$H_{2}: 2 = 2N_{H_{2}O} + 2N_{H_{2}} + N_{H} + N_{OH}$

$O_{2}: 2l = N_{H_{2}O} + N_{OH}$

</div>

Where $N_{H_{2}O}$, $N_{H_{2}}$, $N_{H}$, and $N_{OH}$ are the number of moles of each species. We can then solve for $N_{H_{2}O}$, $N_{H_{2}}$, $N_{H}$, and $N_{OH}$:

<div align="center">

$1 = N_{H_{2}O} + N_{H_{2}}$

$0.760592 = N_{H_{2}O}$

</div>

Resulting in us having 4 unknown values with only 2 equations. This a signal to us prompting the Law of Mass Action. When solving for the Moler Fraction of each species we get:

<div align="center">

$X_{H_{2}O} = \frac{N_{H_{2}O}}{N_{Total}} = \frac{N_{H_{2}O}}{N_{H_{2}O} + N_{H_{2}} + N_{H} + N_{OH}}$

$X_{H} = \frac{N_{H}}{N_{Total}} = \frac{N_{H}}{N_{H_{2}O} + N_{H_{2}} + N_{H} + N_{OH}}$

$X_{OH} = \frac{N_{OH}}{N_{Total}} = \frac{N_{OH}}{N_{H_{2}O} + N_{H_{2}} + N_{H} + N_{OH}}$

$K_{H_{2}O \rightarrow OH+H} = \frac{\frac{N_{OH}}{N_{H_{2}O} + N_{H_{2}}} \frac{N_{H}}{N_{H_{2}O} + N_{H_{2}}}}{\frac{N_{H_{2}O}}{N_{H_{2}O} + N_{H_{2}}}} = \frac{N_{H}N_{OH}}{(N_{H_{2}O} + N_{H_{2}})N_{H_{2}O}}$

</div>

Similarly for the $H_{2}\rightarrow 2H$ Reaction:

<div align="center">


| Species | $\nu_{(i)} '$ | $\nu_{(i)} "$ | $\nu_{(i)} " - \nu_{(i)} '$ |
|---------|-------|-------|---------------|
| $H_{2}$ | 1 | 0 | -1 |
| $H$ | 0 | 2 | 2 |

$X_{H_{2}} = \frac{N_{H_{2}}}{N_{Total}} = \frac{N_{H_{2}}}{N_{H_{2}O} + N_{H_{2}} + N_{H} + N_{OH}}$

$X_{H} = \frac{N_{H}}{N_{Total}} = \frac{N_{H}}{N_{H_{2}O} + N_{H_{2}} + N_{H} + N_{OH}}$

$K_{H_{2}\rightarrow 2H} = \frac{(\frac{N_{H}}{N_{H_{2}}+ N_{H_{2}O}})^{2}}{\frac{N_{H_{2}}}{N_{H_{2}}+ N_{H_{2}O}}} = \frac{N_{H}^{2}}{N_{H_{2}}(N_{H_{2}}+ N_{H_{2}O})}$

</div>

After this we know have to solve a system of coupled equations:

<div align="center">

$1 = N_{H_{2}O} + N_{H_{2}}$

$0.760592 = N_{H_{2}O}$

$K_{H_{2}O \rightarrow OH+H} = \frac{N_{H}N_{OH}}{(N_{H_{2}O} + N_{H_{2}})N_{H_{2}O}}$

$K_{H_{2}\rightarrow 2H} = \frac{N_{H}^{2}}{N_{H_{2}}(N_{H_{2}}+ N_{H_{2}O})}$ 

</div>

We can then solve for $N_{H_{2}O}$, $N_{H_{2}}$, $N_{H}$, and $N_{OH}$:

<div align="center">

| Species | $N_{(i)}$ | $X_{(i)}$ | $m_{(i)}$ | $Y_{(i)}$ | $\hat{c_{p}}_{(i)}$ | $c_{p_{(i)}}$ | $Y_{(i)}c_{p}$ |
|---------|-------|-------|---------------|-------|-------|---------------|-------|
| $H_{2}O$ | $N_{H_{2}O}$ | $X_{H_{2}O} = N_{H_{2}O}$ | $\hat{m} N_{H_{2}O}$ | $\frac{\hat{m}}{m_{(i)}}$ | $\hat{c_{p}} (T_{c})$ | $\frac{\hat{c_{p}}}{\hat{m}}$ | $Y_{H_{2}O}c_{p}$ |
| $H_{2}$ | $N_{H_{2}}$ | $X_{H_{2}} = N_{H_{2}}$ | $\hat{m} N_{H_{2}}$ | $\frac{\hat{m}}{m_{(i)}}$ | $\hat{c_{p}} (T_{c})$ | $\frac{\hat{c_{p}}}{\hat{m}}$ | $Y_{H_{2}}c_{p}$ |
| $H$ | $N_{H}$ | $X_{H} = N_{H}$ | $\hat{m} N_{H}$ | $\frac{\hat{m}}{m_{(i)}}$ | $\hat{c_{p}} (T_{c})$ | $\frac{\hat{c_{p}}}{\hat{m}}$ | $Y_{H}c_{p}$ |
| $OH$ | $N_{OH}$ | $X_{OH} = N_{OH}$ | $\hat{m} N_{OH}$ | $\frac{\hat{m}}{m_{(i)}}$ | $\hat{c_{p}} (T_{c})$ | $\frac{\hat{c_{p}}}{\hat{m}}$ | $Y_{OH}c_{p}$ |


</div>

Now to calculate $c_{p}$ we sum the $Y_{(i)}c_{p}$ for each species:

<div align="center">

$c_{p} = \sum_{i} Y_{(i)}c_{p}$

$m_{(total)} = \sum_{i} X_{(i)}m_{(i)}$

</div>

We perform this step to calculate the gas constant $R$:

<div align="center">

$R = \frac{\hat{R}}{m_{(total)}}$

$c_{v} = c_{p} - R$

</div>

Resulting in use getting a gamma of:

<div align="center">

$\gamma = \frac{c_{p}}{c_{p} - \frac{\hat{R}}{\sum_{i} X_{(i)}m_{(i)}}} = \frac{c_{p}}{c_{p} - \frac{\hat{R}}{m_{(total)}}}=\frac{c_{p}}{c_{p} - R}  =\frac{c_{p}}{c_{v}}$

</div>

## Oxygen Nozzle Performance

Begin by defining our Area Ratio: $A_{r} = \frac{\pi (R_{e})^{2}}{\pi (R_{t})^{2}}$. From here we are to solve for the Mach number:

<div align="center">

$A_{r} = (\frac{\gamma +1}{2})^{-\frac{\gamma+1}{2(\gamma-1)}} \frac{1}{M_{e}} (1+ \frac{\gamma-1}{2} (M_{e})^{2})^{\frac{\gamma+1}{2(\gamma -1)}}$

</div>

We can for Mach number $M_{e}$ by making use of `scipy.optimize.fsolve`, which returns the roots of a non-linear function. After solving for $M_{e}$ we can then solve for the pressure, Temperature, and velocity:


<div align="center">

$T_{e} = T_{c} (1+ \frac{\gamma-1}{2} (M_{e})^{2})^{-1}$

$p_{e} = p_{c} (1+ \frac{\gamma-1}{2} (M_{e})^{2})^{-\frac{\gamma}{\gamma-1}}$

$v_{e} = M_{e} \sqrt{\gamma R_{sp} T_{e}}$

</div>

After Solving for these values we can then calculate the mass flow rate:


<div align="center">

$\dot{m} = \frac{A_{t} p}{\sqrt{T_{c}}} \sqrt{\frac{\gamma}{R_{sp}}} (\frac{2}{\gamma + 1})^{\frac{\gamma + 1}{2(\gamma -1)}}$

</div>

Then onec you have the mass flow rate you can calculate the thrust, specific impulse, and thrust coefficient:

<div align="center">

$Thrust = \dot{m} v_{e}+p_{e}A_{e} - pA_{t}$

$I_{sp} = \frac{Thrust}{\dot{m} g_{0}}$

$C_{T} = \frac{Thrust}{pA_{t}}$

</div>
