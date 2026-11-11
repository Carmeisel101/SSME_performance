# Space Shuttle Main Engine Performance

## Introduction

The Space Shuttle Main Engine (SSME) is a liquid-fueled rocket engine used to power the Space Shuttle Orbiter. This report describes the methods and results of a performance model of the engine: the chamber temperature and equilibrium composition, and the ideal nozzle performance. It covers the original hydrogen/oxygen engine and the same engine run with fluorine as the oxidizer. For this project, we were given the parameters below:

| Combustion Chamber Data | |
|:--|:--|
| Fuel | Hydrogen |
| Oxidizer | $O_{2}$ (and $F_{2}$ for the fluorine variant) |
| Fuel/oxidizer mass ratio, $H_{2}/O_{2}$ | 0.166 |
| Hydrogen injection temperature | 850 K |
| Oxidizer injection temperature | 530 K |
| Chamber pressure | 204 atm |

| Nozzle Geometry | |
|:--|:--|
| Throat radius $R_{t}$ | 5.15 in |
| Exit radius $R_{e}$ | 45.35 in |
| $\theta_{p}$ | 32° (not used by the ideal-nozzle model) |

The fluorine engine keeps the same equivalence ratio, injection temperatures, chamber pressure and nozzle (see [Fluorine](#fluorine)).

## Oxygen

### First Enthalpy Exchanger

We begin by computing the equivalence ratio $\phi$, the fuel/oxidizer ratio divided by its stoichiometric value. For $2H_{2} + O_{2} \rightarrow 2H_{2}O$ the moles of $H_{2}$ per mole of $O_{2}$ are

$$\Phi = (F/O)\frac{\hat{m}_{O_{2}}}{\hat{m}_{H_{2}}} = 2.635, \qquad \phi = \frac{\Phi}{2} = 1.317$$

so the mixture is fuel-rich.

Next we calculate the enthalpy of each injected propellant, including a real-gas correction for the high chamber pressure:

$$\Delta \hat{h}_{H_{2}} = \hat{h}_{H_{2}}(850) + \hat{R}T^{2}\int^{pi}_{pf} \frac{1}{p} \frac{-\hat{b}_{H_{2}}p}{\hat{R}T^{2}} dp$$

$$\Delta \hat{h}_{O_{2}} = \hat{h}_{O_{2}}(530) + \hat{R}T^{2}\int^{pi}_{pf} \frac{1}{p} \frac{-\hat{b}_{O_{2}}p}{\hat{R}T^{2}} dp$$

| Species  | $\hat{b}$ |
|:--------:|:---------:|
| $H_{2}$  |  0.0266   |
| $O_{2}$  |  0.0318   |
| $H_{2}O$ |  0.03049  |
| $N_{2}$  |  0.0391   |
| $F_{2}$  |  0.02896  |
|   $HF$   |  0.0739   |

Where $\hat{b}$ is the Van der Waals co-volume (L/mol), listed above and in `exchanger.py`. $\Delta \hat{h}_{i}$ is calculated from both hydrogen and oxygen. Once we calculate $\Delta \hat{h}_{i}$ for each, we can calculate the enthalpy of the mixture.

$$\Delta \hat{h}_{mix} = \Phi \Delta \hat{h}_{H_{2}} + \Delta \hat{h}_{O_{2}}$$

With $\Phi = 2\phi$ moles of $H_{2}$ per mole of $O_{2}$, complete combustion (no dissociation) gives

$$O_{2} + \Phi H_{2} \rightarrow 2H_{2}O + (\Phi-2)H_{2}$$

The energy available to heat the products is the propellant enthalpy plus the heat released forming two moles of water:

$$\Delta \hat{h}_{total} = \Phi \Delta \hat{h}_{H_{2}} + \Delta \hat{h}_{O_{2}} - 2\Delta \hat{h}^{298}_{f,H_{2}O}$$

The first estimate of the chamber temperature, $T_{FEE}$, is the single temperature at which the products absorb exactly that energy. It is found with `scipy.optimize.brentq`:

$$2\left[\hat{h}_{H_{2}O}(T_{FEE}) - \hat{h}_{H_{2}O}(298)\right] + (\Phi-2)\left[\hat{h}_{H_{2}}(T_{FEE}) - \hat{h}_{H_{2}}(298)\right] = \Delta \hat{h}_{total}$$

This gives $T_{FEE} \approx 4617$ K. Because it ignores dissociation it is an upper bound; the chamber temperature comes from the second enthalpy exchanger in [Chamber Temperature](#chamber-temperature-second-enthalpy-exchanger).

### Combustion Chamber

For the Combustion Chamber, we have the two equations below:

$$H_{2}O \rightarrow OH + H$$

$$H_{2} \rightarrow 2H$$

The Gibbs Free Energy for each species is given by:

- $\hat{g}_{H_{2}O}$ = `h2o_table_g`($T_{C}$)

- $\hat{g}_{H_{2}}$ = `h2_table_g`($T_{C}$)

- $\hat{g}_{OH}$ = `oh_table_g`($T_{C}$)

- $\hat{g}_{H}$ = `h_table_g`($T_{C}$)

We can then calculate the Gibbs Free Energy of the reaction by:

$$\Delta \hat{g}_{H_{2}} = 2\hat{g}_{H} - \hat{g}_{H_{2}}$$

$$\Delta \hat{g}_{H_{2}O} = \hat{g}_{H} + \hat{g}_{OH} - \hat{g}_{H_{2}O}$$

Now is the point where we involve the **Law of Mass Action**:

| Species  | $\nu_{(i)} '$ | $\nu_{(i)} ''$ | $\nu_{(i)} '' - \nu_{(i)} '$ |
|:--------:|:-------------:|:--------------:|:----------------------------:|
| $H_{2}O$ |       1       |       0        |              -1              |
|   $OH$   |       0       |       1        |              1               |
|   $H$    |       0       |       1        |              1               |

We can now use:

$$K_{n} = \prod^{N_{s}}_{i=1} X_{(n)}^{(\nu_{(i)} '' - \nu_{(i)} ')} = \left(\frac{p}{p^{*}}\right)^{-\sigma_{v}}K_{p}$$

where $p^{*} = 1$ atm is the reference pressure of the tables.

Where;

$$\sigma_{v} = \sum^{N_{s}}_{i=1} (\nu_{(i)} '' - \nu_{(i)} ')$$

$$K_{p} = e^{- \frac{\Delta G}{\hat{R}T}}$$

$$\Delta G = \sum^{N_{s}}_{i=1} (\nu_{(i)} '' - \nu_{(i)} ') \hat{g}_{(i)}$$

Both reactions have $\sigma_{v} = +1$, so we can solve for the equilibrium constant $K_{n}$.

$$K_{H_{2}O} = \left(\frac{204 atm}{1 atm}\right)^{-1} e^{-\frac{\Delta \hat{g}_{H_{2}O}}{\hat{R}T_{C}}}$$

$$K_{H_{2}} = \left(\frac{204 atm}{1 atm}\right)^{-1} e^{-\frac{\Delta \hat{g}_{H_{2}}}{\hat{R}T_{C}}}$$

Now we move into the chemistry problem.

**Atom Balance**:

$$H_{2} + l O_{2} \rightarrow N_{H_{2}O} H_{2}O + N_{H_{2}} H_{2} + N_{H}H + N_{OH}OH$$

Per mole of $H_{2}$ fed, the moles of $O_{2}$ are

$$l = \frac{\hat{m}_{H_{2}}}{(F/O)\,\hat{m}_{O_{2}}} = \frac{1}{\Phi} = 0.3795$$

The atom balances for hydrogen and oxygen are:

$$H: \quad 2 = 2N_{H_{2}O} + 2N_{H_{2}} + N_{H} + N_{OH}$$

$$O: \quad 2l = N_{H_{2}O} + N_{OH}$$

That is four unknowns with only two equations, so the Law of Mass Action supplies the other two. The mole fractions are normalized by the full total,

$$X_{(i)} = \frac{N_{(i)}}{N_{Total}}, \qquad N_{Total} = N_{H_{2}O} + N_{H_{2}} + N_{H} + N_{OH},$$

so the two equilibrium relations become

$$K_{H_{2}O} = \frac{X_{H} X_{OH}}{X_{H_{2}O}} = \frac{N_{H} N_{OH}}{N_{H_{2}O} N_{Total}}$$

$$K_{H_{2}} = \frac{X_{H}^{2}}{X_{H_{2}}} = \frac{N_{H}^{2}}{N_{H_{2}} N_{Total}}$$

The four coupled equations (two atom balances and two equilibrium relations) are solved with `scipy.optimize.fsolve` in terms of $\ln N_{(i)}$, which keeps every mole number positive. The starting guess is complete combustion, $N_{H_{2}O} = 2l$ and $N_{H_{2}} = 1 - 2l$. At the chamber temperature $T_{c}$ ([Chamber Temperature](#chamber-temperature-second-enthalpy-exchanger)) the solution is:

$$N_{H_{2}O} = 0.6763, \quad N_{H_{2}} = 0.2559, \quad N_{H} = 0.0530, \quad N_{OH} = 0.0827$$

| Species | $N_{(i)}$ | $X_{(i)}$ | $m_{(i)}$ | $Y_{(i)}$ | $\hat{c}_{p(i)}$ | $c_{p(i)}$ | $Y_{(i)}c_{p}$ |
|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| $H_{2}O$ | $N_{H_{2}O}$ | $X_{H_{2}O} = \frac{N_{H_{2}O}}{N_{Total}}$ | $\hat{m}N_{H_{2}O}$ | $\frac{\hat{m}}{m_{(i)}}$ | $\hat{c}_{p(H_{2}O)}(T_{C})$ | $\frac{\hat{c}_{p}}{\hat{m}}$ | $Y_{H_{2}O}c_{p}$ |
| $H_{2}$ | $N_{H_{2}}$ | $X_{H_{2}} = \frac{N_{H_{2}}}{N_{Total}}$ | $\hat{m}N_{H_{2}}$ | $\frac{\hat{m}}{m_{(i)}}$ | $\hat{c}_{p(H_{2})}(T_{C})$ | $\frac{\hat{c}_{p}}{\hat{m}}$ | $Y_{H_{2}}c_{p}$ |
| $H$ | $N_{H}$ | $X_{H} = \frac{N_{H}}{N_{Total}}$ | $\hat{m}N_{H}$ | $\frac{\hat{m}}{m_{(i)}}$ | $\hat{c}_{p(H)}(T_{C})$ | $\frac{\hat{c}_{p}}{\hat{m}}$ | $Y_{H}c_{p}$ |
| $OH$ | $N_{OH}$ | $X_{OH} = \frac{N_{OH}}{N_{Total}}$ | $\hat{m}N_{OH}$ | $\frac{\hat{m}}{m_{(i)}}$ | $\hat{c}_{p(OH)}(T_{C})$ | $\frac{\hat{c}_{p}}{\hat{m}}$ | $Y_{OH}c_{p}$ |

Now we can calculate $c_{p}$, we sum the $Y_{(i)}c_{p}$ values for each species.

$$c_{p} = \sum_{i} Y_{(i)}c_{p}$$

$$c_{p} = Y_{H_{2}O}c_{p} + Y_{H_{2}}c_{p} + Y_{H}c_{p} + Y_{OH}c_{p}$$

Next we calculate the gas constant $R$.

$$R = \frac{\hat{R}}{m_{Total}}$$

where $m_{Total} = \sum_{i} X_{(i)} m_{(i)}$ is the mixture molar mass. The mole fractions are normalized by $N_{Total} = N_{H_{2}O} + N_{H_{2}} + N_{H} + N_{OH}$, which exceeds one once dissociation produces H and OH.

$$c_{v} = c_{p} - R$$

$$\gamma = \frac{c_{p}}{c_{p} - \frac{\hat{R}}{\sum_{i} X_{(i)}m_{(i)}}} = \frac{c_{p}}{c_{p}- \frac{\hat{R}}{m_{Total}}} = \frac{c_{p}}{c_{p}- R} = \frac{c_{p}}{c_{v}}$$

### Chamber Temperature (Second Enthalpy Exchanger)

Breaking $H_{2}O$ and $H_{2}$ apart into $H$ and $OH$ absorbs energy, so the dissociated products reach a lower temperature than $T_{FEE}$. The composition depends on temperature through $K_{n}$, and the temperature depends on the composition through the energy balance, so the two are solved together. For a trial temperature the combustion-chamber step gives the equilibrium mixture, and the energy balance per mole of $O_{2}$ is checked:

$$\sum_{i} N_{(i)} \left[\Delta \hat{h}^{298}_{f(i)} + \hat{h}_{(i)}(T_{c}) - \hat{h}_{(i)}(298)\right] = \Phi \Delta \hat{h}_{H_{2}} + \Delta \hat{h}_{O_{2}}$$

`brentq` adjusts the temperature between 2000 K and $T_{FEE}$ until both sides match. The sensible enthalpies of $H_{2}$, $O_{2}$ and $H_{2}O$ come straight from the tables, $H$ uses its constant $c_{p}$, and $OH$ is built from its tabulated enthalpy of formation:

$$\hat{h}_{OH}(T) - \hat{h}_{OH}(298) = \Delta \hat{h}_{f,OH}(T) - \Delta \hat{h}_{f,OH}(298) + \frac{1}{2}\left[\hat{h}_{H_{2}}(T) - \hat{h}_{H_{2}}(298)\right] + \frac{1}{2}\left[\hat{h}_{O_{2}}(T) - \hat{h}_{O_{2}}(298)\right]$$

The same construction applied to $H_{2}O$ reproduces its tabulated sensible enthalpy within 0.07%. The result is $T_{c} \approx 3949$ K, which is used for $\gamma$, $R$ and the nozzle. A Cantera equilibrium calculation with the same inputs gives 3881 K; most of the difference is that this model leaves out $O$ and $O_{2}$.

### Nozzle Performance

We begin by defining the Area Ratio:

$$A_{r} = \frac{\pi R_{e}^{2}}{\pi R_{t}^{2}}$$

Allowing for us to solve for the Mach number:

$$A_{r} = \left(\frac{\gamma + 1}{2}\right)^{-\frac{\gamma + 1}{2(\gamma - 1)}} \frac{1}{M_{e}} \left[1 + \frac{\gamma - 1}{2} M_{e}^{2}\right]^{\frac{\gamma + 1}{2(\gamma - 1)}}$$

We can solve for Mach number $M_{e}$ by making use of the `scipy.optimize.fsolve` function, which returns the roots of a non-linear function. After solving for $M_{e}$ we can calculate the exit velocity $v_{e}$, exit Pressure $p_{e}$ and exit Temperature $T_{e}$.

$$T_{e} = T_{c} \left[1 + \frac{\gamma - 1}{2} M_{e}^{2}\right]^{-1}$$

$$p_{e} = p_{c} \left[1 + \frac{\gamma - 1}{2} M_{e}^{2}\right]^{-\frac{\gamma}{\gamma - 1}}$$

$$v_{e} = M_{e} \sqrt{\gamma R T_{e}}$$

With the exit velocity, pressure and temperature known, we can calculate the mass flow rate $\dot{m}$.

$$\dot{m} = \frac{A_{t} p_{c} \sqrt{\gamma}}{\sqrt{R T_{c}}} \left[ \frac{2}{\gamma + 1} \right]^{\frac{\gamma + 1}{2(\gamma - 1)}}$$

Then once you have the mass flow rate you can calculate the thrust $\mathcal{T}$, specific impulse $I_{sp}$ and the thrust coefficient $C_{T}$.

$$\mathcal{T} = \dot{m} v_{e} + (p_{e} - p_{a}) A_{e}$$

where $p_{a}$ is the ambient pressure ($p_{a} = 0$ in vacuum, $p_{a} = 101325$ Pa at sea level).

$$I_{sp} = \frac{\mathcal{T}}{\dot{m} g_{0}}$$

$$C_{T} = \frac{\mathcal{T}}{p_{c} A_{t}}$$

### Oxygen Results

| Output | Vacuum | Sea level |
|:--|:--:|:--:|
| Thrust | 2.151 MN | 1.729 MN |
| $C_{T}$ | 1.936 | 1.556 |
| $I_{sp}$ | 478.8 s | 384.8 s |

The mass flow is 457.9 kg/s, the exit Mach number 4.70 and the exit pressure 18.7 kPa. For reference, the real RS-25 at 109% power is rated at about 2.28 MN and 452 s in vacuum. The model sits above it mainly because of the warm injection temperatures in the givens and the loss-free nozzle.

## Fluorine

The same engine is run with fluorine in place of oxygen, at the same equivalence ratio, injection temperatures ($H_{2}$ at 850 K, $F_{2}$ at 530 K), chamber pressure and nozzle. The stoichiometric reaction is $H_{2} + F_{2} \rightarrow 2HF$, so $\phi = 1.317$ moles of $H_{2}$ per mole of $F_{2}$, an oxidizer/fuel mass ratio of 14.3.

### First Enthalpy Exchanger

Complete combustion, per mole of $F_{2}$, is

$$F_{2} + \phi H_{2} \rightarrow 2HF + (\phi - 1)H_{2}$$

and $T_{FEE}$ follows from

$$2\left[\hat{h}_{HF}(T_{FEE}) - \hat{h}_{HF}(298)\right] + (\phi-1)\left[\hat{h}_{H_{2}}(T_{FEE}) - \hat{h}_{H_{2}}(298)\right] = \phi\Delta\hat{h}_{H_{2}} + \Delta\hat{h}_{F_{2}} - 2\Delta\hat{h}^{298}_{f,HF}$$

which gives $T_{FEE} \approx 7107$ K.

### Combustion Chamber

Both $H_{2}$ and $HF$ dissociate at these temperatures, so the products are $HF$, $H_{2}$, $H$ and $F$:

$$F: \quad 2 = N_{HF} + N_{F}, \qquad H: \quad 2\phi = N_{HF} + 2N_{H_{2}} + N_{H}$$

$$K_{H_{2} \rightarrow 2H} = \frac{N_{H}^{2}}{N_{H_{2}} N_{Total}}, \qquad K_{HF \rightarrow H + F} = \frac{N_{H} N_{F}}{N_{HF} N_{Total}}$$

with $K_{n} = (p/p^{*})^{-1} e^{-\Delta \hat{g}/\hat{R}T}$ for both reactions. $F_{2}$ is left out of the products; at the chamber conditions $X_{F_{2}} \approx 10^{-5}$. The system is solved in $\ln N_{(i)}$ as for oxygen, and the chamber temperature again comes from the energy balance on the equilibrium products.

### Results

|                     | $H_{2}/O_{2}$ | $H_{2}/F_{2}$ |
|:--------------------|:-------------:|:-------------:|
| $T_{FEE}$           |    4617 K     |    7107 K     |
| $T_{c}$             |    3949 K     |    5190 K     |
| Molar mass          |  13.26 g/mol  |  15.43 g/mol  |
| $\gamma$ (frozen)   |     1.199     |     1.319     |
| $\dot{m}$           |  457.9 kg/s   |  445.7 kg/s   |
| $I_{sp}$, vacuum    |    478.8 s    |    459.5 s    |
| $I_{sp}$, sea level |    384.8 s    |    362.8 s    |

The equilibrium mole fractions are $X_{HF} = 0.678$, $X_{H_{2}} = 0.082$, $X_{H} = 0.159$ and $X_{F} = 0.082$. About a quarter of the gas is atoms. Because the nozzle model keeps the chamber composition frozen, the energy stored in those atoms is never released by recombination, and the fluorine $I_{sp}$ comes out below the oxygen one even though the chamber is 1240 K hotter. A shifting-equilibrium nozzle would recover much of that energy, so this result is a lower bound.

### Data

The course tables have no data for atomic fluorine and no Gibbs energy for $HF$, and the $H$ table stopped at 5000 K. These were taken from the [NIST-JANAF thermochemical tables](https://janaf.nist.gov/): the $F$ atom table, the $\Delta \hat{g}_{f}$ column for $HF$, and the 5100–6000 K rows for $H$. A calculation that uses JANAF for every species gives $T_{c} = 5197.5$ K and the same mole fractions to within 0.001. Six $H^{*}-H^{*}_{298}$ entries in the course tables were typos that broke the slope set by $c_{p}$ and were corrected: $H_{2}$ at 1900 and 2400 K, $F_{2}$ at 1700 and 4500 K, and $HF$ at 3400 and 3700 K.
