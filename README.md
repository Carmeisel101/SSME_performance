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

## Running it

Requires Python 3.11+ (tested on 3.15).

```bash
make setup   # create .venv and install requirements.txt
make         # run the model (SSMEMain.py)
make test    # regression checks
```


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

With $\Phi = 2\phi$ moles of $H_{2}$ per mole of $O_{2}$, complete combustion (no dissociation) gives:

<div align="center">

$O_{2} + \Phi H_{2} \rightarrow 2H_{2}O + (\Phi-2)H_{2}$

</div>

The energy available to heat the products is the propellant enthalpy plus the heat released forming two moles of water:

<div align="center">

$\Delta \hat{h}_{total} = \Phi \Delta \hat{h}_{H_{2}} + \Delta \hat{h}_{O_{2}} - 2\Delta \hat{h}_{f,H_{2}O}^{298}$

</div>

The first estimate of the chamber temperature, $T_{FEE}$, is the single temperature at which the products absorb exactly that energy, found with `scipy.optimize.brentq`:

<div align="center">

$2\left[\hat{h}_{H_{2}O}(T_{FEE}) - \hat{h}_{H_{2}O}(298)\right] + (\Phi-2)\left[\hat{h}_{H_{2}}(T_{FEE}) - \hat{h}_{H_{2}}(298)\right] = \Delta \hat{h}_{total}$

</div>

This gives $T_{FEE} \approx 4617$ K. Because it ignores dissociation it is an upper bound; the chamber temperature itself comes from the second enthalpy exchanger below.

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

Where $\bar{g_{fh}}$ is the Gibbs Free Energy of formation of $H_{2}$, $\bar{g_{foh}}$ is the Gibbs Free Energy of formation of OH, and $\bar{g_{fh_{2}o}}$ is the Gibbs Free Energy of formation of $H_{2}O$. These values are read from the `Iterater.py` file.

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

Where $N_{H_{2}O}$, $N_{H_{2}}$, $N_{H}$, and $N_{OH}$ are the number of moles of each species. That is 4 unknowns with only 2 atom-balance equations, so the Law of Mass Action supplies the other two. The mole fractions use the full total:

<div align="center">

$X_{(i)} = \frac{N_{(i)}}{N_{Total}}, \quad N_{Total} = N_{H_{2}O} + N_{H_{2}} + N_{H} + N_{OH}$

$K_{H_{2}O \rightarrow OH+H} = \frac{X_{H} X_{OH}}{X_{H_{2}O}} = \frac{N_{H}N_{OH}}{N_{H_{2}O} N_{Total}}$

</div>

Similarly for the $H_{2}\rightarrow 2H$ Reaction:

<div align="center">


| Species | $\nu_{(i)} '$ | $\nu_{(i)} "$ | $\nu_{(i)} " - \nu_{(i)} '$ |
|---------|-------|-------|---------------|
| $H_{2}$ | 1 | 0 | -1 |
| $H$ | 0 | 2 | 2 |

$K_{H_{2}\rightarrow 2H} = \frac{X_{H}^{2}}{X_{H_{2}}} = \frac{N_{H}^{2}}{N_{H_{2}} N_{Total}}$

</div>

This gives 4 coupled equations in 4 unknowns:

<div align="center">

$2 = 2N_{H_{2}O} + 2N_{H_{2}} + N_{H} + N_{OH}$

$2l = N_{H_{2}O} + N_{OH}$

$K_{H_{2}O \rightarrow OH+H} = \frac{N_{H}N_{OH}}{N_{H_{2}O} N_{Total}}$

$K_{H_{2}\rightarrow 2H} = \frac{N_{H}^{2}}{N_{H_{2}} N_{Total}}$

</div>

The code solves them with `scipy.optimize.fsolve` in terms of $\ln N_{(i)}$, which keeps every mole number positive, starting from complete combustion ($N_{H_{2}O} = 2l$, $N_{H_{2}} = 1 - 2l$). At the chamber temperature $T_{c}$ (see below) this gives $N_{H_{2}O} = 0.6763$, $N_{H_{2}} = 0.2559$, $N_{H} = 0.0530$, $N_{OH} = 0.0827$.

We can then solve for $N_{H_{2}O}$, $N_{H_{2}}$, $N_{H}$, and $N_{OH}$:

<div align="center">

| Species | $N_{(i)}$ | $X_{(i)}$ | $m_{(i)}$ | $Y_{(i)}$ | $\hat{c_{p}}_{(i)}$ | $c_{p_{(i)}}$ | $Y_{(i)}c_{p}$ |
|---------|-------|-------|---------------|-------|-------|---------------|-------|
| $H_{2}O$ | $N_{H_{2}O}$ | $X_{H_{2}O} = \frac{N_{H_{2}O}}{N_{Total}}$ | $\hat{m} N_{H_{2}O}$ | $\frac{\hat{m}}{m_{(i)}}$ | $\hat{c_{p}} (T_{c})$ | $\frac{\hat{c_{p}}}{\hat{m}}$ | $Y_{H_{2}O}c_{p}$ |
| $H_{2}$ | $N_{H_{2}}$ | $X_{H_{2}} = \frac{N_{H_{2}}}{N_{Total}}$ | $\hat{m} N_{H_{2}}$ | $\frac{\hat{m}}{m_{(i)}}$ | $\hat{c_{p}} (T_{c})$ | $\frac{\hat{c_{p}}}{\hat{m}}$ | $Y_{H_{2}}c_{p}$ |
| $H$ | $N_{H}$ | $X_{H} = \frac{N_{H}}{N_{Total}}$ | $\hat{m} N_{H}$ | $\frac{\hat{m}}{m_{(i)}}$ | $\hat{c_{p}} (T_{c})$ | $\frac{\hat{c_{p}}}{\hat{m}}$ | $Y_{H}c_{p}$ |
| $OH$ | $N_{OH}$ | $X_{OH} = \frac{N_{OH}}{N_{Total}}$ | $\hat{m} N_{OH}$ | $\frac{\hat{m}}{m_{(i)}}$ | $\hat{c_{p}} (T_{c})$ | $\frac{\hat{c_{p}}}{\hat{m}}$ | $Y_{OH}c_{p}$ |


</div>

Now to calculate $c_{p}$ we sum the $Y_{(i)}c_{p}$ for each species:

<div align="center">

$c_{p} = \sum_{i} Y_{(i)}c_{p}$

$m_{(total)} = \sum_{i} X_{(i)}m_{(i)}$

$X_{(i)} = \frac{N_{(i)}}{N_{Total}}, \quad N_{Total} = N_{H_{2}O} + N_{H_{2}} + N_{H} + N_{OH}$

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

## Oxygen - Chamber Temperature (Second Enthalpy Exchanger)

Breaking $H_{2}O$ and $H_{2}$ apart into $H$ and $OH$ absorbs energy, so the dissociated products reach a lower temperature than $T_{FEE}$. The composition depends on temperature (through $K_{n}$) and the temperature depends on the composition (through the energy balance), so `ChamberTemperature` solves both together: for a trial $T$ it computes the equilibrium mixture with the combustion-chamber step above and checks the energy balance, per mole of $O_{2}$,

<div align="center">

$\sum_{i} N_{(i)} \left[\Delta \hat{h}_{f(i)}^{298} + \hat{h}_{(i)}(T_{c}) - \hat{h}_{(i)}(298)\right] = \Phi \Delta \hat{h}_{H_{2}} + \Delta \hat{h}_{O_{2}}$

</div>

and `brentq` adjusts $T$ between 2000 K and $T_{FEE}$ until both sides match. The sensible enthalpies come from the tables: $H_{2}$, $O_{2}$ and $H_{2}O$ directly, $H$ from its constant $c_{p}$, and $OH$ from its tabulated enthalpy of formation, $\hat{h}_{OH} = \Delta \hat{h}_{f,OH}(T) + \frac{1}{2}\hat{h}_{H_{2}} + \frac{1}{2}\hat{h}_{O_{2}}$. The result is $T_{c} \approx 3949$ K, which is used for $\gamma$, $R$ and the nozzle. A Cantera equilibrium calculation with the same inputs gives 3881 K; most of the difference is that this model leaves out $O$ and $O_{2}$.

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

Then once you have the mass flow rate you can calculate the thrust, specific impulse, and thrust coefficient:

<div align="center">

$Thrust = \dot{m} v_{e}+(p_{e} - p_{a})A_{e}$

$I_{sp} = \frac{Thrust}{\dot{m} g_{0}}$

$C_{T} = \frac{Thrust}{pA_{t}}$

</div>

The thrust comes from a momentum balance on a control volume around the engine, so the only pressure term is the exit-plane imbalance $(p_{e} - p_{a})A_{e}$. The code reports vacuum ($p_{a} = 0$) and sea-level ($p_{a} = 101325$ Pa) values. Chamber pressure appears only in the normalization of $C_{T}$.

| Output | Vacuum | Sea level |
|--------|--------|-----------|
| Thrust | 2.151 MN | 1.729 MN |
| $C_{T}$ | 1.936 | 1.556 |
| $I_{sp}$ | 478.8 s | 384.8 s |

For reference, the real RS‑25 at 109% power is rated at about 2.28 MN and 452 s in vacuum.


## Code Flow


```mermaid
    flowchart
    A[Givens] --> B[Parameters]
    B --> SSME
    D -->|def| FCC
    FCC1 --> |req| Request
    Request --> |return| del_h_hat
    phi --> FCC2
    T_alg --> |req| Request
    Request --> |return| T_alg
    E -->|def| Comb
    Gibbs --> |req| Request
    Request --> |return| Gibbs
    cpval --> |req| Request
    Request --> |return| cpval
    N --> |def| Nozzle


    subgraph SSME[SSMEMain.py]
    C[Read in Givens] --> phi["equivalence ratio"]
    phi --> D[FCC]
    D--> E[Comb Chamber]
    E --> N[Nozzle]

  
    end
    

    subgraph FCC[FCC]
    FCC1["h_(i)(InjectionTemp)"] 
    b_hat["b_hat_table"]
    del_h_hat["del_h_hat_(i)"]
    b_hat --> del_h_hat
    FCC2["mols_(i)"]
    h_total["del_h_total"]
    T_alg["T_algorithm"]

    del_h_hat --> h_total
    FCC2 --> h_total
    h_total --> T_alg
    end 

    subgraph Iter["Iterater.py"]
    Request["Request_(i) (input,output)"] 

    Request --> |"(T,h)"|Iter1
    Request --> |"(h,T)"|Iter2
    Request --> |"(T,g)"|Iter3
    Request --> |"(T,cp)"|Iter4
    Request --> |"(T,S)"|Iter5
    Iter1["(i)_tableH = (i)_table_f"]
    Iter2["(i)_tableT"]
    Iter3["(i)_tableg"]
    Iter4["(i)_table_cp"]
    Iter5["(i)_tableS"]

    MachS["MachSolve"]
    end

    subgraph Comb["Comb Chamber"]
    Atom["Atom Balance"]
    Gibbs["Gibbs Free Energy"]
    gibbs_r["g_reaction"]
    LMA["LawOfMassAction_df"]
    cpval["cp_values_(i)"]
    Gibbs --> gibbs_r
    gibbs_r --> k_n["K_(n)"]
    T_alg --> k_n

    Atom --> eq["equations"]
    eq --> sol["solution(s)"]
    k_n --> eq
    sol --> LMA
    cpval --> LMA
    LMA --> R["R"]
    R --> cv["cv"]
    cv --> gamma["gamma"]
    end

    subgraph Nozzle["Nozzle"]
    A_r["Area Ratio"]
    gam["gamma"]
    gam --> Mach
    A_r --> Mach
    Mach -->|req| MachS
    MachS --> |return|Mach
    Mach --> mdot["mdot"]
    Mach --> T_e["T_e"]
    Mach --> p_e["p_e"]
    Mach --> v_e["v_e"]

    v_e --> Thrust["Thrust"]
    mdot --> Thrust
    p_e --> Thrust
    c_T["c_T"]
    Thrust --> c_T
    Thrust --> Isp["Isp"]

    end

    subgraph Mach["mach"]
    end
```