# WACNI-rice

## Introduction
Carbon and nitrogen metabolism are the basis for plant growth and development. To predict grain yield formation directly from the molecular processes in different source, sink and transport organs, we develop a mechanistic model of Whole plAnt Carbon Nitrogen Interaction (WACNI). 

In contrast to previous methods that have only considered phloem sucrose transport from leaf to grain, or relied on preset sink growth pattern, the key innovation of WACNI lies in its ability to integrate all essential biochemical and biophysical processes into a single model. This integration enables the kinetic simulation of authentic biological dynamics within a rice plant, spanning from flowering to harvest, and operating at a second-by-second time scale. 

One significant advantage of this model is that it eliminates the need for predefined growth patterns, allowing for a more accurate depiction of the underlying biological processes. Furthermore, this feature opens up possibilities for rational optimization of the fundamental biochemical and biophysical processes for higher grain yield. 

<p align="center">
  <img src="./images/GraphicAbstract.pdf" width = "1000" alt="method" align=center />
</p>


WACNI models both carbon and nitrogen metabolism, and kinetically simulates rates of major basic biochemical and physiological processes in a plant by ordinary differentiation equations. Fourteen types of biochemical and biophysical processes involved in different source, sink and transport organs are mathematically represented. 

These processes include assimilation, transport and utilization of six representative primary metabolites, i.e., triose phosphates (TP), sucrose (Suc), starch, inorganic nitrogen (I-N, including NH4+ and NO3-), free form of organic nitrogen (O-N, including amino acids and amides), and proteins. The 14 groups of simulated processes are listed below: 

> 1, root nitrogen uptake; 
> 
> 2, root nitrogen assimilation; 
> 
> 3, root growth; 
> 
> 4, root senescence; 
> 
> 5, leaf CO2 and nitrogen assimilation; 
> 
> 6, leaf triose phosphate, sucrose and starch interconversion; 
> 
> 7, leaf organic nitrogen and protein interconversion; 
> 
> 8, leaf senescence; 
> 
> 9, grain growth; 
> 
> 10, grain starch and protein synthesis; 
> 
> 11, stem sucrose and starch, organic nitrogen and protein interconversion; 
> 
> 12, phloem transport; 
> 
> 13, transmembrane transport; 
> 
> 14, symplastic diffusion between phloem and stem. 

WACNI further incorporates the interaction between these metabolites and plant developmental processes, e.g. root growth, grain volume expansion (endosperm cell division), grain filling (starch and protein synthesis in endosperm), root senescence and leaf senescence. 

It is worth emphasizing that WACNI simulates a plant in a canopy with user-defined row and column distances to neighboring plants (where plants share incident radiation and soil nitrogen), rather than a single-grown plant. Particularly, a sun-shade model was employed to simulate canopy photosynthesis. Additionally, a custom-built 1-dimensional root metabolism model was developed to simulate various processes, including root nitrogen uptake, nitrogen assimilation, root growth, and root senescence, based on previous experimental and theoretical studies. This enables us to easily scale single plant simulation results to a whole canopy by multiplying a factor of plant density per area.

<p align="center">
  <img src="./images/FourteenProcesses.pdf" width = "1000" alt="method" align=center />
</p>

WACNI follows a principle of modular design. Specifically, WACNI comprises six modules, i.e., root, leaf, grain, stem (including culm and sheath), transport (vascular transport system, including xylem and phloem), and respiration. Metabolites exchange between modules by trans-membrane transport. Each module comprises a few previously well-established sub-models, e.g., the sun-shade canopy photosynthesis model, the stomatal conductance model, the phloem transport model, and organ respiration model. New sub-models are developed where no established sub-models are available, e.g., root metabolism model, grain development and storage model, and stem metabolite homeostasis model.

<p align="center">
  <img src="./images/SixModules.pdf" width = "1000" alt="method" align=center />
</p>

The kinetic parameters for each process, initial metabolite concentrations, carbon and nitrogen mass in different organs at the flowering stage are input variables for the model. Finally, all types of plant-level physiological dynamics from flowering to harvest just emerge as output from this bottom-up simulation.


## System Requirements
### Hardware requirements
`WACNI-rice` package requires only a standard computer with enough RAM to support the in-memory operations. 

### Software requirements
#### OS Requirements
The developmental version of the package has been tested on the following systems:
+ Windows
+ Linux 

#### MATLAB
WACNI is developed under MATLAB (the MathWorks Inc., Natick, MA, USA) (version >= 2012b).


## Usage
WACNI can be run from the command window of MATLAB.

### Files
The `WACNI-rice` package contains the following essential files for running the WACNI model.

> RunModel.m
> 
> simulation_main.m
> 
> MyEventFunction.m
> 
> GA_for_evolutionaryPopulation.m
> 
> weather_input.txt


###  1. Basic simulation (1)
To run basic simulation for once, without showing plant carbon and nitrogen economy during grain filling, type in the following commands in MATLAB command window:

```
c=num2cell(ones(1,49));
[Tt,simulation]=RunModel(c,0,0,0,0);
```

###  2. Basic simulation (2)
To run basic simulation for once, with showing plant carbon and nitrogen economy during grain filling, using the following commands:

```
c=num2cell(ones(1,49));
[Tt,simulation]=RunModel(c,0,0,0,1);
```

###  3. Basic simulation (3)
To change one of the parameters and run simulation for once (e.g. change c{7}: coef_K_protein_oN_leaf_cat to 0.5), using the following commands:

```
c=num2cell(ones(1,49));
c{7}=0.5;
[Tt,simulation]=RunModel(c,0,0,0,0);
```

Note: the mapping between index and parameters are listed below:


|    Abbr    | Description |
| -------- | ----------- |
|c{1}|Prefix of output txt file|
|c{2}|Prefix of output plot file|
|c{3}|coef\_Vc_J|
|c{4}|coef\_PPFD|
|c{5}|coef\_Area_leaf|
|c{6}|coef\_Phi_CO2|
|c{7}|coef\_K\_protein\_oN\_leaf\_cat|
|c{8}|coef\_lower\_car\_photo|
|c{9}|coef\_K\_oN\_leaf\_load\_cat|
|c{10}|coef\_K\_suc\_leaf\_load\_cat|
|c{11}|coef\_K\_N\_assimi\_leaf\_cat|
|c{12}|coef\_K\_TP\_suc\_leaf\_cat|
|c{13}|coef\_K\_TP\_starch\_leaf\_cat|
|c{14}|coef\_K\_starch\_suc\_leaf\_cat|
|c{15}|coef\_K\_oN\_protein\_leaf\_cat|
|c{16}|coef\_conc\_lower\_oN\_to\_leaf\_protein|
|c{17}|coef\_K\_iN\_leaf\_unload\_cat|
|c{18}|coef\_gamma_star| 
|c{19}|coef\_grain_number|
|c{20}|coef\_K\_grain\_growth\_cat|
|c{21}|coef\_K\_grain\_starch\_store\_cat|
|c{22}|coef\_K\_suc\_grain\_unload\_cat|
|c{23}|coef\_K\_grain\_protein\_store\_cat|
|c{24}|coef\_K\_oN\_grain\_unload\_cat|
|c{25}|coef\_grain\_width\_max|
|c{26}|coef\_grain\_length\_width\_ratio|
|c{27}|coef\_Starch\_stem|
|c{28}|coef\_K\_suc\_starch\_stem\_cat|
|c{29}|coef\_K\_starch\_suc\_stem\_cat|
|c{30}|coef\_K\_oN\_protein\_stem\_cat|
|c{31}|coef\_K\_protein\_oN\_stem\_cat|
|c{32}|coef\_Fresh\_root|
|c{33}|coef\_conc\_iN\_soil|    
|c{34}|coef\_K\_iN\_absorb\_cat|
|c{35}|coef\_K\_iN\_root\_load\_cat|
|c{36}|coef\_K\_oN\_root\_load\_cat|
|c{37}|coef\_K\_suc\_root\_unload\_cat|
|c{38}|coef\_K\_oN\_root\_unload\_cat|  
|c{39}|coef\_K\_N\_root\_assimi\_cat|
|c{40}|coef\_K\_root\_growth\_cat|
|c{41}|coef\_R\_l\_s\_phloem0|
|c{42}|coef\_R\_r\_l\_phloem0|
|c{43}|coef\_K\_oN\_xylem\_to\_lphloem\_cat|
|c{44}|coef\_tillerNum|
|c{45}|coef\_CO2a|
|c{46}|protein\_content\_leaf|
|c{47}|grain\_struct\_N|
|c{48}|stem\_struct\_mass|
|c{49}|protein\_stem|


###  4. Simulating different nitrogen regimes

To simulate different nitrogen treatments shown in Zhao et al. (2015), using the following commands:

```
c=num2cell(ones(1,49));

c{1}='MN';c{2}=c{1};c{33}=1;c{46}=1;c{49}=1;c{19}=1;c{44}=1;c{44}=11/12; RunModel(c,0,0,0,0);

c{1}='0N';c{2}=c{1};c{33}=0.2;c{46}=0.66;c{49}=0.66;c{19}=0.43/0.6;c{44}=0.6*11/12; RunModel(c,0,0,0,0);

c{1}='LN';c{2}=c{1};c{33}=0.5;c{46}=0.84;c{49}=0.84;c{19}=0.75/0.9;c{44}=0.9*11/12; RunModel(c,0,0,0,0);

c{1}='HN';c{2}=c{1};c{33}=1.5;c{46}=0.98;c{49}=0.98;c{19}=1.18/1.23;c{44}=1.23*11/12; RunModel(c,0,0,0,0);

```
The meaning of the abbreviations used in the code is listed as follows:

|    Abbr    | Description |
| -------- | ----------- |
|  MN  | Medium nitrogen |
|  0N  | No nitrogen |
|  LN   | Low nitrogen |
|  HN   | High nitrogen |


###  5. Simulating different light regimes
To simulate different lighting regimes shown in Kobata et al. (2000), using the following commands:

```
c=num2cell(ones(1,49));

c{1}='CK';c{2}=c{1}; RunModel(c,0,0,0,0);

c{1}='S10_25';c{2}=c{1}; RunModel(c,25,1,0,0);

c{1}='S10_50';c{2}=c{1}; RunModel(c,50,1,0,0);

c{1}='S10_75';c{2}=c{1}; RunModel(c,75,1,0,0);

c{1}='S10_25_Sp50';c{2}=c{1}; RunModel(c,25,1,50,0);

c{1}='S10_50_Sp50';c{2}=c{1}; RunModel(c,50,1,50,0);

c{1}='S10_75_Sp50';c{2}=c{1}; RunModel(c,75,1,50,0);
```

The meaning of the abbreviations used in the code is listed as follows:

|    Abbr    | Description |
| -------- | ----------- |
|  CK  | Normal condition |
|  S10_25  | Shading 25% at the first 10 days|
|  S10_50   | Shading 50% at the first 10 days|
|  S10_75   | Shading 75% at the first 10 days|
|  S10_25_Sp50  | Shading 25% followed by 50% spacing|
|  S10_50_Sp50   | Shading 50% followed by 50% spacing|
|  S10_75_Sp50   | Shading 75% followed by 50% spacing|


###  6. Simulating different soil-nitrogen & air-CO2 concentration combinations
To simulate different soil nitrogen & air CO2 concentrations shown in Kim et al. (2001), using the following commands:

```
c=num2cell(ones(1,49));c{1}='MN';c{2}=c{1};c{19}=145/160;c{45}=390/400; RunModel(c,0,0,0,0);

c=num2cell(ones(1,49));c{1}='LN';c{2}=c{1};c{19}=145/160;c{25}=1.07^(1/3);c{44}=0.87;c{5}=0.87/c{44};c{32}=1/c{44};c{48}=0.87/c{44};c{27}=0.87/c{44};c{4
9}=0.95;c{19}=c{19}*0.96;c{46}=0.95;c{33}=0.5;c{45}=390/400; RunModel(c,0,0,0,0);

c=num2cell(ones(1,49));c{1}='HN';c{2}=c{1};c{19}=145/160;c{25}=0.98^(1/3);c{44}=1.07;c{5}=1.06/c{44};c{32}=1.03/c{44};c{48}=1.06/c{44};c{27}=1.06/c{44};
c{49}=1.17;c{19}=c{19}*1.05;c{46}=1.18;c{33}=1.5;c{45}=390/400; RunModel(c,0,0,0,0);

c=num2cell(ones(1,49));c{1}='LN-
FACE';c{2}=c{1};c{19}=145/160;c{25}=1.09^(1/3);c{22}=1;c{44}=0.89;c{5}=1.21/c{44};c{32}=1.18/c{44};c{48}=1.21/c{44};c{27}=1.21/c{44};c{49}=0.96;c{19}=c{
19}*0.95;c{46}=0.70;c{33}=0.5;c{45}=690/400; RunModel(c,0,0,0,0);

c=num2cell(ones(1,49));c{1}='MN-
FACE';c{2}=c{1};c{19}=145/160;c{25}=1.04^(1/3);c{22}=1;c{44}=1.09;c{5}=1.45/c{44};c{32}=1.46/c{44};c{48}=1.45/c{44};c{27}=1.45/c{44};c{49}=1.06;c{19}=c{
19}*1.01;c{46}=0.79;c{33}=1;c{45}=690/400; RunModel(c,0,0,0,0);

c=num2cell(ones(1,49));c{1}='HN-
FACE';c{2}=c{1};c{19}=145/160;c{25}=0.93^(1/3);c{22}=1;c{44}=1.15;c{5}=1.51/c{44};c{32}=1.35/c{44};c{48}=1.51/c{44};c{27}=1.51/c{44};c{49}=1.33;c{19}=c{
19}*1.13;c{46}=1.01;c{33}=1.5;c{45}=690/400; RunModel(c,0,0,0,0);

```


To simulate `vm_Nupt=64%` under FACE condition, using the following commands: 

```
c=num2cell(ones(1,49));c{1}='LN-FACE-
vm_Nup64';c{2}=c{1};c{19}=145/160;c{25}=1.09^(1/3);c{22}=1;c{34}=0.64;c{44}=0.89;c{5}=1.21/c{44};c{32}=1.18/c{44};c{48}=1.21/c{44};c{27}=1.21/c{44};c{49
}=0.96;c{19}=c{19}*0.95;c{46}=0.70;c{33}=0.5;c{45}=690/400; RunModel(c,0,0,0,0);

c=num2cell(ones(1,49));c{1}='MN-FACE-
vm_Nup64';c{2}=c{1};c{19}=145/160;c{25}=1.04^(1/3);c{22}=1;c{34}=0.64;c{44}=1.09;c{5}=1.45/c{44};c{32}=1.46/c{44};c{48}=1.45/c{44};c{27}=1.45/c{44};c{49
}=1.06;c{19}=c{19}*1.01;c{46}=0.79;c{33}=1;c{45}=690/400; RunModel(c,0,0,0,0);

c=num2cell(ones(1,49));c{1}='HN-FACE-
vm_Nup64';c{2}=c{1};c{19}=145/160;c{25}=0.93^(1/3);c{22}=1;c{34}=0.64;c{44}=1.15;c{5}=1.51/c{44};c{32}=1.35/c{44};c{48}=1.51/c{44};c{27}=1.51/c{44};c{49}=1.33;c{19}=c{19}*1.13;c{46}=1.01;c{33}=1.5;c{45}=690/400; RunModel(c,0,0,0,0);
```


The meaning of the abbreviations used in the code is listed as follows:

|    Abbr    | Description |
| -------- | ----------- |
| MN| Medium nitrogen |
| LN| Low nitrogen|
| HN| High nitrogen|
| MN-FACE| Medium nitrogen with FACE|
| LN-FACE| Low nitrogen with FACE|
| HN-FACE| High nitrogen with FACE|
| MN-FACE-vm_Nup64| Medium nitrogen with FACE by assuming maximum root nitrogen uptake rate under FACE being 64% of that under ambient conditions|
| LN-FACE-vm_Nup64| Low nitrogen with FACE by assuming maximum root nitrogen uptake rate under FACE being 64% of that under ambient conditions|
| HN-FACE-vm_Nup64| High nitrogen with FACE by assuming maximum root nitrogen uptake rate under FACE being 64% of that under ambient conditions|


###  7. Simulating effect of different leaf protein degradation activity 

To simulate effect of the OsNAP gene alteration shown in Liang et al. (2014), using the following commands:

```
c=num2cell(ones(1,49));c{1}='OsNAP-WT';c{2}=c{1};c{19}=1;c{7}=0.5;c{48}=1;c{27}=1;c{49}=1;c{5}=1;c{32}=1; RunModel(c,0,0,0,0);

c=num2cell(ones(1,49));c{1}='OsNAP-mutant';c{2}=c{1};c{19}=0.735;c{7}=2;c{48}=0.735;c{27}=0.735;c{49}=0.735;c{5}=0.735;c{32}=0.735; RunModel(c,0,0,0,0);

c=num2cell(ones(1,49));c{1}='OsNAP-RNAi';c{2}=c{1};c{19}=1.11;c{7}=0.25;c{48}=1.11;c{27}=1.11;c{49}=1.11;c{5}=1.11;c{32}=1.11; RunModel(c,0,0,0,0);
```

The meaning of the abbreviations used in the code is listed as follows:

|    Abbr    | Description |
| -------- | ----------- |
| OsNAP-WT| Wild type |
| OsNAP-mutant| ps1-D mutant|
| OsNAP-RNAi| OsNAP gene RNAi|


###  8. Simulating effect of different grain sucrose unloading activity

To simulate effect of the GIF1 gene alteration shown in Wang et al. (2008), using the following commands:

```
c=num2cell(ones(1,49));c{1}='GIF1-WT';c{2}=c{1};c{19}=90/160;c{48}=1/2;c{27}=1/2;c{49}=1/2;c{5}=1/2;c{32}=1/2;c{22}=2; RunModel(c,0,0,0,0);

c=num2cell(ones(1,49));c{1}='gif1';c{2}=c{1};c{19}=90/160;c{48}=1/2;c{27}=1/2;c{49}=1/2;c{5}=1/2;c{32}=1/2;c{22}=1;c{25}=0.92; RunModel(c,0,0,0,0);

c=num2cell(ones(1,49));c{1}='GIF1-OE';c{2}=c{1};c{19}=90/160;c{48}=1/2;c{27}=1/2;c{49}=1/2;c{5}=1/2;c{32}=1/2;c{22}=4;c{25}=1.04; RunModel(c,0,0,0,0);
```

To simulate effect of increasing the maximal grain sucrose unloading rate (c{22}) with sufficient large sink size, using the following commands:

```
c=num2cell(ones(1,49));c{1}='WT_s190';c{2}=c{1};c{19}=190/160;c{22}=1; RunModel(c,0,0,0,0);

c=num2cell(ones(1,49));c{1}='gif1_s190';c{2}=c{1};c{19}=190/160;c{22}=1/2; RunModel(c,0,0,0,0);

c=num2cell(ones(1,49));c{1}='GIF1_OE_s190';c{2}=c{1};c{19}=190/160;c{22}=2; RunModel(c,0,0,0,0);
```


The meaning of the abbreviations used in the code is listed as follows:

|    Abbr    | Description |
| -------- | ----------- |
| GIF1-WT| Wild type |
| gif1 | gif1 gene mutant|
| GIF1-OE| GIF-1 gene overexpression|
| WT_s190| Wild type with large sink (190 grains per ear)|
| gif1_s190| Hypothetical gif1 gene mutant with large sink (190 grains per ear)|
| GIF1_OE_s190| Hypothetical GIF-1 gene overexpression with large sink (190 grains per ear)|



###  9. Simulating effect of different grain starch synthesis activity

To simulate effect of Sh2r6hs transforming shown in Smidansky et al. (2003), using the following commands:

```
c=num2cell(ones(1,49));c{1}='Sh2r6hs-WT';c{2}=c{1};c{19}=90/160;c{48}=1/2;c{27}=1/2;c{49}=1/2;c{5}=1/2;c{32}=1/2;c{21}=1; RunModel(c,0,0,0,0);

c=num2cell(ones(1,49));c{1}='Sh2r6hs-
OE';c{2}=c{1};c{19}=1.13*90/160;c{44}=1.05;c{25}=1.007;c{5}=1.22/1.05*1/2;c{32}=1.22/1.05*1/2;c{27}=1.22/1.05*1/2;c{48}=1.22/1.05*1/2;c{49}=1.22/1.05*1/
2;c{21}=1.5; RunModel(c,0,0,0,0);
```


To simulate effect of increasing the maximal grain starch synthesis rate (c{21}) without change in plant size at flowering, using the following commands:

```
c=num2cell(ones(1,49));c{1}='Sh2r6hs-OE-AGPonly';c{2}=c{1};c{19}=90/160;c{48}=1/2;c{27}=1/2;c{49}=1/2;c{5}=1/2;c{32}=1/2;c{21}=1.5; RunModel(c,0,0,0,0);
```


To simulate effect of changing plant size at flowering without increase in the maximal grain starch synthesis rate (c{21}) , using the following commands:

```
c=num2cell(ones(1,49));c{1}='Sh2r6hs-OE-
plantSizeOnly';c{2}=c{1};c{19}=1.13*90/160;c{44}=1.05;c{25}=1.007;c{5}=1.22/1.05*1/2;c{32}=1.22/1.05*1/2;c{27}=1.22/1.05*1/2;c{48}=1.22/1.05*1/2;c{49}=1
.22/1.05*1/2;c{21}=1; RunModel(c,0,0,0,0);
```

The meaning of the abbreviations used in the code is listed as follows:

|    Abbr    | Description |
| -------- | ----------- |
| Sh2r6hs-WT| Wild type |
| Sh2r6hs-OE| Sh2r6hs gene overexpression |
| Sh2r6hs-OE-AGPonly| Sh2r6hs gene overexpression with assuming only grain AGP activity altered|
| Sh2r6hs-OE-plantSizeOnly| Sh2r6hs gene overexpression with assuming only plant size at flowering altered|



###  10. Perturbing all model parameters one by one

To test model behavior by perturbing all model parameters one by one, using the following commands:

```
para_list=...

    [[1,1];[1,1];[1,1];[1,1];[0.5,2];... % invariant: A_Jmax, PPFD,leaf_area, Phi_CO2

    [1,1];[0.5,2];[0.5,2];[0.5,2];[0.5,2];... % invariant: [NSC]inhibit_A_low
    
    [0.5,2];[0.5,2];[0.5,2];[1,1];[0.5,2];... % invariant: [O-N]leaf_protein_syn_low
    
    [1,1];[1,1];[0.5,2];[0.5,2];[0.5,2];... % invariant: gamma_star, seed_number
    
    [0.5,2];[0.5,2];[1,1];[1,1];[1,1];... % invariant: seed_width_max, seed_length_width_ratio, stem_starch_content
    
    [0.5,2];[0.5,2];[0.5,2];[0.5,2];[1,1];... % invariant: root_weight
    
    [1,1];[0.5,2];[0.5,2];[0.5,2];[0.5,2];... % invariant: soil_iN_conc
    
    [0.5,2];[0.5,2];[0.5,2];[0.5,2];[0.5,2];... % invariant: 
    
    [0.5,2];[1,1];[1,1];[1,1];[1,1];[1,1];[1,1]]; % invariant: tiller_number; CO2a; leaf_protein_content; seed_struct_N; stem_struct; stem_protein_content

seedNum=190/160;

candi_pos=[5,7:13,15,18:22,26:29,32:41]; % 28 variables

candi_name={'leafP2ON','leafONload','leafSUCload','leafIN2ON','leafTP2SUC','leafTP2STAR','leafSTAR2SUC','leafON2P','leafINunload',...
    'seedGrow','seedSTARstore','seedSUCunload','seedPstore','seedONunload','stemSUC2STAR','stemSTAR2SUC','stemON2P','stemP2ON',...
    'rootINabsorb','rootINload','rootONload','rootSUCunload','rootONunload','rootIN2ON','rootGrow','phloemRLS','phloemRLR','xylemONload'};

count=0;

for i=11:28

    c_list=ones(1,47);
    
	c_list(17)=seedNum;
	
    for val=[0.1:0.1:0.9,2:1:10]
    
        c_list(candi_pos(i))=val;
	
        fnOutNA=[candi_name{i},'_',num2str(val)];
	
        if exist([fnOutNA,'.txt'], 'file')==2 % when there are existing files
	
            disp('****************************************************************');
	    
            disp([fnOutNA,' has been done! Skipping it ...']);
	    
            continue
	    
        end
	
        input_info=num2cell([1,1,c_list]);
	
        input_info{1}=fnOutNA;
	
        input_info{2}=fnOutNA;
	
		count=count+1;
		
        fprintf('echo the %s (%d) is running...\n',fnOutNA,count);
	
        RunModel(input_info,0,0,0,0);
	
    end
    
end
```


###  11. In silico design of rice grain filling ideotype for super-high yield

To identify optimal parameter combinations maximizing grain yield with a genetic algorithm (on a Linux server), using the following commands (on a Linux server):

```
matlab -singleCompThread -nodesktop -nojvm -nosplash -nodisplay -r GA_for_evolutionaryPopulation;exit &
```

## How to Cite WACNI
Please cite the following manuscript:
>A critical role of stable grain filling rate in maximizing rice yield revealed by whole plant carbon nitrogen interaction modeling. [https://www.biorxiv.org/content/10.1101/2020.10.06.329029v2.full](https://www.biorxiv.org/content/10.1101/2020.10.06.329029v2.full). <br>
Tian-Gen Chang, Zhong-Wei Wei, Zai Shi, Yi Xiao, Honglong Zhao, Shuo-Qi Chang, Mingnan Qu, Qingfeng Song, Faming Chen, Fenfen Miao, Xin-Guang Zhu



## License

WACNI is licensed under the GNU General Public License v3.0. <br>
If you have any questions, please submit them on the [GitHub issues page](https://github.com/rootchang/WACNI-rice/issues).

## Contact

Tiangen Chang: changtiangen@gmail.com 
