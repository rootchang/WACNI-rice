# WACNI-rice
Carbon and nitrogen metabolism are the basis for plant growth and development. To predict grain yield formation directly from the molecular processes in different source, sink and transport organs, we develop a mechanistic model of Whole plAnt Carbon Nitrogen Interaction (WACNI). In contrast to previous methods that have only considered phloem sucrose transport from leaf to grain, or relied on preset sink growth pattern, our model considers both carbon and nitrogen, and kinetically simulates rates of major basic biochemical and physiological processes in a plant by ordinary differentiation equations. 

The model comprises five modules, i.e., root, leaf, grain, stem (including culm and sheath) and vascular transport system (xylem and phloem). Metabolites exchange between modules by trans-membrane transport. Fourteen types of biochemical and biophysical processes involved in different source, sink and transport organs are mathematically represented. These processes include assimilation, transport and utilization of six representative primary metabolites, i.e., triose phosphates (TP), sucrose (Suc), starch, inorganic nitrogen (I-N, including NH4+ and NO3-), free form of organic nitrogen (O-N, including amino acids and amides), and proteins. The model further incorporates the interaction between these metabolites and plant developmental processes, e.g. root growth, grain volume expansion (endosperm cell division), grain filling (starch and protein synthesis in endosperm), root senescence and leaf senescence. The apparent kinetic parameters for each process, initial metabolite concentrations, carbon and nitrogen mass in different organs at the flowering stage are input variables for the model. Finally, all types of plant-level physiological dynamics from flowering to harvest just emerge from this bottom-up simulation.

%%%%%%%% This folder contains the following essential files for running the WACNI model and generating in silico evolution populations %%%%%%%%%

%% RunModel.m

%% simulation_main.m

%% MyEventFunction.m

%% GA_for_evolutionaryPopulation.m

%% weather_input.txt

%% README.txt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%######  1. To run basic simulation for once, type in the following commands in MATLAB commands window:  ######%

c=num2cell(ones(1,49));

[Tt,simulation]=RunModel(c,0,0,0,0);

%######  2. To run basic simulation for once, showing plant carbon and nitrogen economy during grain filling, using:  ######%

c=num2cell(ones(1,49));

[Tt,simulation]=RunModel(c,0,0,0,1);

%######  3. To change one of the following parameters and run simulation for once (e.g. change c{7}: coef_K_protein_oN_leaf_cat to 0.5), using:  ######%

% c{1}: prefix of raw result file

% c{2}: prefix of file for ploting

% c{3}: coef_Vc_J

% c{4}: coef_PPFD   

% c{5}: coef_Area_leaf

% c{6}: coef_Phi_CO2

% c{7}: coef_K_protein_oN_leaf_cat

% c{8}: coef_lower_car_photo    

% c{9}: coef_K_oN_leaf_load_cat

% c{10}: coef_K_suc_leaf_load_cat

% c{11}: coef_K_N_assimi_leaf_cat

% c{12}: coef_K_TP_suc_leaf_cat

% c{13}: coef_K_TP_starch_leaf_cat    

% c{14}: coef_K_starch_suc_leaf_cat

% c{15}: coef_K_oN_protein_leaf_cat

% c{16}: coef_conc_lower_oN_to_leaf_protein

% c{17}: coef_K_iN_leaf_unload_cat

% c{18}: coef_gamma_star    

% c{19}: coef_grain_number

% c{20}: coef_K_grain_growth_cat

% c{21}: coef_K_grain_starch_store_cat

% c{22}: coef_K_suc_grain_unload_cat

% c{23}: coef_K_grain_protein_store_cat    

% c{24}: coef_K_oN_grain_unload_cat

% c{25}: coef_grain_width_max

% c{26}: coef_grain_length_width_ratio

% c{27}: coef_Starch_stem

% c{28}: coef_K_suc_starch_stem_cat    

% c{29}: coef_K_starch_suc_stem_cat

% c{30}: coef_K_oN_protein_stem_cat

% c{31}: coef_K_protein_oN_stem_cat

% c{32}: coef_Fresh_root

% c{33}: coef_conc_iN_soil    

% c{34}: coef_K_iN_absorb_cat

% c{35}: coef_K_iN_root_load_cat

% c{36}: coef_K_oN_root_load_cat

% c{37}: coef_K_suc_root_unload_cat

% c{38}: coef_K_oN_root_unload_cat    

% c{39}: coef_K_N_root_assimi_cat

% c{40}: coef_K_root_growth_cat

% c{41}: coef_R_l_s_phloem0

% c{42}: coef_R_r_l_phloem0

% c{43}: coef_K_oN_xylem_to_lphloem_cat      

% c{44}: coef_tillerNum

% c{45}: coef_CO2a

% c{46}: protein_content_leaf

% c{47}: grain_struct_N

% c{48}: stem_struct_mass

% c{49}: Protein_stem

c=num2cell(ones(1,49));

c{7}=0.5;

[Tt,simulation]=RunModel(c,0,0,0,0);


%######  4. To simulate different nitrogen treatments shown in Zhao et al. (2015), using:  ######%

c=num2cell(ones(1,49));

c{1}='MN';c{2}=c{1};c{33}=1;c{46}=1;c{49}=1;c{19}=1;c{44}=1;c{44}=11/12; RunModel(c,0,0,0,0);

c{1}='0N';c{2}=c{1};c{33}=0.2;c{46}=0.66;c{49}=0.66;c{19}=0.43/0.6;c{44}=0.6*11/12; RunModel(c,0,0,0,0);

c{1}='LN';c{2}=c{1};c{33}=0.5;c{46}=0.84;c{49}=0.84;c{19}=0.75/0.9;c{44}=0.9*11/12; RunModel(c,0,0,0,0);

c{1}='HN';c{2}=c{1};c{33}=1.5;c{46}=0.98;c{49}=0.98;c{19}=1.18/1.23;c{44}=1.23*11/12; RunModel(c,0,0,0,0);


%######  5. To simulate different lighting regimes shown in Kobata et al. (2000), using:  ######%

c=num2cell(ones(1,49));

c{1}='CK';c{2}=c{1}; RunModel(c,0,0,0,0);

c{1}='S10_25';c{2}=c{1}; RunModel(c,25,1,0,0);

c{1}='S10_50';c{2}=c{1}; RunModel(c,50,1,0,0);

c{1}='S10_75';c{2}=c{1}; RunModel(c,75,1,0,0);

c{1}='S10_25_Sp50';c{2}=c{1}; RunModel(c,25,1,50,0);

c{1}='S10_50_Sp50';c{2}=c{1}; RunModel(c,50,1,50,0);

c{1}='S10_75_Sp50';c{2}=c{1}; RunModel(c,75,1,50,0);


%######  6. To simulate different soil-nitrogen & air-CO2 concentrations shown in Kim et al. (2001), using:  ######%

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

% to simulate vm_Nupt=64% under FACE condition, using: 

c=num2cell(ones(1,49));c{1}='LN-FACE-
vm_Nup64';c{2}=c{1};c{19}=145/160;c{25}=1.09^(1/3);c{22}=1;c{34}=0.64;c{44}=0.89;c{5}=1.21/c{44};c{32}=1.18/c{44};c{48}=1.21/c{44};c{27}=1.21/c{44};c{49
}=0.96;c{19}=c{19}*0.95;c{46}=0.70;c{33}=0.5;c{45}=690/400; RunModel(c,0,0,0,0);

c=num2cell(ones(1,49));c{1}='MN-FACE-
vm_Nup64';c{2}=c{1};c{19}=145/160;c{25}=1.04^(1/3);c{22}=1;c{34}=0.64;c{44}=1.09;c{5}=1.45/c{44};c{32}=1.46/c{44};c{48}=1.45/c{44};c{27}=1.45/c{44};c{49
}=1.06;c{19}=c{19}*1.01;c{46}=0.79;c{33}=1;c{45}=690/400; RunModel(c,0,0,0,0);

c=num2cell(ones(1,49));c{1}='HN-FACE-
vm_Nup64';c{2}=c{1};c{19}=145/160;c{25}=0.93^(1/3);c{22}=1;c{34}=0.64;c{44}=1.15;c{5}=1.51/c{44};c{32}=1.35/c{44};c{48}=1.51/c{44};c{27}=1.51/c{44};c{49}=1.33;c{19}=c{19}*1.13;c{46}=1.01;c{33}=1.5;c{45}=690/400; RunModel(c,0,0,0,0);


%######  7. To simulate effect of the OsNAP gene shown in Liang et al. (2014), using:  ######%

c=num2cell(ones(1,49));c{1}='OsNAP-WT';c{2}=c{1};c{19}=1;c{7}=0.5;c{48}=1;c{27}=1;c{49}=1;c{5}=1;c{32}=1; RunModel(c,0,0,0,0);

c=num2cell(ones(1,49));c{1}='OsNAP-mutant';c{2}=c{1};c{19}=0.735;c{7}=2;c{48}=0.735;c{27}=0.735;c{49}=0.735;c{5}=0.735;c{32}=0.735; RunModel(c,0,0,0,0);

c=num2cell(ones(1,49));c{1}='OsNAP-RNAi';c{2}=c{1};c{19}=1.11;c{7}=0.25;c{48}=1.11;c{27}=1.11;c{49}=1.11;c{5}=1.11;c{32}=1.11; RunModel(c,0,0,0,0);


%######  8. To simulate effect of the GIF1 gene shown in Wang et al. (2008), using:  ######%

c=num2cell(ones(1,49));c{1}='GIF1-WT';c{2}=c{1};c{19}=90/160;c{48}=1/2;c{27}=1/2;c{49}=1/2;c{5}=1/2;c{32}=1/2;c{22}=2; RunModel(c,0,0,0,0);

c=num2cell(ones(1,49));c{1}='gif1';c{2}=c{1};c{19}=90/160;c{48}=1/2;c{27}=1/2;c{49}=1/2;c{5}=1/2;c{32}=1/2;c{22}=1;c{25}=0.92; RunModel(c,0,0,0,0);

c=num2cell(ones(1,49));c{1}='GIF1-OE';c{2}=c{1};c{19}=90/160;c{48}=1/2;c{27}=1/2;c{49}=1/2;c{5}=1/2;c{32}=1/2;c{22}=4;c{25}=1.04; RunModel(c,0,0,0,0);

% to simulate effect of increasing the maximal grain sucrose unloading rate (c{22}) with sufficient large sink size, using: 
c=num2cell(ones(1,49));c{1}='WT_s190';c{2}=c{1};c{19}=190/160;c{22}=1; RunModel(c,0,0,0,0);

c=num2cell(ones(1,49));c{1}='gif1_s190';c{2}=c{1};c{19}=190/160;c{22}=1/2; RunModel(c,0,0,0,0);

c=num2cell(ones(1,49));c{1}='GIF1_OE_s190';c{2}=c{1};c{19}=190/160;c{22}=2; RunModel(c,0,0,0,0);

%######  9. To simulate effect of Sh2r6hs transforming shown in Smidansky et al. (2003), using:  ######%

c=num2cell(ones(1,49));c{1}='Sh2r6hs-WT';c{2}=c{1};c{19}=90/160;c{48}=1/2;c{27}=1/2;c{49}=1/2;c{5}=1/2;c{32}=1/2;c{21}=1; RunModel(c,0,0,0,0);

c=num2cell(ones(1,49));c{1}='Sh2r6hs-
OE';c{2}=c{1};c{19}=1.13*90/160;c{44}=1.05;c{25}=1.007;c{5}=1.22/1.05*1/2;c{32}=1.22/1.05*1/2;c{27}=1.22/1.05*1/2;c{48}=1.22/1.05*1/2;c{49}=1.22/1.05*1/
2;c{21}=1.5; RunModel(c,0,0,0,0);

% to simulate effect of increasing the maximal grain starch synthesis rate (c{21}) without change in plant size at flowering, using: 

c=num2cell(ones(1,49));c{1}='Sh2r6hs-OE-AGPonly';c{2}=c{1};c{19}=90/160;c{48}=1/2;c{27}=1/2;c{49}=1/2;c{5}=1/2;c{32}=1/2;c{21}=1.5; RunModel(c,0,0,0,0);

% to simulate effect of changing plant size at flowering without increase in the maximal grain starch synthesis rate (c{21}) , using: 

c=num2cell(ones(1,49));c{1}='Sh2r6hs-OE-
plantSizeOnly';c{2}=c{1};c{19}=1.13*90/160;c{44}=1.05;c{25}=1.007;c{5}=1.22/1.05*1/2;c{32}=1.22/1.05*1/2;c{27}=1.22/1.05*1/2;c{48}=1.22/1.05*1/2;c{49}=1
.22/1.05*1/2;c{21}=1; RunModel(c,0,0,0,0);

%######  10. To perform sensitivity analysis of model parameters, using:  ######%

para_list=[[1,1];[1,1];[1,1];[1,1];[0.5,2];... % invariant: A_Jmax, PPFD,leaf_area, Phi_CO2
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


%######  11. To identify optimal parameter combinations maximizing grain yield with a genetic algorithm (on a Linux server), using:  ######%

matlab -singleCompThread -nodesktop -nojvm -nosplash -nodisplay -r GA_for_evolutionaryPopulation;exit &

