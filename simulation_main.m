function mc = simulation_main(t,c)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulation_main.m: This function defines kinetic property of reactions and diffusions, and simulate plant physiolocial change from flowering to harvest.
%
%  Input:
%   t: Current time after flowering
%   c: input numerical array (len=600)
%
%  Output:
%   mc: differential of parameters in c at time point t
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% global
global shading_span;
global SPACING;

    coef_Vc_J=c(458);
    coef_PPFD=c(459);
%     coef_Area_leaf=c(460);
    coef_Phi_CO2=c(461);
    coef_K_Pro_oN_leaf_cat=c(462);
    coef_conc_lower_NSC_inhibit_photosyn_leaf=c(463); % NSC, carbohydrate
    coef_K_oN_leaf_load_cat=c(464);
    coef_K_Suc_leaf_load_cat=c(465);
    coef_K_N_assimi_leaf_cat=c(466);
    coef_K_TP_Suc_leaf_cat=c(467);
    coef_K_TP_Star_leaf_cat=c(468);
    coef_K_Star_Suc_leaf_cat=c(469);
    coef_K_oN_Pro_leaf_cat=c(470);
    coef_conc_lower_oN_to_leaf_Pro=c(471);
    coef_K_iN_leaf_unload_cat=c(472);
    coef_gamma_star=c(473);
    
    % grain
%     coef_grainNum=c(474);
    coef_K_grain_grow_cat=c(475);
    coef_K_grain_Star_store_cat=c(476);
    coef_K_Suc_grain_unload_cat=c(477);
    coef_K_grain_Pro_store_cat=c(478);
    coef_K_oN_grain_unload_cat=c(479);
%     coef_grain_wid_max=c(480);
%     coef_grain_len_wid_ratio=c(481);
    
    % stem
%     coef_Star_stem=c(482);
    coef_K_Suc_Star_stem_cat=c(483);
    coef_K_Star_Suc_stem_cat=c(484);
    coef_K_oN_Pro_stem_cat=c(485);
    coef_K_Pro_oN_stem_cat=c(486);
    
    % Root
%     coef_Fresh_root=c(487);
%     coef_conc_iN_soil=c(488);
    coef_K_iN_absorb_cat=c(489);
    coef_K_iN_root_load_cat=c(490);
    coef_K_oN_root_load_cat=c(491);
    coef_K_Suc_root_unload_cat=c(492);
    coef_K_oN_root_unload_cat=c(493);
    coef_K_N_root_assimi_cat=c(494);
    coef_K_root_grow_cat=c(495);
    
    % flow resistance
%     coef_R_l_s_phloem0=c(496);
%     coef_R_r_l_phloem0=c(497);
    coef_K_oN_xylem_to_lphloem_cat=c(498);
    
    % others
%     coef_tillerNum=c(499);
%     coef_CO2a=c(500);
    
    % input variables
%     global weather_time_list;
%     global T_air_list;
    global PPFD_list;
    global VPD_leaf_list;
    
    % initial variables
    Suc_molWeight=c(56);
    N_molWeight=c(57);
    Pro_molWeight=c(58);
    N_to_Pro_mass=c(59);
    root_struc_N=c(60);
    grain_struc_N=c(61);
	grow_efficiency=c(62);
    grow_efficiency_grain=c(63);
    tillerNum=c(64);
    plantDens=c(65);
    coef_Suc_maintain_leaf=c(66);
    coef_Suc_maintain_grain=c(67);
    coef_Suc_maintain_root=c(68);
    coef_Suc_maintain_stem=c(69);
    Vol_Fresh_ratio_leaf=c(70);
    Thick_leaf=c(71);
    CO2_a=c(72);
    lat_Shanghai=c(73);
    DOY_start=c(74);
    Area_leaf0=c(75);
%     Vol_leaf0=c(76);
%     photosyn_Pro_leaf0=c(77);
%     Vol_soil=c(78);
%     Rad_root=c(79);
    Vol_Fresh_ratio_root=c(80);
    Fresh_dry_ratio_root=c(81);
%     Absorb_area_ratio=c(82);
%     Soil_water_content=c(83);
%     Vol_root0=c(84);
%     Dry_root0=c(85);
    Vol_xylem=c(86);
    Vol_rphloem=c(87);
    Vol_lphloem=c(88);
    Vol_sphloem=c(89);
    R_l_s_phloem0=c(90);
    R_r_l_phloem0=c(91);
    Vol_Fresh_ratio_grain=c(92);
    Fresh_new_dry_ratio_grain=c(93);
    pho_Pro=c(94);
    pho_Star=c(95);
    grainNum=c(96);
    grain_len_wid_ratio=c(97);
    Vol_grain_max=c(98);
%     Area_all_grain_max=c(99);
    Vol_stem=c(100);
    
    %% variants set
    Vol_root=c(1);
    Dry_root=c(2);
    conc_Suc_root=c(3);
    conc_iN_root=c(4);
    conc_oN_root=c(5);
    conc_iN_soil=c(6);
%     Vol_soil_water=c(7);
    
    Area_leaf=c(101);
    Vol_leaf=c(102);
    conc_Suc_leaf=c(103);
    conc_iN_leaf=c(104);
    conc_oN_leaf=c(105);
    totMol_Pro_leaf=c(106);
    Star_leaf=c(107);
    conc_TP_leaf=c(108);
    
    Vol_grain=c(201);
    Area_grain=c(202);
    conc_Suc_grain=c(203);
    conc_oN_grain=c(204);
    store_Pro_grain=c(205);
    store_Star_grain=c(206);
    
    conc_iN_xylem=c(301);
    conc_Suc_rphloem=c(302);
    conc_Suc_lphloem=c(303);
    conc_Suc_sphloem=c(304);
    conc_oN_rphloem=c(305);
    conc_oN_lphloem=c(306);
    conc_oN_sphloem=c(307);
    conc_oN_xylem=c(308);
    
    conc_Suc_stem=c(401);
    conc_oN_stem=c(402);
    Star_stem=c(403);
    Pro_stem=c(404);
    Vol_soluble_stem=c(405);
%% local general
    % coef.
    coef_unload=0.06; %mol Suc/mol metabolite transported
    coef_load=0.06; %mol Suc/mol metabolite transported
    coef_Suc_to_Star=0.03; %mol Suc/mol Suc-Star inter-convertion
    coef_iN_to_oN=0.11; % mol Suc/mol N
    coef_C_skeleton_iN2oN=0.417; % 5/12
    coef_oN_to_Pro=0.3; % 17.8/60 mol Suc/mol N->Pro, 17.8mol ATP/mol N
    coef_Pro_decomp=0.05; % mol Suc/mol N
    coef_iN_absorb=0.025; %mol Suc/mol N
    
    %% root
    conc_Suc_grow_root=10; % lower than this, grow ceases.
    K_root_grow_cat=coef_K_root_grow_cat*0.2/86400; %g/g/s, assume max relative grow rate is 20% (m3/m3/d)
    Km_root_grow=50-conc_Suc_grow_root; % mol Suc/m3

    conc_Suc_sene_root=8; % lower than this, root senescence was accelerated
    relative_sene_rate_normal=0.03;
    K_root_sene_cat=relative_sene_rate_normal/86400; % assume regular relative senescence rate is 3% (m3/m3/d)
    coef_exp_sene_root=0.24;

    conc_Suc_absorb_root=0; % lower than this, root nitrogen absorption ceases
    conc_lower_all_N_inhibit_absorb_root=0; % mol/m3
    conc_upper_all_N_inhibit_absorb_root=100; % mol/m3
    K_absorb_HATS_cat=coef_K_iN_absorb_cat*3.65*10^-9; % mol/gDW/s
    Km_absorb_HATS=0.075; % mol/m3
    K_absorb_LATS_cat=coef_K_iN_absorb_cat*5.5*10^-10; %m3/gDW/s
    Km_Suc_for_root_absorb=10; % mol/m3
    
    K_iN_root_load_cat=coef_K_iN_root_load_cat*2.5*10^-9; % mol/g/s
    Ke_iN_root_load=1; % mol/m3
    Km_iN_root_load=10; % mol/m3
    
    K_oN_root_transfer_coef=3.2;
    K_oN_root_load_cat=coef_K_oN_root_load_cat*K_oN_root_transfer_coef*10^-9; % mol/g/s
	Ke_oN_root_load=1; % mol/m3
    Km_oN_root_load=10; % mol/m3
    
    K_Suc_root_unload_cat=coef_K_Suc_root_unload_cat*1*10^-8; % mol Suc/g/s
    Km_Suc_root_unload=300; % mol/m3
    Ke_Suc_root_unload=0.1; % mol/m3
    
    K_oN_root_unload_cat=coef_K_oN_root_unload_cat*K_oN_root_transfer_coef*10^-9; % mol/g/s
    Km_oN_root_unload=100; % mol/m3
    Ke_oN_root_unload=0.1; % mol/m3

    K_N_root_assimi_cat=coef_K_N_root_assimi_cat*1.25*10^-9; % mol N/g/s
    Km_N_root_assimi_N=20; % mol/m3
    Km_N_root_assimi_Suc=10; % mol/m3
    
    %% leaf
    T_hour=t/3600;
    timestamp=max(1,ceil(T_hour)); % timestamp is 0,1,2,3,..., can be used as index
    coef_k_Beer=1;
    if SPACING && T_hour>=24*shading_span
        plantDens=plantDens*(100-SPACING)/100;
    end

    k_Beer=coef_k_Beer*0.37; 
    LAI0=tillerNum*plantDens*Area_leaf0;
    
    LAI=max(0.01,LAI0*Area_leaf/(Area_leaf0));
    PPFD=coef_PPFD*PPFD_list(timestamp)*10^-6; % mol/m2 ground/s, start from 6:00 AM
    VPD_leaf=VPD_leaf_list(timestamp);
    Vcmax0=coef_Vc_J*80*10^-6; %mol CO2/m2/s
    Jmax0=coef_Vc_J*160*10^-6; %mol CO2/m2/s
    gamma_star=coef_gamma_star*38.6; % ubar=ppm
    Km_CO2_assimi_leaf=418; % ubar=ppm
    K_N_assimi_leaf_cat=coef_K_N_assimi_leaf_cat*6.18*10^-7; % *0.17/22 mol N/m2/s, set as 0.17/22*Vcmax (NO3 assimi rate), CO2 based
    Km_N_assimi_leaf=150; % mol N/m3
    
    f_a=0.425; % diffusion property (The proportion of attenuated radiation that reaches the surface as diffuse radiation), 0.4-0.45
    air_trans=0.75; % air transparantness
    pho_cb=0.029; % canopy reflectrum coef. for beam PAR
    pho_cd=0.036; % canopy reflectrum coef. for diffusion PAR
    sigma=0.15; % leaf scatter coef. of PAR
    k_d=0.78; % diffsuion PAR extinction coef.
    k_d1=0.9*k_d; % diffusion and scattered PAR extinction coef.
    Area_ground=1/(tillerNum*plantDens); % m2; area of ground that one tiller can cover
    
    local_time_in_24=mod(t/3600+6,24);
    DOY_current=floor((t/3600+6)/24)+DOY_start;
    A_delta=-0.4084*cos(2*pi*(DOY_current+10)/365); % -23.4*pi/180*cos(2*pi*(DOY_current+10)/365)
    sin_beta=sin(lat_Shanghai)*sin(A_delta)+cos(lat_Shanghai)*cos(A_delta)*cos(pi*(local_time_in_24-12)/12);
    if PPFD>0 && sin_beta>0
        k_b=k_Beer/sin_beta; % beam PAR extinction coef.
        k_b1=0.9*k_b; % beam and scatter PAR extinction coef.
        
        f_d=(1-air_trans^(1/sin_beta))/(1+air_trans^(1/sin_beta)*(1/f_a-1));
        I_d=PPFD*f_d;
        I_b=PPFD*(1-f_d);
        I_canopy=(1-pho_cb)*I_b*(1-exp(-k_b1*LAI))+(1-pho_cd)*I_d*(1-exp(-k_d1*LAI)); % total canopy absorbed light, mol/m2 ground/s
        I_sun=I_b*(1-sigma)*(1-exp(-k_b*LAI))+...
            I_d*(1-pho_cd)*(1-exp(-(k_d1+k_b)*LAI))*k_d1/(k_d1+k_b)+...
            I_b*((1-pho_cb)*(1-exp(-(k_b1+k_b)*LAI))*k_b1/(k_b1+k_b)-(1-sigma)*(1-exp(-2*k_b*LAI))/2); % sunlit leaf absorbed light, mol/m2 ground/s
        I_shade=I_canopy-I_sun; % shaded leaf absorbed light, mol/m2 ground/s
        LAI_sun=(1-exp(-k_b*LAI))/k_b; % sunlit leaf area index
        LAI_shade=LAI-LAI_sun;
        Vcmax0_sun=Vcmax0*LAI_sun;
        Jmax0_sun=Jmax0*LAI_sun;
        Vcmax0_shade=Vcmax0*LAI_shade;
        Jmax0_shade=Jmax0*LAI_shade;
        K_N_assimi_leaf_cat_sun=K_N_assimi_leaf_cat*LAI_sun;
        K_N_assimi_leaf_cat_shade=K_N_assimi_leaf_cat*LAI_shade;
    end
    conc_Pro_max_promote_leaf_photosyn=2200; %mol N/m3
    conc_lower_NSC_inhibit_photosyn_leaf=coef_conc_lower_NSC_inhibit_photosyn_leaf*400; % mol Suc equiv. carbohydrate/m3
    conc_upper_NSC_inhibit_photosyn_leaf=1900;
    
    g0=0.01; %  mol CO2/m2/s, from Leuning 1992
    alpha_gs=20; % coef, ave. value from Leuning 1992
    VPD_0=0.35; % kPa, from Leuning 1992
    
    K_TP_Suc_leaf_cat=coef_K_TP_Suc_leaf_cat*3.3*10^-6; % mol Suc/m3/s
    Ke_TP_Suc_leaf=100;
    Km_TP_Suc_leaf=2; % mol Suc/m3
    K_TP_Star_leaf_cat=coef_K_TP_Star_leaf_cat*8.3*10^-7; % mol Suc/m3/s
    Km_TP_Star_leaf=2; % mol Suc/m3
    
    K_Star_Suc_leaf_cat=coef_K_Star_Suc_leaf_cat*8.3*10^-6;
    conc_Suc_leaf_Star_deg_low=120; % mol Suc/m3
    Km2_Star_Suc_leaf=15; % mol Suc/m3
    
    K_oN_Pro_leaf_cat=coef_K_oN_Pro_leaf_cat*6.21*10^-7; % mol N/m3/s
    conc_lower_oN_to_leaf_Pro=coef_conc_lower_oN_to_leaf_Pro*100; % mol N/m3
    Km_oN_Pro_leaf=45; % mol N/m3
    K_Pro_oN_leaf_cat=coef_K_Pro_oN_leaf_cat*9.19*10^-8; % mol N/mol N equiv. Pro/s
    Km2_Pro_oN_leaf=500; % mol N/m3
    
    K_Pro_turnover_leaf_cat=0.1/86400; % assume leaf Pro turnover rate 10%/d
    
    conc_all_N_sene_leaf=2000; % mol N/m3
    conc_iN_oN_Retrieve_upper_leaf=200; % mol N/m3, if free N concentration exceeds this value, N from sene. leaf will NOT be retrieved
    conc_iN_oN_Retrieve_lower_leaf=150; % mol N/m3, if free N concentration lower than this value, N from sene. leaf will ALL be retrieved
    K_leaf_sene_cat=0.015/86400; % assume normal relative senescence rate is 1.5% (m3/m3/d)
    coef_exp_sene_leaf_new=0.0015; 
    
    K_iN_leaf_unload_cat=coef_K_iN_leaf_unload_cat*6*10^-8;
    Ke_iN_leaf_unload=10; % mol N/m3
    Km_iN_leaf_unload=10; % mol N/m3
    
    K_oN_leaf_unload_cat=7.5*10^-8;
    Ke_oN_leaf_unload=15; % mol N/m3
    Km_oN_leaf_unload=10; % mol N/m3

    K_oN_leaf_load_cat=coef_K_oN_leaf_load_cat*4.8*10^-7; % mol/m2/s
    Ke_oN_leaf_load=3; % mol N/m3
    Km_oN_leaf_load=50; % mol N/m3
    
    K_Suc_leaf_load_cat=coef_K_Suc_leaf_load_cat*1.2*10^-6; % mol Suc/m2/s
    Ke_Suc_leaf_load=7; % mol Suc/m3
    Km_Suc_leaf_load=100; % mol Suc/m3
    
    %% ear
    K_grain_grow_cat=coef_K_grain_grow_cat*4*10^-9; % m3/m2 surface area/s, single grain
    Km_grain_grow_oN=50;
    Km_grain_grow_Suc=80;
    
    K_grain_Star_store_cat=coef_K_grain_Star_store_cat*0.0095; % mol Suc equiv. Star/m3 soluble space/s, all grains
    Km_grain_Star_store_Suc=2*Km_grain_grow_Suc;
    K_grain_Pro_store_cat=coef_K_grain_Pro_store_cat*0.0034; % mol N equiv. Pro/m3 soluble space/s
    Km_grain_Pro_store_oN=2*Km_grain_grow_oN;
    
    K_oN_grain_unload_cat=coef_K_oN_grain_unload_cat*2*10^-6; % mol N/m2 surface area/s
    Ke_oN_grain_unload=1; 
    Km_oN_grain_unload=100;
    
    K_Suc_grain_unload_cat=coef_K_Suc_grain_unload_cat*5*10^-6; % mol Suc/m2/s
    Ke_Suc_grain_unload=1;
    Km_Suc_grain_unload=400;
    e_glume_restrict=0.3;
    
    %% stem
    K_Suc_Star_stem_cat=coef_K_Suc_Star_stem_cat*0.0125; % mol Suc/m3 soluble space/s
    conc_lower_Suc_to_stem_Star=500;
    Km_Suc_Star_stem=160-conc_lower_Suc_to_stem_Star;
    K_Star_Suc_stem_cat=coef_K_Star_Suc_stem_cat*0.01; % mol Suc/mol Suc equiv. Star/s
    Km2_Star_Suc_stem=600;
    
    K_oN_Pro_stem_cat=coef_K_oN_Pro_stem_cat*3.4*10^-3; % mol N/m3 soluble space/s
    conc_lower_oN_to_stem_Pro=300;
    Km_oN_Pro_stem=500-conc_lower_oN_to_stem_Pro;
    K_Pro_oN_stem_cat=coef_K_Pro_oN_stem_cat*2.1*10^-3; % mol N/mol N equiv. Pro/s
    Km2_Pro_oN_stem=Km_oN_Pro_stem;
    
    R_Suc_lphloem_stem=7.4*10^7; %s/m3, diffusion resistance, [Suc] gradient based
    R_oN_lphloem_stem=3.2*10^7; %s/m3, diffusion resistance, [O-N] gradient based
    
    %% flow
    K_oN_xylem_to_lphloem_cat=coef_K_oN_xylem_to_lphloem_cat*1*10^-9;
    Ke_oN_xylem_to_lphloem=45;
    Km_oN_xylem_to_lphloem=10;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% metabolism in each organ and ODE development %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% root metabolism
    % N absorption (Suc promote both, organic N inhibit HATS), mol N/s
    conc_all_N_root=conc_oN_root+conc_iN_root;
    coef_all_N_inhibit_root_absorb=1;
    if conc_all_N_root>conc_lower_all_N_inhibit_absorb_root && conc_all_N_root<conc_upper_all_N_inhibit_absorb_root
        coef_all_N_inhibit_root_absorb=-conc_all_N_root/conc_upper_all_N_inhibit_absorb_root+1;
    elseif conc_all_N_root>=conc_upper_all_N_inhibit_absorb_root
        coef_all_N_inhibit_root_absorb=0;
    end
    v_iN_absorb=0;
    if conc_Suc_root>conc_Suc_absorb_root
        v_iN_absorb=Dry_root*(conc_Suc_root-conc_Suc_absorb_root)/(Km_Suc_for_root_absorb+conc_Suc_root) ...
            *(coef_all_N_inhibit_root_absorb*K_absorb_HATS_cat*conc_iN_soil/(Km_absorb_HATS+conc_iN_soil) ...
            +K_absorb_LATS_cat*conc_iN_soil);
    end
    % N assimi
    v_N_root_assimi=Dry_root*K_N_root_assimi_cat*conc_iN_root/(Km_N_root_assimi_N+conc_iN_root) ...
        *conc_Suc_root/(Km_N_root_assimi_Suc+conc_Suc_root); % mol N/s
    % grow (Suc promote grow), m3/s
    if conc_Suc_root>conc_Suc_grow_root 
        v_root_grow=Vol_root*K_root_grow_cat*(conc_Suc_root-conc_Suc_grow_root)/(Km_root_grow+(conc_Suc_root-conc_Suc_grow_root)); %m3/m3/s
    else
        v_root_grow=0;
    end
    % sene.(Suc starvation accelerate sene.), m3/s
    if conc_Suc_root<conc_Suc_sene_root
        v_root_sene=Vol_root*K_root_sene_cat*exp(-coef_exp_sene_root*(conc_Suc_root-conc_Suc_sene_root)); %m3/m3/s
    else
        v_root_sene=Vol_root*K_root_sene_cat;
    end
    % load, mol N/s
    v_iN_root_load=Dry_root*K_iN_root_load_cat...
        *max(0,conc_iN_root-conc_iN_xylem/Ke_iN_root_load)/(Km_iN_root_load+conc_iN_root);
    v_oN_root_load=Dry_root*K_oN_root_load_cat...
        *max(0,(conc_oN_root-conc_oN_xylem/Ke_oN_root_load))/(Km_oN_root_load+conc_oN_root);
    % unload, mol/s
    v_oN_root_unload=Dry_root*K_oN_root_unload_cat...
        *max(0,(conc_oN_rphloem-conc_oN_root/Ke_oN_root_unload))/(Km_oN_root_unload+conc_oN_rphloem);
    v_Suc_root_unload=Dry_root*K_Suc_root_unload_cat...
        *max(0,(conc_Suc_rphloem-conc_Suc_root/Ke_Suc_root_unload))/(Km_Suc_root_unload+conc_Suc_rphloem);

    % ODE in root
    delta_Vol_root=v_root_grow-v_root_sene;
    delta_Dry_root=Fresh_dry_ratio_root*delta_Vol_root/Vol_Fresh_ratio_root*10^6; % g
    
    grow_Weight_root=Fresh_dry_ratio_root*v_root_grow/Vol_Fresh_ratio_root*10^6; % g
    grow_Rd_root=grow_Weight_root*(1/grow_efficiency-1)/Suc_molWeight;
    grow_Suc_cost_root=grow_Rd_root + grow_Weight_root*(1-root_struc_N*N_to_Pro_mass)/Suc_molWeight;
    maintain_cost=coef_Suc_maintain_root*Vol_root;
    if conc_Suc_root<0.001
        maintain_cost=0;
    end
    delta_conc_Suc_root=(v_Suc_root_unload+conc_Suc_root*v_root_sene... % unload + senescence release
        -coef_unload*(v_oN_root_unload+v_Suc_root_unload)-coef_load*(v_iN_root_load+v_oN_root_load) ... % (un)load
        -maintain_cost-coef_iN_absorb*v_iN_absorb-(coef_iN_to_oN+coef_C_skeleton_iN2oN)*v_N_root_assimi ... 
        -conc_Suc_root*v_root_grow-grow_Suc_cost_root)/Vol_root; % grow cost (volume increase dilution + struc built cost)
    delta_conc_iN_root=(v_iN_absorb+conc_iN_root*v_root_sene ... % absorption + senescence release
        -v_iN_root_load-v_N_root_assimi ... % load + assimi
        -conc_iN_root*v_root_grow)/Vol_root; % grow volume increase dilution
    delta_conc_oN_root=(v_oN_root_unload-v_oN_root_load+conc_oN_root*v_root_sene ... % (un)load + senescence release
        +v_N_root_assimi ... % N assimi
        -conc_oN_root*v_root_grow-grow_Weight_root*root_struc_N/N_molWeight)/Vol_root; % grow volume increase dilution + struc built cost
    delta_Vol_soil_water=0;
    delta_conc_iN_soil=0;
    
    % ODE for root record paras.
    delta_Suc_in_root=v_Suc_root_unload;
    delta_oN_shoot_root_netFlux=v_oN_root_unload-max(0,v_oN_root_load-v_N_root_assimi);
    delta_N_shoot_root_netFlux=v_oN_root_unload - v_iN_root_load - v_oN_root_load;
    delta_iN_out_root=v_iN_root_load;
    delta_Vol_sene_root=v_root_sene;
    delta_Vol_grow_root=v_root_grow;
    delta_Rd_root=coef_unload*(v_oN_root_unload+v_Suc_root_unload) ... % unload cost
        +coef_load*(v_iN_root_load+v_oN_root_load) ... % load cost
        +coef_Suc_maintain_root*Vol_root+coef_iN_absorb*v_iN_absorb ...% % maintain and absorption cost
        +coef_iN_to_oN*v_N_root_assimi... % N assimi cost
        +grow_Rd_root; % grow Rd, % mol Suc/s
    delta_N_assimi_root=v_N_root_assimi;
    delta_N_absorb_root=v_iN_absorb;
    
    delta_oN_import_root=v_oN_root_unload;
    delta_oN_export_root=v_oN_root_load;
    delta_iN_export_root=v_iN_root_load;
    
    %% leaf metabolism
    % effect of leaf Pro and carbohydrate on photosynsynthesis
    conc_Pro_leaf=totMol_Pro_leaf/Vol_leaf; %molN/m3
    coef_Pro_promote_leaf_photosyn=min(conc_Pro_leaf,conc_Pro_max_promote_leaf_photosyn)/conc_Pro_max_promote_leaf_photosyn;
    conc_Star_leaf=Star_leaf/Vol_leaf; % mol Suc/m3
    conc_NSC_leaf=conc_Suc_leaf+conc_Star_leaf;
    if conc_NSC_leaf>conc_lower_NSC_inhibit_photosyn_leaf && conc_NSC_leaf<conc_upper_NSC_inhibit_photosyn_leaf
        coef_NSC_inhibit_leaf_photosyn=1-(conc_NSC_leaf-conc_lower_NSC_inhibit_photosyn_leaf)/ ...
            (conc_upper_NSC_inhibit_photosyn_leaf-conc_lower_NSC_inhibit_photosyn_leaf); % range: 0-1
    elseif conc_NSC_leaf>=conc_upper_NSC_inhibit_photosyn_leaf
        coef_NSC_inhibit_leaf_photosyn=0;
    else
        coef_NSC_inhibit_leaf_photosyn=1;
    end
    % night-t stored Star remobilization
    if conc_Suc_leaf<conc_Suc_leaf_Star_deg_low && PPFD<20*10^-6 % Star decomposes in the night, when Suc concentration is low
        v_Star_Suc_leaf=Area_leaf*K_Star_Suc_leaf_cat*(conc_Suc_leaf_Star_deg_low-conc_Suc_leaf)/...
            (conc_Suc_leaf_Star_deg_low+conc_Suc_leaf)*conc_Star_leaf/(conc_Star_leaf+Km2_Star_Suc_leaf); 
    else
        v_Star_Suc_leaf=0;
    end
    % photosyn_Pro and organic N conversion
    if conc_oN_leaf>conc_lower_oN_to_leaf_Pro && t<20*86400
        v_oN_Pro_leaf=Area_leaf*K_oN_Pro_leaf_cat*(conc_oN_leaf-conc_lower_oN_to_leaf_Pro)/...
            (conc_oN_leaf+Km_oN_Pro_leaf); 
    elseif conc_oN_leaf<conc_lower_oN_to_leaf_Pro
        v_oN_Pro_leaf=Area_leaf*K_Pro_oN_leaf_cat*(conc_oN_leaf-conc_lower_oN_to_leaf_Pro)/...
            (conc_oN_leaf+Km_oN_Pro_leaf)...
            *conc_Pro_leaf/(Km2_Pro_oN_leaf+conc_Pro_leaf);
    else
            v_oN_Pro_leaf=0;
    end
    % senescence
    conc_all_N_leaf=conc_Pro_leaf+conc_oN_leaf+conc_iN_leaf;
    conc_iN_oN_leaf=conc_oN_leaf+conc_iN_leaf;
    if conc_all_N_leaf<conc_all_N_sene_leaf
        v_leaf_sene=Area_leaf*K_leaf_sene_cat*exp(-coef_exp_sene_leaf_new*(conc_all_N_leaf-conc_all_N_sene_leaf)); % m2/s
    else
        v_leaf_sene=Area_leaf*K_leaf_sene_cat;
    end
    delta_Area_leaf=-v_leaf_sene;
    delta_Vol_leaf=delta_Area_leaf*Thick_leaf*Vol_Fresh_ratio_leaf;
    % load and unload
    v_oN_leaf_load=Area_leaf*K_oN_leaf_load_cat*max(0,(conc_oN_leaf-conc_oN_lphloem/Ke_oN_leaf_load))/ ...
        (Km_oN_leaf_load+conc_oN_leaf);
    v_Suc_leaf_load=Area_leaf*K_Suc_leaf_load_cat*max(0,(conc_Suc_leaf-conc_Suc_lphloem/Ke_Suc_leaf_load))/ ...
        (Km_Suc_leaf_load+conc_Suc_leaf);
    
    v_iN_leaf_unload=Area_leaf*K_iN_leaf_unload_cat...
        * max(0,conc_iN_xylem-conc_iN_leaf/Ke_iN_leaf_unload)/(Km_iN_leaf_unload+conc_iN_xylem);
    v_oN_leaf_unload=Area_leaf*K_oN_leaf_unload_cat...
        * max(0,conc_oN_xylem-conc_oN_leaf/Ke_oN_leaf_unload)/(Km_oN_leaf_unload+conc_oN_xylem);
    
    % Rd
    v_Rd_leaf=coef_unload*(v_iN_leaf_unload+v_oN_leaf_unload)+coef_load*(v_oN_leaf_load+v_Suc_leaf_load) ... % (un)load cost
        +coef_Suc_maintain_leaf*Vol_leaf ... % maintain cost
        +coef_oN_to_Pro*max(0,v_oN_Pro_leaf)+coef_Pro_decomp*(max(0,-v_oN_Pro_leaf) + conc_Pro_leaf*delta_Vol_leaf)... % AA <-> pro inter-conversion
        +(coef_oN_to_Pro+coef_Pro_decomp)*totMol_Pro_leaf*K_Pro_turnover_leaf_cat... % Pro turnover
        +coef_Suc_to_Star*(v_Star_Suc_leaf + conc_Star_leaf*delta_Vol_leaf); % Star <-> Suc inter-conversion cost, mol Suc/s; 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % photosynsynthesis, Rd and stomata conductance
    % using sun-shade model, calculating sunlit and shaded part separately
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if PPFD>0 && sin_beta>0
        alpha=coef_Phi_CO2*0.361*coef_Pro_promote_leaf_photosyn;  % 0.361=0.85*0.85*0.5
        theta=0.7*coef_Pro_promote_leaf_photosyn;
        Jmax_sun=Jmax0_sun*coef_Pro_promote_leaf_photosyn*coef_NSC_inhibit_leaf_photosyn;
        Vcmax_sun=Vcmax0_sun*coef_Pro_promote_leaf_photosyn*coef_NSC_inhibit_leaf_photosyn;
        J_photosyn_electron_sun=(alpha*I_sun+Jmax_sun-sqrt((alpha*I_sun+Jmax_sun)^2-4*theta*alpha*I_sun*Jmax_sun))/(2*theta);
        itera_err=1;
        A=20; % umol CO2/m2 ground/s, initial gross photosynsynthesis rate
        Rd_sun=v_Rd_leaf*12*10^6/Area_ground*LAI_sun/LAI; % umol CO2/m2 ground/s
        CO2i=CO2_a*0.7;
        gs_ini=g0+alpha_gs*A/((CO2i-gamma_star)*(1+VPD_leaf/VPD_0));% Leuning 1992: g0=0.01 mol/m2/s; VPD_0=0.35 kPa; alpha_gs=20; gm_leaf=0.3 mol/m2/s
        while itera_err>0.00001
            v_CO2_assimi_vc=Vcmax_sun*max(0,(CO2i-gamma_star))/(Km_CO2_assimi_leaf+CO2i);
            v_N_assimi_vc=coef_Pro_promote_leaf_photosyn ...
                *K_N_assimi_leaf_cat_sun*conc_iN_leaf/(Km_N_assimi_leaf+conc_iN_leaf);
            v_NADPH_produce=J_photosyn_electron_sun/2;
            v_NADPH_consume=(Vcmax_sun*CO2i/(Km_CO2_assimi_leaf+CO2i)+Vcmax_sun*gamma_star/(Km_CO2_assimi_leaf+CO2i))*2+v_N_assimi_vc*2.5; 
                % CBC, PR, (2NADPH per cycle) and N assimi. (NO3 -> AA needs 4.5NADPH+0.5ATP, NH4 -> AA needs 0.5NADPH+0.5ATP, so ave. N -> AA needs 2.5NADPH)
            discount_ratio=min(1,v_NADPH_produce*10^10./(v_NADPH_consume*10^10+0.00001));
            A=discount_ratio*v_CO2_assimi_vc*10^6; % umol CO2/m2 ground/s
            CO2i=CO2_a-(A-Rd_sun)/gs_ini;
            gs_now=g0+alpha_gs*A/((CO2i-gamma_star)*(1+VPD_leaf/VPD_0));
            itera_err=abs(gs_now-gs_ini);
            gs_ini=gs_now;
            CO2i=CO2_a-(A-Rd_sun)/gs_ini;
        end
        v_CO2_assimi_leaf_sun=Area_ground*discount_ratio*v_CO2_assimi_vc; % mol CO2/s
        v_N_assimi_leaf_sun=Area_ground*discount_ratio*v_N_assimi_vc; % mol N/s
        
        Jmax_shade=Jmax0_shade*coef_Pro_promote_leaf_photosyn*coef_NSC_inhibit_leaf_photosyn;
        Vcmax_shade=Vcmax0_shade*coef_Pro_promote_leaf_photosyn*coef_NSC_inhibit_leaf_photosyn;
        J_photosyn_electron_shade=(alpha*I_shade+Jmax_shade-sqrt((alpha*I_shade+Jmax_shade)^2-4*theta*alpha*I_shade*Jmax_shade))/(2*theta);
        itera_err=1;
        Rd_shade=v_Rd_leaf*12*10^6/Area_ground*LAI_shade/LAI; % umol CO2/m2 ground/s
        while itera_err>0.00001
            v_CO2_assimi_vc=Vcmax_shade*max(0,(CO2i-gamma_star))/(Km_CO2_assimi_leaf+CO2i); % mol CO2/m2 ground/s
            v_N_assimi_vc=coef_Pro_promote_leaf_photosyn ...
                *K_N_assimi_leaf_cat_shade*conc_iN_leaf/(Km_N_assimi_leaf+conc_iN_leaf);
            v_NADPH_produce=J_photosyn_electron_shade/2;
            v_NADPH_consume=(Vcmax_shade*CO2i/(Km_CO2_assimi_leaf+CO2i)+Vcmax_shade*gamma_star/(Km_CO2_assimi_leaf+CO2i))*2+v_N_assimi_vc*2.5; 
                % CBC, PR, (2NADPH per cycle) and N assimi. (NO3 -> AA needs 4.5NADPH+0.5ATP, NH4 -> AA needs 0.5NADPH+0.5ATP, so ave. N -> AA needs 2.5NADPH)
            discount_ratio=min(1,v_NADPH_produce*10^10./(v_NADPH_consume*10^10+0.00001));
            A=discount_ratio*v_CO2_assimi_vc*10^6; % umol CO2/m2 ground/s
            CO2i=CO2_a-(A-Rd_shade)/gs_ini;
            gs_now=g0+alpha_gs*A/((CO2i-gamma_star)*(1+VPD_leaf/VPD_0));
            itera_err=abs(gs_now-gs_ini);
            gs_ini=gs_now;
            CO2i=CO2_a-(A-Rd_shade)/gs_ini;
        end
        v_CO2_assimi_leaf_shade=Area_ground*discount_ratio*v_CO2_assimi_vc; % mol CO2/s
        v_N_assimi_leaf_shade=Area_ground*discount_ratio*v_N_assimi_vc; % mol N/s
            
        v_CO2_assimi_leaf=(v_CO2_assimi_leaf_sun+v_CO2_assimi_leaf_shade)/12; % convert to mol Suc/s, whole canopy
        v_N_assimi_leaf=v_N_assimi_leaf_sun+v_N_assimi_leaf_shade;
    else
        v_CO2_assimi_leaf=0;
        v_N_assimi_leaf=0;
        gs_now=g0;
        CO2i=CO2_a;
    end
    
    % daily leaf triose phosphate -> Suc and Star (or vac. Suc) storage
    v_TP_Suc_leaf=Area_leaf*K_TP_Suc_leaf_cat*max(0,conc_TP_leaf-conc_Suc_leaf/Ke_TP_Suc_leaf)/(conc_TP_leaf+Km_TP_Suc_leaf);
    v_TP_Star_leaf=Area_leaf*K_TP_Star_leaf_cat*conc_TP_leaf/(conc_TP_leaf+Km_TP_Star_leaf);
    
    % ODE in leaf
    retrieve_ratio_C=1;% retrieve ratio of C from senescence leaf
    maintain_cost=coef_Suc_maintain_leaf*Vol_leaf;
    if conc_Suc_leaf<0.001
        maintain_cost=0;
    end
    delta_conc_Suc_leaf=(v_TP_Suc_leaf+v_Star_Suc_leaf ... % CO2 assimi + Star degradation
        -coef_unload*(v_iN_leaf_unload+v_oN_leaf_unload)-coef_load*v_oN_leaf_load ...% unload and load cost
        -(1+coef_load)*v_Suc_leaf_load... % Suc load
        -retrieve_ratio_C*(conc_Suc_leaf+conc_Star_leaf)*delta_Vol_leaf... % senescence retrieve
        -coef_C_skeleton_iN2oN*v_N_assimi_leaf... % 1 mol iN N assimi needs 5 mol C skeleton
        -coef_oN_to_Pro*max(0,v_oN_Pro_leaf)-coef_Pro_decomp*(max(0,-v_oN_Pro_leaf) + conc_Pro_leaf*delta_Vol_leaf)... % AA <-> pro Rd
        -(coef_oN_to_Pro+coef_Pro_decomp)*totMol_Pro_leaf*K_Pro_turnover_leaf_cat... % Pro turnover Rd
        -coef_Suc_to_Star*(v_Star_Suc_leaf + conc_Star_leaf*delta_Vol_leaf) ... % Star decomposition
        -maintain_cost)/Vol_leaf; % maintain cost
    retrieve_ratio_N=1-min(conc_iN_oN_Retrieve_upper_leaf-conc_iN_oN_Retrieve_lower_leaf,max(0,conc_iN_oN_leaf-conc_iN_oN_Retrieve_lower_leaf))/100; 
            % retrieve ratio of N from senescence leaf decreased with higher N concentration
    delta_conc_iN_leaf=(v_iN_leaf_unload-retrieve_ratio_N*conc_iN_leaf*delta_Vol_leaf... % unload + retrieve
        -v_N_assimi_leaf)/Vol_leaf; % assimi
    delta_conc_oN_leaf=(v_oN_leaf_unload+v_N_assimi_leaf-retrieve_ratio_N*(conc_Pro_leaf+conc_oN_leaf)*delta_Vol_leaf... % unload, N assimi + retrieve
        -v_oN_leaf_load-v_oN_Pro_leaf)/Vol_leaf; % load, convert to Pro
    delta_photosyn_Pro_leaf=v_oN_Pro_leaf+conc_Pro_leaf*delta_Vol_leaf; % Pro <-> O-N and loss by senescence
    delta_Star_leaf=v_TP_Star_leaf-v_Star_Suc_leaf+conc_Star_leaf*delta_Vol_leaf;
    retrieve_ratio_TP=0;
    delta_conc_TP_leaf=(v_CO2_assimi_leaf-v_TP_Suc_leaf-v_TP_Star_leaf-retrieve_ratio_TP*conc_TP_leaf*delta_Vol_leaf)/Vol_leaf; % mol Suc/m3
    
    % ODE in leaf for record
    delta_CO2_assimi_leaf=v_CO2_assimi_leaf; % mol Suc/s
    delta_Suc_out_leaf=v_Suc_leaf_load; %mol Suc/s
    delta_oN_out_leaf=v_oN_leaf_load;
    delta_iN_in_leaf=v_iN_leaf_unload;
    delta_oN_in_leaf=v_oN_leaf_unload;
    delta_N_in_leaf=v_iN_leaf_unload+v_oN_leaf_unload-v_oN_leaf_load;
    delta_N_leak_leaf=-(1-retrieve_ratio_N)*(conc_Pro_leaf+conc_oN_leaf+conc_iN_leaf)*delta_Vol_leaf;
    delta_C_leak_leaf=-(1-retrieve_ratio_C)*(conc_Suc_leaf+conc_Star_leaf)*delta_Vol_leaf - (1-retrieve_ratio_TP)*conc_TP_leaf*delta_Vol_leaf;
    delta_Area_sene_leaf=v_leaf_sene;
    delta_Rd_leaf=v_Rd_leaf;
    delta_gs_leaf=gs_now;
    delta_CO2i=CO2i;
    delta_N_ass_leaf=v_N_assimi_leaf;
    %% grain metabolism
    % grow
    coef_glume_restrict=min(1,1/e_glume_restrict*(1-Vol_grain/Vol_grain_max));% physical restriction takes effect when current volume >70% total grain volume
    v_grain_grow=coef_glume_restrict*Area_grain*K_grain_grow_cat*conc_oN_grain/(Km_grain_grow_oN+conc_oN_grain) ...
        *conc_Suc_grain/(Km_grain_grow_Suc+conc_Suc_grain); %m3/s, single grain
    % Pro and Star store
    Vol_ear_solid=(store_Star_grain/pho_Star+store_Pro_grain/pho_Pro)*10^-6; % from g to m3, all grains
    Vol_ear_soluble=Vol_grain*grainNum-Vol_ear_solid; % all grains
    v_ear_Star_store=Vol_ear_soluble*K_grain_Star_store_cat*...
        conc_Suc_grain/(Km_grain_Star_store_Suc+conc_Suc_grain); %mol Suc/s, all grains
    v_ear_Pro_store=Vol_ear_soluble*K_grain_Pro_store_cat*...
        conc_oN_grain/(Km_grain_Pro_store_oN+conc_oN_grain); %mol N/s, all grains
    
    % unload
    v_oN_ear_unload=Area_grain*grainNum*K_oN_grain_unload_cat ...
        *max(0,(conc_oN_sphloem-conc_oN_grain/Ke_oN_grain_unload))/(Km_oN_grain_unload+conc_oN_sphloem); % all grains
    v_Suc_ear_unload=Area_grain*grainNum*K_Suc_grain_unload_cat ...
        *max(0,(conc_Suc_sphloem-conc_Suc_grain/Ke_Suc_grain_unload))/(Km_Suc_grain_unload+conc_Suc_sphloem); % all grains
    
    % ODE in grain
    delta_store_Pro_ear=v_ear_Pro_store*Pro_molWeight; %  convert mol N to g Pro
    delta_store_Star_ear=v_ear_Star_store*Suc_molWeight; % convert mol Suc to g Star
    delta_Vol_ear_solid=(delta_store_Pro_ear/pho_Pro+delta_store_Star_ear/pho_Star)*10^-6; % from g to m3, all grains
    delta_Vol_grain=v_grain_grow; % single grain
    delta_Area_grain=3.2240*((1+2*grain_len_wid_ratio^1.6075)/3)^0.6221*grain_len_wid_ratio^(-0.6667)*Vol_grain^(-1/3)*delta_Vol_grain; 
        %    given grain_len_wid_ratio=L; S=4pi*((1+2*L^1.6075)/3)^(1/1.6075)*r2, V=4/3*pi*r3*L
        % -> dS=4*pi*((1+2*L^1.6075)/3)^(1/1.6075)*(3/(4*pi*L))^(2/3)*2/3*V^(-1/3)*dV = 3.6376*V^(-1/3)*dV, single grain
    grow_Weight_ear=Fresh_new_dry_ratio_grain*v_grain_grow*grainNum/Vol_Fresh_ratio_grain*10^6; % g
    grow_Rd_ear=grow_Weight_ear*(1/grow_efficiency_grain-1)/Suc_molWeight; % mol C
    grow_Suc_cost_ear=grow_Rd_ear + grow_Weight_ear*(1-grain_struc_N*N_to_Pro_mass)/Suc_molWeight;
    maintain_cost=coef_Suc_maintain_grain*Vol_ear_soluble;
    if conc_Suc_grain<0.001
        maintain_cost=0;
    end
    delta_conc_Suc_grain=(v_Suc_ear_unload+conc_Suc_grain*delta_Vol_ear_solid... % Suc unload + concentration effect 
        -coef_unload*(v_oN_ear_unload+v_Suc_ear_unload)... % unload cost
        -conc_Suc_grain*v_grain_grow*grainNum-grow_Suc_cost_ear ... % grow cost
        -(coef_Suc_to_Star+1)*v_ear_Star_store-coef_oN_to_Pro*v_ear_Pro_store... % storage cost
        -maintain_cost)/Vol_ear_soluble; % maintain cost
    delta_conc_oN_grain=(v_oN_ear_unload+conc_oN_grain*delta_Vol_ear_solid... % unload + concentration effect 
        -v_ear_Pro_store... % storage
        -conc_oN_grain*v_grain_grow*grainNum ... % dilution effect due to new grown aqueous space
        -grow_Weight_ear*grain_struc_N/N_molWeight)/Vol_ear_soluble; % grow struc built cost

    % ODE in grain for record
    delta_Suc_in_ear=v_Suc_ear_unload;
    delta_oN_in_ear=v_oN_ear_unload;
    delta_Vol_grow_ear=v_grain_grow*grainNum;
    delta_Biomass_grow_ear=grow_Weight_ear + delta_store_Star_ear + delta_store_Pro_ear;
    delta_Rd_ear=coef_unload*(v_oN_ear_unload+v_Suc_ear_unload)... % unload Rd
        +grow_Rd_ear ... % grow Rd
        +coef_Suc_to_Star*v_ear_Star_store+coef_oN_to_Pro*v_ear_Pro_store... % storage Rd
        +maintain_cost; % maintain Rd. unit: mol Suc/s
    
    %% flow metabolism
    v_oN_xylem_to_lphloem=K_oN_xylem_to_lphloem_cat*max(0,conc_oN_xylem-conc_oN_lphloem/Ke_oN_xylem_to_lphloem)/ ...
        (Km_oN_xylem_to_lphloem+conc_oN_xylem);
    
    % flux between rphloem, sphloem and lphloem
    R_r_l_phloem=R_r_l_phloem0;
    R_l_s_phloem=R_l_s_phloem0;
    Flow_lphloem_rphloem=(conc_Suc_lphloem+conc_oN_lphloem-conc_Suc_rphloem-conc_oN_rphloem)/R_r_l_phloem; % m3/s, bi-directional
    Flow_lphloem_sphloem=(conc_Suc_lphloem+conc_oN_lphloem-conc_Suc_sphloem-conc_oN_sphloem)/R_l_s_phloem;
    % diffusion between stem and lphloem
    v_Suc_lphloem_stem=(conc_Suc_lphloem-conc_Suc_stem)/R_Suc_lphloem_stem;
    v_oN_lphloem_stem=(conc_oN_lphloem-conc_oN_stem)/R_oN_lphloem_stem;
    delta_conc_Suc_lphloem=-v_Suc_lphloem_stem/Vol_lphloem;
    delta_conc_oN_lphloem=(v_oN_xylem_to_lphloem-v_oN_lphloem_stem)/Vol_lphloem;
    % diffusison of oN or Suc from sphloem/rphloem to lphloem if sphloem/rphloem has higher concentration
    v_diffuse_oN_lphloem_sphloem=min(0,conc_oN_lphloem-conc_oN_sphloem)*Vol_sphloem;
    v_diffuse_oN_lphloem_rphloem=min(0,conc_oN_lphloem-conc_oN_rphloem)*Vol_rphloem;
    v_diffuse_Suc_lphloem_sphloem=min(0,conc_Suc_lphloem-conc_Suc_sphloem)*Vol_sphloem;
    v_diffuse_Suc_lphloem_rphloem=min(0,conc_Suc_lphloem-conc_Suc_rphloem)*Vol_rphloem;

    % ODE in flow
    delta_conc_iN_xylem=(v_iN_root_load-v_iN_leaf_unload)/Vol_xylem;
    delta_conc_oN_xylem=(v_oN_root_load-v_oN_leaf_unload-v_oN_xylem_to_lphloem)/Vol_xylem;

    if Flow_lphloem_rphloem>0
        delta_conc_Suc_rphloem=(Flow_lphloem_rphloem*conc_Suc_lphloem-v_Suc_root_unload)/Vol_rphloem;
        delta_conc_oN_rphloem=(Flow_lphloem_rphloem*conc_oN_lphloem-v_oN_root_unload)/Vol_rphloem;
        delta_conc_Suc_lphloem=delta_conc_Suc_lphloem+(-Flow_lphloem_rphloem*conc_Suc_lphloem+v_Suc_leaf_load)/Vol_lphloem;
        delta_conc_oN_lphloem=delta_conc_oN_lphloem+(-Flow_lphloem_rphloem*conc_oN_lphloem+v_oN_leaf_load)/Vol_lphloem;
    else
        delta_conc_Suc_rphloem=(Flow_lphloem_rphloem*conc_Suc_rphloem-v_Suc_root_unload)/Vol_rphloem;
        delta_conc_oN_rphloem=(Flow_lphloem_rphloem*conc_oN_rphloem-v_oN_root_unload)/Vol_rphloem;
        delta_conc_Suc_lphloem=delta_conc_Suc_lphloem+(-Flow_lphloem_rphloem*conc_Suc_rphloem+v_Suc_leaf_load)/Vol_lphloem;
        delta_conc_oN_lphloem=delta_conc_oN_lphloem+(-Flow_lphloem_rphloem*conc_oN_rphloem+v_oN_leaf_load)/Vol_lphloem;
    end
    if Flow_lphloem_sphloem>0
        delta_conc_Suc_sphloem=(Flow_lphloem_sphloem*conc_Suc_lphloem-v_Suc_ear_unload)/Vol_sphloem;
        delta_conc_oN_sphloem=(Flow_lphloem_sphloem*conc_oN_lphloem-v_oN_ear_unload)/Vol_sphloem;
        delta_conc_Suc_lphloem=delta_conc_Suc_lphloem-(Flow_lphloem_sphloem*conc_Suc_lphloem)/Vol_lphloem;
        delta_conc_oN_lphloem=delta_conc_oN_lphloem-(Flow_lphloem_sphloem*conc_oN_lphloem)/Vol_lphloem;
    else
        delta_conc_Suc_sphloem=(Flow_lphloem_sphloem*conc_Suc_sphloem-v_Suc_ear_unload)/Vol_sphloem;
        delta_conc_oN_sphloem=(Flow_lphloem_sphloem*conc_oN_sphloem-v_oN_ear_unload)/Vol_sphloem;
        delta_conc_Suc_lphloem=delta_conc_Suc_lphloem-(Flow_lphloem_sphloem*conc_Suc_sphloem)/Vol_lphloem;
        delta_conc_oN_lphloem=delta_conc_oN_lphloem-(Flow_lphloem_sphloem*conc_oN_sphloem)/Vol_lphloem;
    end
    delta_conc_Suc_sphloem=delta_conc_Suc_sphloem + v_diffuse_Suc_lphloem_sphloem/Vol_sphloem;
    delta_conc_oN_sphloem=delta_conc_oN_sphloem + v_diffuse_oN_lphloem_sphloem/Vol_sphloem;
    delta_conc_Suc_lphloem=delta_conc_Suc_lphloem - (v_diffuse_Suc_lphloem_sphloem+v_diffuse_Suc_lphloem_rphloem)/Vol_lphloem;
    delta_conc_oN_lphloem=delta_conc_oN_lphloem - (v_diffuse_oN_lphloem_sphloem+v_diffuse_oN_lphloem_rphloem)/Vol_lphloem;
    delta_conc_Suc_rphloem=delta_conc_Suc_rphloem + v_diffuse_Suc_lphloem_rphloem/Vol_rphloem;
    delta_conc_oN_rphloem=delta_conc_oN_rphloem + v_diffuse_oN_lphloem_rphloem/Vol_rphloem;
    % ODE for flow record
    delta_Flow_lphloem_rphloem=Flow_lphloem_rphloem;
    delta_Flow_lphloem_sphloem=Flow_lphloem_sphloem;
    delta_oN_xylem_to_lphloem=v_oN_xylem_to_lphloem;
    
    %% stem metabolism
    conc_star_CSP=Star_stem/Suc_molWeight/Vol_stem;
    conc_pro_CSP=Pro_stem/Pro_molWeight/Vol_stem;
    
    if conc_Suc_stem>conc_lower_Suc_to_stem_Star
        v_Suc_Star_stem=Vol_soluble_stem*K_Suc_Star_stem_cat ...
            *(conc_Suc_stem-conc_lower_Suc_to_stem_Star)/(Km_Suc_Star_stem+conc_Suc_stem); 
    else
        v_Suc_Star_stem=Vol_soluble_stem*K_Star_Suc_stem_cat ...
            *(conc_Suc_stem-conc_lower_Suc_to_stem_Star)/(conc_lower_Suc_to_stem_Star+conc_Suc_stem) ...
            *conc_star_CSP/(Km2_Star_Suc_stem+conc_star_CSP); 
    end % mol Suc/s, Star decomposition and synthesis
    
    if conc_oN_stem>conc_lower_oN_to_stem_Pro
        v_oN_Pro_stem=Vol_soluble_stem*K_oN_Pro_stem_cat ...
            *(conc_oN_stem-conc_lower_oN_to_stem_Pro)/(Km_oN_Pro_stem+conc_oN_stem); 
    else
        v_oN_Pro_stem=Vol_soluble_stem*K_Pro_oN_stem_cat ...
            *(conc_oN_stem-conc_lower_oN_to_stem_Pro)/(Km_oN_Pro_stem+conc_oN_stem)...
            *conc_pro_CSP/(Km2_Pro_oN_stem+conc_pro_CSP); 
    end % mol N/s, Pro decomposition and synthesis
    
    % ODE for stem
    delta_Star_stem=v_Suc_Star_stem*Suc_molWeight; % g/s
    delta_Pro_stem=v_oN_Pro_stem*Pro_molWeight; % g/s
    delta_Vol_soluble_stem=-delta_Star_stem*10^-6/pho_Star - delta_Pro_stem*10^-6/pho_Pro;
    delta_conc_Suc_stem=(v_Suc_lphloem_stem... % diffusion 
        -conc_Suc_stem*delta_Vol_soluble_stem... % concentration/dilution effect
        -v_Suc_Star_stem... % convert to Star
        -coef_Suc_maintain_stem*Vol_soluble_stem... % maintain cost
        -coef_Suc_to_Star*abs(v_Suc_Star_stem)... % Suc<->Star inter-conversion cost
        -coef_oN_to_Pro*max(0,v_oN_Pro_stem)-coef_Pro_decomp*max(0,-v_oN_Pro_stem))/Vol_soluble_stem; % AA<->Pro inter-conversion cost
    delta_conc_oN_stem=(v_oN_lphloem_stem... % diffusion
        -conc_oN_stem*delta_Vol_soluble_stem... % concentration/dilution effect
        -v_oN_Pro_stem)/Vol_soluble_stem; % to Pro
    % ODE for stem record
    delta_Rd_stem=coef_Suc_maintain_stem*Vol_soluble_stem+coef_Suc_to_Star*abs(v_Suc_Star_stem)+...
        coef_oN_to_Pro*max(0,v_oN_Pro_stem)+coef_Pro_decomp*max(0,-v_oN_Pro_stem);
    
    mc = [delta_Vol_root,delta_Dry_root,delta_conc_Suc_root,delta_conc_iN_root,delta_conc_oN_root,delta_conc_iN_soil,...
            delta_Vol_soil_water,zeros(1,100-7),...
          delta_Area_leaf,delta_Vol_leaf,delta_conc_Suc_leaf,delta_conc_iN_leaf,delta_conc_oN_leaf,delta_photosyn_Pro_leaf,...
            delta_Star_leaf,delta_conc_TP_leaf,zeros(1,100-8),...
          delta_Vol_grain,delta_Area_grain,delta_conc_Suc_grain,delta_conc_oN_grain,delta_store_Pro_ear,delta_store_Star_ear,...
            zeros(1,100-6),...
          delta_conc_iN_xylem,delta_conc_Suc_rphloem,delta_conc_Suc_lphloem,delta_conc_Suc_sphloem,delta_conc_oN_rphloem,...
            delta_conc_oN_lphloem,delta_conc_oN_sphloem,delta_conc_oN_xylem,zeros(1,100-8),...
          delta_conc_Suc_stem,delta_conc_oN_stem,delta_Star_stem,delta_Pro_stem,delta_Vol_soluble_stem,zeros(1,100-5),...
          delta_Suc_in_root,delta_oN_shoot_root_netFlux,delta_iN_out_root,delta_Vol_sene_root,delta_Vol_grow_root,delta_Rd_root,...
            delta_N_assimi_root,delta_N_absorb_root,delta_N_shoot_root_netFlux,delta_oN_import_root,delta_oN_export_root,delta_iN_export_root,zeros(1,20-12),...
            delta_CO2_assimi_leaf,delta_Suc_out_leaf,delta_oN_out_leaf,delta_iN_in_leaf,delta_oN_in_leaf,...
            delta_N_in_leaf,delta_N_ass_leaf,delta_N_leak_leaf,delta_C_leak_leaf,delta_Area_sene_leaf,delta_Rd_leaf,delta_gs_leaf,...
            delta_CO2i,zeros(1,20-13),...
            delta_Suc_in_ear,delta_oN_in_ear,delta_Vol_grow_ear,delta_Biomass_grow_ear,delta_Rd_ear,zeros(1,20-5),...
            delta_Flow_lphloem_rphloem,delta_Flow_lphloem_sphloem,delta_oN_xylem_to_lphloem,zeros(1,20-3),...
            delta_Rd_stem,zeros(1,20-1)]';
    