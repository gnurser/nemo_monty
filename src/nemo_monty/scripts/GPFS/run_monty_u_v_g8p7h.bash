#!/usr/bin/env bash
export R12=/noc/msm/working/poles_apart/GO8p7/eORCA12
time 3dfiles_M2_mp.py --meshdir $R12/../../stefryn --meshfile mesh_mask_CLASS-MEDUSA.nc --infile $R12/monthly/nemo_g8p7ho_1m_2021??-2021??_grid-T.nc  --tracers ssh T S u v --xtracers z_s sigma_s sigma_med_s T_s S_s outcrop_s incrop_s mont  --density 26.5 26.8 27.0 --iterate_TS0 all --deltaTS 0.01 0.004 --TS0 2. 34.5 --nthreads 16 --neos 2 --dims t y x --ylimits 2782 3605 --xlimits 1 3823 --passive_s u_s v_s
