#!/usr/bin/env bash
export D025=/noc/msm/scratch/nemo2/agn/NUCLEAR/NPD_Domains/025
time 3dfiles_M2_mp.py --meshdir $D025 --meshfile mesh_mask.nc --infile eORCA025_1m_grid_T_200506-200506.nc \
--tracers ssh T S age u v --xtracers z_s sigma_s sigma_med_s T_s S_s outcrop_s incrop_s age_s mont \
--passive_s u_s v_s age_s --density 26.5 26.8 27.0 --iterate_TS0 all --deltaTS 0.01 0.004 --TS0 2. 34.5 \
     --nthreads 16 --neos 2 --dims t y x --ylimits 700 1205 --xlimits 1 1439
    # --nthreads 16 --neos 2 --dims t y x --ylimits 700 1206 --xlimits 1 1440 \
