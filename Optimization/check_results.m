% This script investigates the results and creates stacked plots

frac_RE_wec= 0.25;
RE_penetration= 0.75; %0.25;

P_Capacity_Storage_frac= 0.5;
sol_50= run_iterate_storage_and_wave(frac_RE_wec, P_Capacity_Storage_frac, RE_penetration);

P_Capacity_Storage_frac= 1;
sol_100= run_iterate_storage_and_wave(frac_RE_wec, P_Capacity_Storage_frac, RE_penetration);