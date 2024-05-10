%example of running Hannah and Maxx's ECG delineation code

lead = 'V1';
width = 3;
load('all_1_table.mat')
activationTime = 67;


[QRS_start_idx, QRS_end_idx, QRS_end_flag, QRS_duration, t_wave_duration, QT_duration, ...
    t_peak_end, t_peak_idx, t_start_peak, t_wave_start_idx, t_wave_end_idx, ...
    inverse, t_magnitude_true] = delineate_ecg(ECG_table, lead, width, ...
    activationTime)
