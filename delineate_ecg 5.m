function [QRS_start_idx, QRS_end_idx, QRS_end_flag, QRS_duration, t_wave_duration, QT_duration, ...
    t_peak_end, t_peak_idx, t_start_peak, t_wave_start_idx, t_wave_end_idx, ...
    inverse, t_magnitude_true] = delineate_ecg(ECG_table, lead, width, ...
    activationTime)
    % Calculates and extracts the following biomarkers:
    % - QRS Start Point, QRS End Point, QRS warning Flag, QRS Duration,
    % T wave duration, QT interval
    % - T peak to end duration, T Peak Point, T start to end duration, 
    % T Start Point, T End Point, T wave inversion flag, T magnitude

    % inputs - ECG table with columns for time, and voltage for each lead
    % width of window for sliding window, 3 is recommended
    % maximum activation time for the simulation - can put 0 if not known
    % though the QRS end flag will probably be true
    
    %% 0) Setup

    T = ECG_table.Time;
    T = T - T(1);
    V = ECG_table{:,lead};
    dV = gradient(V);
    ddV = gradient(dV);

    T_ex = [(-3:1)'; T];
    V_ex = [repmat(0, 3, 1); V];
    dV_ex = gradient(V_ex);
    ddV_ex = gradient(dV_ex);

    for i = width+1:300 
        dV_windowed(i-width) = abs(get_window(dV_ex, i, width));
    end

    % Tolerance values
    QRS_start_tol = 0.01*max(abs(V))/30;
    QRS_end_tol_ddV = 0.1*max(abs(V))/(30*width+2);
    QRS_end_tol_dV = 0.07*max(abs(dV_windowed));
    T_start_tol = 0.12*max(abs(V))/30;
    T_end_tol = 0.01*max(abs(V))/30;
    
    %% 1) Determine QRS start time 


    for i = width+1:100
        QRS_window = abs(get_window(ddV_ex, i, width));
        if QRS_window > QRS_start_tol
            QRS_start_idx = i-width;
            break
        end
    end 
    

    %% 2) Determing QRS end time, and QRS duration

    for i = 30:500
        QRS_window2(i-30+1) = get_window(abs(ddV_ex), i, 2);
    end


    QRS_end_idx_lst = 30 + find(QRS_window2<QRS_end_tol_ddV)-2;
    
    for idx = QRS_end_idx_lst
        dV_window = get_window(abs(dV), idx, width+2);
        if abs(dV_window) < QRS_end_tol_dV
            QRS_end_idx = idx-(width-2);
            break
        end
    end

    QRS_start_time = T(QRS_start_idx);
    QRS_end_time = T(QRS_end_idx);
    QRS_duration = QRS_end_time - QRS_start_time;
    
    %% 3) Check QRS against activation times, and flag simulation if this is
    %    is not compatible.
    QRS_end_flag = abs(QRS_end_time-activationTime)>=25;

    %% 4) Find T-wave Peak and T-wave peak time

    segment = V(QRS_end_idx+100:999);
    t_magnitude = max(abs(segment));
    peak_idx = find(abs(segment)==abs(t_magnitude), 1, "last");
    t_peak_idx = QRS_end_idx + 100 + peak_idx - 1;
    t_sign = sign(V(t_peak_idx));
    t_magnitude_true = t_sign * t_magnitude;

    t_wave_peak_time = T(t_peak_idx);

    %% 5) Find T-wave end point
    
    for i = length(V)-width:-1:t_peak_idx
        tend_window = abs(get_window(dV, i, width));
        if tend_window > T_end_tol
            t_wave_end_idx = i;
            break
        end
    end

    t_wave_end_time = T(t_wave_end_idx);
    
    %% 6) Find T-wave start point and calculate QT
    try
        for i = t_peak_idx-20:-1:QRS_end_idx
            tstart_window = abs(get_window(dV, i, width));
            if tstart_window < T_start_tol
                i
                tstart_window
                if max(abs(V(i-30:i-1))) < abs(V(i)) & abs(V(i)) < 0.5 * t_magnitude
                    t_wave_start_idx = i;
                    break
                end
            end
        end
        t_wave_start_time = T(t_wave_start_idx);
    catch
        t_wave_start_idx = nan;
        t_wave_start_time = nan;
    end
    
        
    
    t_wave_duration = t_wave_end_time - t_wave_start_time;
    QT_duration = t_wave_end_time - QRS_start_time;
    t_peak_end = t_wave_end_time - t_wave_peak_time;
    t_start_peak = t_wave_peak_time - t_wave_start_time;

    segment = V(t_wave_end_idx-10:t_wave_end_idx);
    if max(abs(segment)) > t_magnitude * 0.1
        if max(segment) > V(t_wave_end_idx)
            inverse = 1;
        else
            inverse = 0;
        end
    else
        inverse = 0;
    end

end
