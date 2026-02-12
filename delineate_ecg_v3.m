function [QRS_start_idx, QRS_end_idx, QRS_end_flag, QRS_duration, t_wave_duration, QT_duration, ...
    t_peak_end, t_peak_idx, t_start_peak, t_wave_start_idx, t_wave_end_idx, ...
    inverse, t_magnitude_true, diag] = delineate_ecg_v3(ECG_table, lead, width, activationTime, duration)

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
    T = ECG_table.t;
    T = T - T(1);
    V = ECG_table{:, lead};

    n = numel(V);
    duration = min(duration, n); % duration treated as index (dt=1ms)
    if duration < 200
        duration = n;
    end

    dV  = gradient(V);
    ddV = gradient(dV);

    V_ex = [repmat(V(1), 3, 1); V];
    dV_ex = gradient(V_ex);
    ddV_ex = gradient(dV_ex);

    %% Estimate baseline noise (sigma)
    n = numel(V);
    
    % Baseline indices: early + late
    baseA = 1:min(80,n);
    baseB = max(n-150,1):max(n-20,1);
    baseIdx = unique([baseA, baseB]);
    
    sigma_dV  = 1.4826 * mad(dV(baseIdx),  1);
    sigma_ddV = 1.4826 * mad(ddV(baseIdx), 1);
    
    % If MAD collapses, keep it nonzero
    sigma_dV  = max(sigma_dV,  eps);
    sigma_ddV = max(sigma_ddV, eps);
    
    sigma_dV  = max(sigma_dV,  0.01 * max(abs(dV))  + eps);   % 1% of peak slope
    sigma_ddV = max(sigma_ddV, 0.01 * max(abs(ddV)) + eps);   % 1% of peak curvature

    % Thresholds
    k_qrs_start   = 0.5;
    k_qrs_end     = 1;
    k_qrs_end_dV  = 1;
    k_t_start     = 1;
    k_t_end       = 0.6;
    
    QRS_start_tol      = k_qrs_start  * sigma_ddV;
    QRS_end_tol_ddV    = k_qrs_end    * sigma_ddV;
    T_start_tol        = k_t_start    * sigma_dV;
    T_end_tol          = k_t_end      * sigma_dV;
    QRS_end_tol_dV     = k_qrs_end_dV * sigma_dV;   % start here; tune 6–12

    % Diagnosis struct to extract information for tuning.
    diag = struct();
    diag.sigma_dV  = sigma_dV;
    diag.sigma_ddV = sigma_ddV;
    diag.thresholds = struct( ...
        'QRS_start_tol',   QRS_start_tol, ...
        'QRS_end_tol_ddV', QRS_end_tol_ddV, ...
        'QRS_end_tol_dV',  QRS_end_tol_dV, ...
        'T_start_tol',     T_start_tol, ...
        'T_end_tol',       T_end_tol ...
    );
    diag.meta = struct( ...
        'lead', lead, ...
        'beat_length', numel(V), ...
        'activationTime', activationTime ...
    );

    %% 1) Determine QRS start index
    QRS_start_idx = NaN;
    qrsStartSearchStart = max(width+1, 15);
    qrsStartSearchEnd   = min(n, 200);           

    % Add a gentle dV gate to prevent V1 early false triggers
    k_qrs_start_dV = 0.2;                     % start here
    QRS_start_tol_dV = k_qrs_start_dV * sigma_dV;
        
    for i = qrsStartSearchStart:qrsStartSearchEnd
        w_ddV = get_window(abs(ddV_ex), i, width);
        w_dV  = get_window(abs(dV_ex),  i, width);
    
        if (w_ddV > QRS_start_tol) && (w_dV > QRS_start_tol_dV)
            QRS_start_idx = i - width;
    
            % Refinement: walk back to earliest sustained slope rise
            refBackMax = 30;   % ms to look back
            refSustain = 6;    % ms sustained
            refFrac    = 0.15; % fraction of local peak slope
    
            j0 = max(1, QRS_start_idx - refBackMax);
            j1 = min(n, QRS_start_idx + 20);
            local = abs(dV(j0:j1));
            localPeak = max(local);
    
            if localPeak > 0
                thr = refFrac * localPeak;
                for jj = j0:(QRS_start_idx - refSustain)
                    if all(abs(dV(jj:jj+refSustain)) > thr)
                        QRS_start_idx = jj;
                        break
                    end
                end
            end
    
            break
        end
    end

    % Fallback: choose strongest ddV event in early window if threshold never crossed
    if isnan(QRS_start_idx)
        win = qrsStartSearchStart:qrsStartSearchEnd;
        [~, rel] = max(abs(ddV(win)));
        imax = win(1) + rel - 1;
        QRS_start_idx = max(1, imax - 20);
    end


    %% 2) Determine QRS end index
    QRS_end_idx = NaN;
    
    % Search window after QRS start
    qrsEndSearchStart = max(1, QRS_start_idx + 20);
    qrsEndSearchEnd   = min(n, QRS_start_idx + 160);
    
    holdMs = 8;   % require ddV to remain low for this many ms
    
    % Scan forward for first sustained ddV "quiet" region
    for idx = qrsEndSearchStart : (qrsEndSearchEnd - holdMs)
        w = abs(ddV(idx:idx+holdMs));
        if max(w) < QRS_end_tol_ddV
            QRS_end_idx = idx;
            break
        end
    end
    
    % Fallback: choose point of minimum |ddV| in the search window
    if isnan(QRS_end_idx)
        win = qrsEndSearchStart:qrsEndSearchEnd;
        [~, rel] = min(abs(ddV(win)));
        QRS_end_idx = win(1) + rel - 1;
    end
    
    QRS_end_idx = max(1, min(n, round(QRS_end_idx)));
    
    QRS_start_time = T(QRS_start_idx);
    QRS_end_time   = T(QRS_end_idx);
    QRS_duration   = QRS_end_time - QRS_start_time;

    %% 3) Activation time compatibility flag
    QRS_end_flag = abs(QRS_end_time - activationTime) >= 25;

    %% 4) Find T-wave peak
    tPeakStart = QRS_end_idx + 40;
    tPeakEnd   = duration;

    if tPeakEnd <= tPeakStart + 5
        tPeakStart = min(n, QRS_end_idx + 20);
        tPeakEnd   = n;
    end

    segment = V(tPeakStart:tPeakEnd);
    [t_magnitude, relIdx] = max(abs(segment));
    t_peak_idx = tPeakStart + relIdx - 1;

    t_sign = sign(V(t_peak_idx));
    if t_sign == 0
        t_sign = 1;
    end
    t_magnitude_true = t_sign * t_magnitude;

    t_wave_peak_time = T(t_peak_idx);

    %% 5) Find T-wave end
    t_wave_end_idx = NaN;
    
    sustain  = 35;          % ms
    maxDelay = 260;         % ms
    minDelay = 30;          % ms
    
    % Baseline for T-end
    tailStart = min(n, t_peak_idx + 120);
    tailEnd   = min(n, t_peak_idx + maxDelay);
    
    if tailEnd > tailStart + 10
        baseLevel = median(V(tailStart:tailEnd));
    else
        % fallback to ST-based baseline if beat is short
        b1 = min(n, QRS_end_idx + 20);
        b2 = min(n, QRS_end_idx + 60);
        if b2 <= b1
            b1 = min(n, QRS_end_idx + 5);
            b2 = min(n, QRS_end_idx + 30);
        end
        baseLevel = median(V(b1:b2));
    end
    
    Apeak = abs(V(t_peak_idx) - baseLevel);
    Apeak = max(Apeak, 1e-12);
    
    % Amplitude tolerance
    fracEnd = 0.15;
    ampTolEnd = fracEnd * Apeak;
    
    ampTolFloor = 0.01 * max(abs(V));     % floor for low-amplitude tails
    ampTolEnd   = max(ampTolEnd, ampTolFloor);
    
    iStart = min(n, t_peak_idx + minDelay);
    iStop  = min(n - sustain - 1, t_peak_idx + maxDelay);
    
    V_s   = movmean(V, 7);
    dV_s  = gradient(V_s);
    
    if iStop > iStart
        for i = iStart:iStop
            flatSlope = max(abs(dV_s(i:i+sustain))) < T_end_tol;
            smallAmp  = max(abs(V(i:i+sustain) - baseLevel)) < ampTolEnd;
    
            if flatSlope && smallAmp
                t_wave_end_idx = i;
                break
            end
        end
    end
    
    if isnan(t_wave_end_idx) && (iStop > iStart)
        for i = iStart:iStop
            if max(abs(dV_s(i:i+sustain))) < T_end_tol
                t_wave_end_idx = i;
                break
            end
        end
    end
    
    if isnan(t_wave_end_idx)
        t_wave_end_idx = min(n, t_peak_idx + 140);
    end
    
    % Clamp
    t_wave_end_idx = max(1, min(n, round(t_wave_end_idx)));
    t_wave_end_time = T(t_wave_end_idx);
    
    % Fallback: if T end pins to end of beat, re-estimate using relaxed smoothed slope ----
    pinTol = 2;
    if t_wave_end_idx >= (n - pinTol)
    
        V_s2  = movmean(V, 11);
        dV_s2 = gradient(V_s2);
    
        sustain2  = 25;
        iStart2   = min(n, t_peak_idx + 30);
        iStop2    = min(n - sustain2 - 1, t_peak_idx + 260);
    
        t2 = NaN;
        if iStop2 > iStart2
            for i = iStart2:iStop2
                if max(abs(dV_s2(i:i+sustain2))) < (1.5 * T_end_tol)
                    t2 = i;
                    break
                end
            end
        end
    
        if isfinite(t2)
            t_wave_end_idx = max(1, min(n, round(t2)));
            t_wave_end_time = T(t_wave_end_idx);
        end
    end
    
    % Fallback: post-end bump check
    % If there's a meaningful secondary bump after current end, extend end past it.
    postWinStart = min(n, t_wave_end_idx + 5);
    postWinEnd   = min(n, t_peak_idx + maxDelay);
    
    if postWinEnd > postWinStart + 10
    
        postSeg = V(postWinStart:postWinEnd);
    

        bumpAmp = max(abs(postSeg - median(postSeg)));

        bumpFrac = 0.12;
        if bumpAmp > bumpFrac * Apeak
    
            [~, rel] = max(abs(postSeg - median(postSeg)));
            bumpIdx = postWinStart + rel - 1;
    
            sustain3 = 20;
            searchStart = min(n, bumpIdx + 10);
            searchEnd   = min(n - sustain3 - 1, bumpIdx + 180);
    
            for i = searchStart:searchEnd
                if max(abs(dV_s(i:i+sustain3))) < (1.5 * T_end_tol)
                    t_wave_end_idx = i;
                    t_wave_end_time = T(t_wave_end_idx);
                    break
                end
            end
        end
    end


    %% 6) Find T-wave start
    t_wave_start_idx = NaN;
    
    % Define an ST plateau reference shortly after QRS end
    st1 = min(n, QRS_end_idx + 15);
    st2 = min(n, QRS_end_idx + 60);
    if st2 <= st1
        st1 = min(n, QRS_end_idx + 5);
        st2 = min(n, QRS_end_idx + 30);
    end
    
    stLevel = median(V(st1:st2));
    
    % T amplitude relative to ST plateau
    Apeak = abs(V(t_peak_idx) - stLevel);
    Apeak = max(Apeak, 1e-12);
    
    % Threshold: require departure of X% of peak, sustained for Y ms
    fracStart  = 0.16;
    sustainST  = 12;
    
    ampTolStart = fracStart * Apeak;
    
    % Polarity-aware slope direction: +1 for upright T, -1 for inverted T
    dV = gradient(V);
    tSign = sign(V(t_peak_idx) - stLevel);
    if tSign == 0, tSign = 1; end
    
    iStart = min(n, QRS_end_idx + 10);
    iStop  = min(n - sustainST - 1, t_peak_idx - 10);
    
    if iStop > iStart
        for i = iStart:iStop
            ampOK   = all(abs(V(i:i+sustainST) - stLevel) > ampTolStart);
            slopeOK = mean(tSign * dV(i:i+sustainST)) > 0;   % slope toward T peak
    
            if ampOK && slopeOK
                t_wave_start_idx = i;
                break
            end
        end
    end
    
    % Fallback: if not found, place it a conservative distance before peak
    if isnan(t_wave_start_idx)
        t_wave_start_idx = max(1, t_peak_idx - 80);
    end
    
    t_wave_start_idx = max(1, min(n, round(t_wave_start_idx)));
    t_wave_start_time = T(t_wave_start_idx);



    %% Durations
    t_wave_duration = t_wave_end_time - t_wave_start_time;
    QT_duration     = t_wave_end_time - QRS_start_time;
    t_peak_end      = t_wave_end_time - t_wave_peak_time;
    t_start_peak    = t_wave_peak_time - t_wave_start_time;

    %% Inversion flag
    inverse = 0;
    w0 = max(1, t_wave_end_idx-10);
    segment2 = V(w0:t_wave_end_idx);

    if max(abs(segment2)) > t_magnitude * 0.1
        if max(segment2) > V(t_wave_end_idx)
            inverse = 1;
        else
            inverse = 0;
        end
    end
end
