%% ECG Delineation Script - MonoAlg3D / Pfizer Project Version

close all
clear

%% Setup
transUnits   = 1;      % scale units to ms if needed
cycleLength  = 1000;   % assuming 1 KHz sampling
analyseBeat  = 3;   % First beat you wish to analyse
analyseBeatEnd = analyseBeat; % Last beat you wish to analyse

leads        = ["V1","V2","V3","V4","V5","V6"];
leadNames    = ["V1","V2","V3","V4","V5","V6"];
width         = 3;     % ms window
activationTime= 67;    % ms until activation of QRS
duration      = cycleLength;   % ms of data examined by delineation

primaryQTLeads = ["V5", "V6"];
QT_sanity = struct( ... 
    'minQT',    150, ... % ms
    'maxQT',    750, ... % ms
    'minQRS',   20,  ... % ms
    'maxQRS',   200 ... % ms
    );

parentDir = uigetdir(pwd, 'Select Parent Directory');
if parentDir == 0
    error('No directory selected.');
end

resultsDir = uigetdir(pwd, 'Select Results Save Directory');
if resultsDir == 0
    error('No results directory selected.');
end

% --- dictionary-like container: simName -> simResult struct
allResults = containers.Map('KeyType','char','ValueType','any');

% Obtain list of subdirectories (exclude '.' and '..')
d = dir(parentDir);
isSub = [d.isdir] & ~ismember({d.name}, {'.', '..'});
subdirs = d(isSub);

%% Loop through subdirectories
for k = 1:numel(subdirs)
    thisDir = fullfile(parentDir, subdirs(k).name);
    simName = string(subdirs(k).name);

    fname = fullfile(thisDir, 'ecg.txt');
    if exist(fname, 'file') ~= 2
        continue
    end

    % Load ecg.txt robustly
    try
        ecgData = load(fname, '-ascii');
    catch
        ecgData = readmatrix(fname);
    end

    if isempty(ecgData) || size(ecgData,2) < 11
        warning('Skipping %s: ecg.txt empty or has < 11 columns.', simName);
        continue
    end

    n = size(ecgData, 1);

    % Build processed ECG table
    ECG_table = buildECGTable(ecgData, transUnits);

    % Split into beats

    minRequired = cycleLength/2;

    if n < minRequired
        warning('Skipping %s: only %d ms available (need ≥%d ms).', simName, n, minRequired);
        continue
    end

    % Always analyse at least one beat; last beat may be partial
    completeBeats = floor(n / cycleLength);
    if completeBeats < 1
        totalBeats = 1;
    else
        % include an additional partial beat if remainder exists
        totalBeats = completeBeats;
        remainder = n - completeBeats * cycleLength;
        if remainder > minRequired; totalBeats = completeBeats + (remainder > 0); end
    end
    % ------------------------------

    simResult = struct();
    simResult.simName      = simName;
    simResult.sourceDir    = thisDir;
    simResult.ecgFile      = fname;
    simResult.cycleLength  = cycleLength;
    simResult.nSamples     = n;
    simResult.nBeats       = totalBeats;
    simResult.ECG_table    = ECG_table;
    simResult.beats        = repmat(struct(), 1, totalBeats);

    for j = analyseBeat:analyseBeatEnd
        s = (j-1)*cycleLength + 1;
        e = min(j*cycleLength, n);
        isPartial = (e - s + 1) < cycleLength;

        if s > n
            break
        end

        beatTbl = ECG_table(s:e, :);
        beatTbl.t = beatTbl.t - beatTbl.t(1);  % reset time to 0

        % Create a new results object per beat
        resObj = ECGResults(char(simName), {}); 

        beatAnalysis = struct();   % per-lead analysis stored here

        % Plot figure for this beat
        fig = figure(1);
        clf(fig);
        set(fig, 'Name', sprintf('%s - Beat %d', simName, j), 'NumberTitle', 'off');

        primaryQT = NaN;
        primaryQTLead = "";
        primaryQTok = false;

        for i = 1:numel(leads)
            lead = leads(i);

            time = beatTbl.t;
            V    = beatTbl.(lead);

            duration_eff = min(duration, numel(time));

            [QRS_start_idx, QRS_end_idx, QRS_end_flag, QRS_duration, ...
             t_wave_duration, QT_duration, t_peak_end, t_peak_idx, ...
             t_start_peak, t_wave_start_idx, t_wave_end_idx, inverse, ...
             t_magnitude_true, diag] = delineate_ecg_v3(beatTbl, lead, width, activationTime, duration_eff);

            % Clamp indices defensively (in case delineator returns out-of-range)
            QRS_start_idx    = clampIndex(QRS_start_idx, numel(time));
            QRS_end_idx      = clampIndex(QRS_end_idx,   numel(time));
            t_wave_start_idx = clampIndex(t_wave_start_idx, numel(time));
            t_peak_idx       = clampIndex(t_peak_idx, numel(time));
            t_wave_end_idx   = clampIndex(t_wave_end_idx, numel(time));

            entry = struct( ...
                'lead', lead, ...
                'QRS_start_idx', QRS_start_idx, ...
                'QRS_end_idx', QRS_end_idx, ...
                'QRS_end_flag', QRS_end_flag, ...
                'QRS_duration', QRS_duration, ...
                't_wave_duration', t_wave_duration, ...
                'QT_duration', QT_duration, ...
                't_peak_idx', t_peak_idx, ...
                't_wave_start_idx', t_wave_start_idx, ...
                't_wave_end_idx', t_wave_end_idx, ...
                't_magnitude_true', t_magnitude_true, ...
                'isInvertedT', logical(inverse), ... 
                'QRS_start_time', time(QRS_start_idx), ...
                'QRS_start_amp',  V(QRS_start_idx), ...
                'QRS_end_time',   time(QRS_end_idx), ...
                'QRS_end_amp',    V(QRS_end_idx), ...
                't_start_time',   time(t_wave_start_idx), ...
                't_start_amp',    V(t_wave_start_idx), ...
                't_peak_time',    time(t_peak_idx), ...
                't_peak_amp',     V(t_peak_idx), ...
                't_end_time',     time(t_wave_end_idx), ...
                't_end_amp',      V(t_wave_end_idx), ...
                'isPartialBeat',  isPartial ...
            );

    
            % Store into struct
            resObj = resObj.addEntry(lead, entry);
            beatAnalysis.(char(lead)) = entry;
            % Choose QT from primary lead (V5, fallback V6)
            if any(lead == primaryQTLeads)
                ok = isfinite(entry.QT_duration) && isfinite(entry.QRS_duration) && ...
                     entry.QT_duration >= QT_sanity.minQT && entry.QT_duration <= QT_sanity.maxQT && ...
                     entry.QRS_duration >= QT_sanity.minQRS && entry.QRS_duration <= QT_sanity.maxQRS;
            
                % Prefer first lead in primaryQTLeads (V5) then V6
                if ok
                    if primaryQTLead == "" || lead == primaryQTLeads(1)
                        primaryQT = entry.QT_duration;
                        primaryQTLead = lead;
                        primaryQTok = true;
                    end
                end
            end

            subplot(3,2,i);
            plotCustomDelineation(time, V, ...
                entry.QRS_start_time, entry.QRS_start_amp, ...
                entry.QRS_end_time,   entry.QRS_end_amp, ...
                entry.t_start_time,   entry.t_start_amp, ...
                entry.t_peak_time,    entry.t_peak_amp, ...
                entry.t_end_time,     entry.t_end_amp, ...
                leadNames(i));
        end

        % Save outputs for this beat
        base = sprintf('%s_ECG_Beat_%03d', simName, j);
        matPath = fullfile(resultsDir, base + ".mat");
        pngPath = fullfile(resultsDir, base + ".png");
        figPath = fullfile(resultsDir, base + ".fig");

        save(matPath, 'resObj', 'beatTbl', 'beatAnalysis');
        saveas(fig, pngPath);
        savefig(fig, figPath);

        % Store into simResult
        simResult.beats(j).index     = j;
        simResult.beats(j).table     = beatTbl;
        simResult.beats(j).analysis  = beatAnalysis;
        simResult.beats(j).resObj    = resObj;
        simResult.beats(j).files     = struct('mat',matPath,'png',pngPath,'fig',figPath);
        simResult.beats(j).isPartial = isPartial;
        simResult.beats(j).startIdx  = s;
        simResult.beats(j).endIdx    = e;
        simResult.beats(j).QT_primary_ms   = primaryQT;
        simResult.beats(j).QT_primary_lead = primaryQTLead;
        simResult.beats(j).QT_primary_ok   = primaryQTok;
    end

    % Put this simulation in the Map
    allResults(char(simName)) = simResult;
end

save(fullfile(resultsDir, "AllSim_ECG_Results_Map.mat"), 'allResults');
T = allResultsToTable(allResults, analyseBeat, analyseBeatEnd);
writetable(T, fullfile(resultsDir, 'AllResults.csv'));
disp("Done. Results stored in `allResults` (containers.Map) and 'AllResults.csv' and saved to disk.");

%% Functions

function ECG_table = buildECGTable(ecgData, transUnits)
    t  = transUnits * ecgData(:, 1);
    LA = transUnits * ecgData(:, 2);
    RA = transUnits * ecgData(:, 3);
    LL = transUnits * ecgData(:, 4);
    RL = transUnits * ecgData(:, 5);
    V1 = transUnits * ecgData(:, 6);
    V2 = transUnits * ecgData(:, 7);
    V3 = transUnits * ecgData(:, 8);
    V4 = transUnits * ecgData(:, 9);
    V5 = transUnits * ecgData(:, 10);
    V6 = transUnits * ecgData(:, 11);

    Vw  = (RA + LA + LL) / 3;
    I   = LA - RA;
    II  = LL - RA;
    III = LL - LA;
    aVR = 3/2 * (RA - Vw);
    aVL = 3/2 * (LA - Vw);
    aVF = 3/2 * (LL - Vw);

    V1 = V1 - Vw;  V2 = V2 - Vw;  V3 = V3 - Vw;
    V4 = V4 - Vw;  V5 = V5 - Vw;  V6 = V6 - Vw;

    max_leads = max(abs([I;II;III;aVR;aVL;aVF;V1;V2;V3;V4;V5;V6]));
    if max_leads == 0 || isnan(max_leads)
        max_leads = 1;
    end

    I_means   = I   / max_leads;
    II_means  = II  / max_leads;
    III_means = III / max_leads;
    aVR_means = aVR / max_leads;
    aVL_means = aVL / max_leads;
    aVF_means = aVF / max_leads;

    V1  = V1  / max_leads;
    V2  = V2  / max_leads;
    V3  = V3  / max_leads;
    V4  = V4  / max_leads;
    V5  = V5  / max_leads;
    V6  = V6  / max_leads;

    ECG_table = table(t, I_means, II_means, III_means, aVR_means, aVL_means, aVF_means, ...
                      V1, V2, V3, V4, V5, V6);

    ECG_table.Properties.VariableNames{1} = 't';
end

function idx = clampIndex(idx, n)
    if isempty(idx) || ~isfinite(idx)
        idx = 1;
        return
    end
    idx = round(idx);
    idx = max(1, min(n, idx));
end

function T = allResultsToTable(allResults, analyseBeat, analyseBeatEnd)

    simNames = keys(allResults);

    rows = {};
    r = 0;

    for s = 1:numel(simNames)
        simName = simNames{s};
        sim = allResults(simName);

        RR_s = sim.cycleLength / 1000;

        for b = analyseBeat:analyseBeatEnd
            beat = sim.beats(b);
            leadNames = fieldnames(beat.analysis);

            % Primary QT + QTcF
            QT_primary_ms   = NaN;
            QT_primary_lead = "";
            QT_primary_ok   = false;
            QTcF_primary_ms = NaN;

            if isfield(beat, 'QT_primary_ms')
                QT_primary_ms   = beat.QT_primary_ms;
                QT_primary_lead = string(beat.QT_primary_lead);
                QT_primary_ok   = beat.QT_primary_ok;

                if isfinite(QT_primary_ms) && isfinite(RR_s) && RR_s > 0
                    QTcF_primary_ms = QT_primary_ms / (RR_s)^(1/3);
                end
            end

            for l = 1:numel(leadNames)
                leadName = leadNames{l};
                A = beat.analysis.(leadName);

                r = r + 1;

                % Per-lead QTcF
                QTcF_ms = NaN;
                if isfinite(A.QT_duration) && isfinite(RR_s) && RR_s > 0
                    QTcF_ms = A.QT_duration / (RR_s)^(1/3);
                end

                rows(r,:) = { ...
                    simName, ...
                    b, ...
                    leadName, ...
                    A.isPartialBeat, ...
                    QT_primary_ms, ...
                    QTcF_primary_ms, ...
                    A.isInvertedT, ...
                    QT_primary_lead, ...
                    QT_primary_ok, ...
                    A.QRS_duration, ...
                    A.QT_duration, ...
                    QTcF_ms, ...
                    A.t_wave_duration, ...
                    A.QRS_start_time, ...
                    A.QRS_end_time, ...
                    A.t_start_time, ...
                    A.t_peak_time, ...
                    A.t_end_time, ...
                    A.t_magnitude_true ...
                };
            end
        end
    end

    T = cell2table(rows, 'VariableNames', { ...
        'Simulation', ...
        'Beat', ...
        'Lead', ...
        'IsPartialBeat', ...
        'QT_primary_ms', ...
        'QTcF_primary_ms', ...
        'T_wave_inversion', ...
        'QT_primary_lead', ...
        'QT_primary_ok', ...
        'QRS_duration_ms', ...
        'QT_duration_ms', ...
        'QTcF_ms', ...
        'Twave_duration_ms', ...
        'QRS_start_time_ms', ...
        'QRS_end_time_ms', ...
        'T_start_time_ms', ...
        'T_peak_time_ms', ...
        'T_end_time_ms', ...
        'T_peak_amplitude' ...
    });
end
