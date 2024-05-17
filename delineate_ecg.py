import numpy as np

def delineate_ecg(ECG_table, lead, width, activationTime ):
    T = ECG_table['Time'][0][0].transpose()[0]
    T = T - T[0]
    V = ECG_table[lead][0][0].transpose()[0]

    def get_window(signal, i, width):
        window = signal[(i - width):(i + width + 1)]
        return np.mean(window)

    dV = np.gradient(V)
    ddV = np.gradient(dV)

    T_ex = np.concatenate([[-3, -2, -1, 0, 1], T])
    V_ex = np.concatenate([[0, 0, 0], V])
    dV_ex = np.gradient(V_ex)
    ddV_ex = np.gradient(dV_ex)
    dV_windowed = np.zeros(300)
    for i in range(width, 301):
        dV_windowed[i - width] = abs(get_window(dV_ex, i, width))

    QRS_start_tol = 0.01 * max(abs(V)) / 30
    QRS_end_tol_ddV = 0.1 * max(abs(V)) / (30 * width + 2)
    QRS_end_tol_dV = 0.07 * max(abs(dV_windowed))
    T_start_tol = 0.12 * max(abs(V)) / 30
    T_end_tol = 0.01 * max(abs(V)) / 30

    # Determine QRS start time
    for i in range(width, 101):
        QRS_window = abs(get_window(ddV_ex, i, width))
        if (QRS_window > QRS_start_tol):
            QRS_start_idx = i - width
            break
            
    # # Determining QRS end time, and QRS duration
    QRS_window2 = np.zeros(501-30)
    for i in range(30, 501):
        QRS_window2[i - 30] = get_window(abs(ddV_ex), i, 2)
    QRS_end_idx_lst = 30 + np.where(QRS_window2<QRS_end_tol_ddV)[0] - 1

    for idx in QRS_end_idx_lst:
        dV_window = get_window(abs(dV), idx, width+2)
        if (abs(dV_window) < QRS_end_tol_dV):
            QRS_end_idx = idx - (width - 2)
            break

    QRS_start_time = T[QRS_start_idx]
    QRS_end_time = T[QRS_end_idx]
    QRS_duration = QRS_end_time - QRS_start_time

    QRS_end_flag = abs(QRS_end_time - activationTime)>=25

    segment = V[QRS_end_idx+100:1000]  # Assuming ~100 of ST segment
    t_magnitude = max(abs(segment))
    peak_idx = np.where(abs(segment) == t_magnitude)[-1][0]
    t_peak_idx = QRS_end_idx + 100 + peak_idx
    t_sign = np.sign(V[t_peak_idx])
    t_magnitude_true = t_sign * t_magnitude

    t_wave_peak_time = T[t_peak_idx]

    # Find T-wave end point
    for i in range(len(V)-width-1, t_peak_idx, -1):
        tend_window = abs(get_window(dV, i, width))
        if tend_window > T_end_tol:
            t_wave_end_idx = i
            break

    t_wave_end_time = T[t_wave_end_idx]

    # Find T-wave start point and calculate QT
    try:
        for i in range(t_peak_idx-20, QRS_end_idx, -1):
            tstart_window = abs(get_window(dV, i, width))
            if tstart_window < T_start_tol:
                if (max(abs(V[i-30:i])) < abs(V[i])) & (abs(V[i]) < (0.5 * t_magnitude)):
                    t_wave_start_idx = i
                    break
        t_wave_start_time = T[t_wave_start_idx]
    except:
        t_wave_start_idx = np.nan
        t_wave_start_time = np.nan

    t_wave_duration = t_wave_end_time - t_wave_start_time
    QT_duration = t_wave_end_time - QRS_start_time
    t_peak_end = t_wave_end_time - t_wave_peak_time
    t_start_peak = t_wave_peak_time - t_wave_start_time

    segment = V[t_wave_end_idx-10:t_wave_end_idx]
    if max(abs(segment)) > t_magnitude * 0.1:
        if max(segment) > V(t_wave_end_idx):
            inverse = True
        else:
            inverse = False
    else:
        inverse = False
    return QRS_start_idx, QRS_end_idx, QRS_end_flag, QRS_duration, t_wave_duration, QT_duration, \
        t_peak_end, t_peak_idx, t_start_peak, t_wave_start_idx, t_wave_end_idx, \
        inverse, t_magnitude_true



