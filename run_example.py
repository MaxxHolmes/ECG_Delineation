from scipy import io
from delineate_ecg import delineate_ecg

lead = 'V1'
width = 3
activationTime = 67
ECG_table = io.loadmat('ecg_table.mat')  # This is the same table as all_1_table.mat, however, it is python compatible.
delineate_ecg(ECG_table=ECG_table['ecg_table'], lead=lead, width=width, activationTime=activationTime)
