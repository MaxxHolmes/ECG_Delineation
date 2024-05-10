# ECG_Delineation
Simulated ECG Delineation Code. Hannah Smith and Maxx Holmes; 2024.

Custom ECG delineation script that works on gradient analysis of sliding windows. ECG delineation is a tricky problem, particularly for ECGs with significant pathology. Off the shelf solutions did not work, particularly on chronic infarct cases, so we made a solution that capitalises on various features of these ECGs (like no P waves, no recording noise and a T wave that goes back to zero with no DADs). 

How this solution works:
1. Firstly, determine QRS start time
    Objective: find the time when the gradient is changing enough to be a wave
      Find the second derivative of the voltage vector with respect to time (gradient function twice)
      Calculate the mean of ddV in a window with centre of timestep 4 and width of 3 in either direction (ie values 1 to 7)
      Check if the absolute value of this is above a threshold (that depends on the maximum absolute value of V), if it is then the gradient of V is changing quickly enough that the QRS complex has started – this would have a start time of the first timestamp
      If not, move the centre of the sliding window in the positive t direction until this criteria has been met

2. Then determine QRS end time (harder)
    Objective: determine the point at which the gradient is relatively constant and low (but note that the ST segment is not guaranteed to be straight or constant)
      For timestamps of 30ms onwards, find the mean sliding window (width 2) of the absolute of ddV
      All window centres with values below a threshold (that depends on the maximum absolute value of the voltage) are taken as candidates for the QRS end time
      For each candidate, the sliding window of the absolute of dV is taken (width 5) and if it is lower than a threshold value then take this as the QRS end time
      A flag is raised if the QRS end time is 25ms or more different from the maximum activation time. The maximum activation time was calculated as any activation time more than three scaled median absolute deviations (MAD) from the median (wording taken   from matlab website) – this is the same however for every lead so this is likely to throw up discrepancies when there is significant dispersion between leads from an infarction
      A flag is also raised when all of the 43 translations and rotations are finished if the dispersion in the QRS end time is 25ms or higher. Again, this can happen for genuine reasons, but all the traces need to be checked for if there is a situation where small changes in the ECG cause unreasonably large jumps in the QRS end delineation

3. Then determine T peak:
    Objective: just the peak (either positive or negative) after the QRS and ST segment are finished
      The peak T wave is the maximum value of the absolute of the voltage from 100ms after the QRS end to the end of the beat – this gives some margin for failures in the QRS end time delineation and screwy ST segments
 
4. Then find T wave end:
    Objective: find the last point where there is significant change in the voltage
      Go backwards from the end of the signal, working out the sliding window of absolute of dV (width 3)
      The first point (latest timestamp) where the gradient is above a threshold is the end of the T
      Note this works because the model is incapable of producing DADs
      If the algorithm failed to find a T end (happens when the T wave is incredibly small) then the T wave was manually delineated
 
Composite values:
  QRS duration = QRS end – QRS start
  T peak end = T end – T peak
  QT interval = T end – QRS start
