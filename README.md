# Causal-Inference
data/code for 'Causal inference in the multisensory brain'

`Data_15subjs_22Trls_MEGextract.mat`: behavioural data for 4-choice rate categorisation
- `Data_runs`: [1672 x 11 x 15] = [Trials x Vars x Subjects]

   Column Vars represent:
  - task (Aud = 0/Vis = 1)
  - Aud reliability (High = 1/Low = 2/Uni = 3)
  - Vis rate (1,2,3,4 increasing)
  - Aud rate (1,2,3,4 increasing)
  - Response (1,2,3,4 increasing)
  - reaction time (s)
  - onset sample (Fs = 1017.25Hz)
  - offset sample (Fs = 1017.25Hz) 1 [1/60] blank frame after 0.55s
  - sessionID
  - blockID
 - `Data_sum`: [76 x 9 x 15] = [Conditions x Vars x Subjects]
   Last 12 conditions are unisensory
   Column Vars represent:
  - task (Aud = 0/Vis = 1)
  - Aud reliability (High = 1/Low = 2/Uni = 3)
  - Vis rate (9.09, 12.7, 16.36, 20 Hz)
  - Aud rate (9.09, 12.7, 16.36, 20 Hz)
  - Count of slowest choice (9.09)
  - Count of slow choice (12.7)
  - Count of fast choice (16.36)
  - Count of fastest choice (20Hz)
  - Total number of trials per condition
