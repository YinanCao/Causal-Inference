## Causal-Inference
Data/code for 'Causal inference in the multisensory brain' (https://doi.org/10.1016/j.neuron.2019.03.043)

## Abstract
When combining information across different senses, humans need to flexibly select cues of a common origin while avoiding distraction from irrelevant inputs. The brain could solve this challenge using a hierarchical principle by deriving rapidly a fused sensory estimate for computational expediency and, later and if required, filtering out irrelevant signals based on the inferred sensory cause(s). Analyzing time- and source-resolved human magnetoencephalographic data, we unveil a systematic spatiotemporal cascade of the relevant computations, starting with early segregated unisensory representations, continuing with sensory fusion in parietal-temporal regions, and culminating as causal inference in the frontal lobe. Our results reconcile previous computational accounts of multisensory perception by showing that prefrontal cortex guides flexible integrative behavior based on candidate representations established in sensory and association cortices, thereby framing multisensory integration in the generalized context of adaptive behavior.

### `Data_15subjs_22Trls_MEGextract.mat`: behavioural data for 4-choice rate categorisation
- `Data_runs`: [1672 x 11 x 15] = [Trials x Vars x Subjects]

   Column Vars represent:
  - Task (Aud = 0/Vis = 1)
  - Aud reliability (High = 1/Low = 2/Uni = 3)
  - Vis rate (1,2,3,4 increasing)
  - Aud rate (1,2,3,4 increasing)
  - Response (1,2,3,4 increasing)
  - Reaction time (s)
  - Onset sample (Fs = 1017.25Hz)
  - Offset sample (Fs = 1017.25Hz) 1 [1/60] blank frame after 0.55s
  - SessionID
  - BlockID
  
- `Data_sum`: [76 x 9 x 15] = [Conditions x Vars x Subjects]; Last 12 conditions are unisensory
   
   Column Vars represent:
  - Task (Aud = 0/Vis = 1)
  - Aud reliability (High = 1/Low = 2/Uni = NaN)
  - Vis rate (9.09, 12.7, 16.36, 20 Hz; NaN for respective unisensory)
  - Aud rate (9.09, 12.7, 16.36, 20 Hz; NaN for respective unisensory)
  - Count of slowest choice (9.09)
  - Count of slow choice (12.7)
  - Count of fast choice (16.36)
  - Count of fastest choice (20 Hz)
  - Total number of trials per condition


### `Model.mat`: Modelling results
- `Modelpred`: [15 x 64 x 4] = [Subjects x Conditions x Models (CI-MA, Fusion, SA, SV)], units of Hz
- `modelRDM`: [15 x 2016 x 4] = [Subjects x Pairs x Models (CI-MA, Fusion, SA, SV)], ranked Euclidean distance
- `parameter`: 6 cells, model order: (CI-MA, CI-PM, CI-MS, Likelihood, Fusion, Seg) see Table S1.
   
   Columns: 
   - SD_p
   - mu_p
   - SD(a1 high)
   - SD(a1 low)
   - SD(v1)
   - SD(a4 high)
   - SD(v4)
   - ka
   - kv
   - p(common) or pc
   - Log-likelihood
   - R-sq


### `MEG_ROI_RDM.mat`: cross-validated Mahalanobis distance (5 folds)
- `MEG_ROI_RDM`: [2016 x 24 x 15] = [Pairs x ROIs x Subjects]
- `MEG_ROI_name`: Anatomical label, stim or resp-locked
