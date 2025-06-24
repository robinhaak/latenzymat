# The *latenZy* repository

Welcome to the *latenZy* repository — a Python and MATLAB toolbox containing two novel, binning-free methods for estimating the onset of neural spiking activity with high temporal precision: `latenZy` and `latenZy2`. - **`latenZy`** estimates **when neural responses begin following discrete events** by detecting event-locked changes in spiking rates. - **`latenZy2`** identifies the time point at which neural spiking **begins to diverge between experimental conditions**.

Our preprint describing these methods is now online: ...

## Rationale



## Estimating Response Latency with `latenZy`

`latenZy` estimates when neural responses start relative to an event, given spike times, event times, and a time window. Below are usage examples in MATLAB and Python.

**Python example:**
```python
from latenzy import latenzy

L, s_latenzy = latenzy(spike_times, event_times, use_max_dur)
print(f"Estimated latency: {L:.2f} ms")
```

**MATLAB example:**
```matlab
[L, sLatenzy] = latenzy(spikeTimes, eventTimes, useMaxDur);
fprintf('Estimated latency: %.2f ms\n', L);
```

**Example, V1 neuron stimulated with drifting gratings:**
