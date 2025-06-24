# The *latenZy* repository

Welcome to the *latenZy* repository — a Python and MATLAB toolbox containing two novel, binning-free methods for estimating the onset of neural spiking activity: `latenZy` and `latenZy2`.

Our preprint describing these methods is now online: ...

## Overview

To address key limitations in existing latency estimation approaches, we developed `latenZy` and `latenZy2`, two methods that estimate the onset of neural responses with high temporal precision—without relying on arbitrary binning or thresholds.

- **`latenZy`** estimates **when neural responses begin following discrete events** by detecting event-locked changes in spiking rates.

- **`latenZy2`** identifies the time point at which neural spiking **begins to diverge between experimental conditions**.

### Estimating Response Latency with `latenZy`

`latenZy` estimates when neural responses start after an event, given spike times and event times. Below are usage examples in MATLAB and Python.

**Python example:**
```python

```

**MATLAB example:**
```matlab
[L, sLatenzy] = latenzy(spikeTimes, eventTimes, useMaxDur);
fprintf('Estimated latency: %.2f ms\n', L);
```
