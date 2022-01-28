# Satellite Graph Lists
These files are the outputs from dataToGraph.m. Each file has one container variable. The keys are ints representing each timestep, or **the number of seconds since the start of the scenario.** The values are **adjacency matrices**. For each adjacency matrix *A*, *A<sub>i,j</sub>* is equal to the **slant range** between satellites *i* and *j* in **km** if there is an unobstructed line of sight between them. *A<sub>i,j</sub>* = 0 otherwise.

To load satellite data, make sure the file you want is added to the MATLAB path, and type **`load("list name")`**.

**Example: `load("LunaNetWithLLOGraphList");`**

## File descriptions
- **ATrainCTrainGraphList**
  - NASA's Afternoon constellation, composed of A-Train satellites (OCO-2, GCOM-W1, Aqua,and Aura) and C-Train satellites (CALIPSO and CloudSat). Propagated over 1 day (86400 seconds). Timesteps in this data set are 60 seconds.
- **GPSWalkerApproxGraphList**
  - Walker Delta constellation (similar to GPS, but not actual GPS data) composed of 24 Earth-orbiting satellites in 6 equally spaced orbital planes. Each orbit is circular and has an altitude of 20,200 km. Propagated over 1 day (86400 seconds). Timesteps in this data set are 1800 seconds.
- **MMSGraphList**
  - NASA's Magnetospheric Multiscale Mission (MMS). Propagated over 1 year (31536000 seconds). Timesteps in this data set are 42300 seconds.
- **LunaNetGraphList**
  - NASA's LunaNet concept, composed of five lunar satellites. "Two relay orbiters phased 180° apart on the 12-hour frozen elliptical orbit with its line of apsides liberating over the North Pole, plus [t]wo relay orbiters phased 180° apart on the 12-hour frozen elliptical orbit with its line of apsides liberating over the South Pole, plus [o]ne relay orbiter on the 12-hour circular orbit around the equator."[^1] Propagated over 1 day (86400 seconds). Timesteps in this data set are 60 seconds.
- **LunaNetWithLLOGraphList**
  - The same as LunaNet above, but with an additional satellite representing a user in low lunar orbit (LLO). Nodes 1-5 are LunaNet sats, and node 6 is the LLO sat.

[^1]: National Aeronautics and Space Administration, “LunaNet Concept of Operations and Architecture,” 2020.
