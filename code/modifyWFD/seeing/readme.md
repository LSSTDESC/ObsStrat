The goal here is to understand the impacts of seasonal variation in seeing, as implemented by Eric N.'s model [`owsee`](https://github.com/LSSTDESC/ObsStrat/tree/owsee/code/simsee ). `owsee` essentially overwrites the seeing column in the OpSim database and updates the 5$\sigma$ point-source depth column. There are four types of owsee outputs, with tags `ow6`, `ow6c`, `ow7`, `ow7c`:
- The "6" and "7" are different (but overlappig) sets of years on Pachon.
- The "c" indicates including a rough estimate of extinction due to clouds.

We look at the seeing, 5$\sigma$ point-source depth, and 5$\sigma$ coadded depth (after Y10) based on (all four types) `owsee` outputs vs. Opsim's (only for baseline2018a and pontus_2002). Then, we take all `ow6` outputs<sup>1</sup> and implement depth+ebv cuts [as for the static forcasts](https://github.com/LSSTDESC/ObsStrat/tree/static/static/depth_cuts ) and compare the results with the OpSim dbs.

<sup>1</sup> Note that here we are only working with ten cadences: `baseline2018a`, `kraken_2026`, `kraken_2035`, `colossus_2665`, `colossus_2664`, `colossus_2667`, `pontus_2002`, `pontus_2489`, `pontus_2502`, and `mothra_2045`.