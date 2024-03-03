- [ ] Build tests for testMotifEnrich. (Probably in longtests for the iterative model)
- [x] Produce examples & example datasets
    + AR/ER Peaks. Can get seq from BSgenome
    + ESR1, ANDR + 3 motifs
- [ ] Write vignette
- [X] Scale up for lists of motifs using `unversalmotif`
- [ ] Manage duplicate rownames when using a universalmotif list
- [ ] Find points to speed up matches

(Running the vignette on 689 motifs took 21 mins to get the best matches with 4 cores. This is 30 matches/min which is too slow)
