## To install libraries:

Run

```
Rscript config.R
```

## To run the code:

Run

```
sh src/web.sh INPUT_FOLDER OUTPUT_FOLDER
```

to read from the `INPUT_FOLDER` and write to the `OUTPUT_FOLDER`. 

INPUT_FOLDER must contain:
1. tree.nwk, newick tree
2. info.csv, comma delimited file with the columns:
    - ID (tips names of the newick tree)
    - Date (collection date of sequences)
    - Query (0 if training and 1 if censored)
3. runid.txt, file containing a run identifier

OUTPUT_FOLDER will contain: 
1. rooted_tree.nwk, newick tree rooted by RTT
2. stats.csv, one row (plus header) comma delimted file with columns:
    - RunID (run identifier)
    - dAIC (difference between null AIC and linear model AIC)
    - EstimatedRootDate (estimated date of the root as per the linear model)
    - EstimatedRootDate95(Low|High) (95% confidence interval of estimated root date)
    - EstimatedEvolutionaryRate (estimated evolutionary rate as per the linear model)
    - Fit (1 if the linear model fit (dAIC > 10 and EstimatedRootDate95Low < all dates) and 0 otherwise)
3. data.csv, comma delimited file with columns:
    - ID (sequence ID)
    - EstimatedDate (estimated date of the sequence as per the linear model)
    - EstimatedDate(Low|High) (95% confidence interval of date of the sequence)
4. divergence_versus_time.pdf, PDF file of divergence versus time plot of the sequences with linear regression and ancestral traces
5. divergence_versus_time.png, PNG version of (4)

## Note

The data folder does not contain anything that is  usable. 
