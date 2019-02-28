## HST Bin and Sim

This takes results from the online ETC and bins and simulates time series.

### Instructions
Steps:

 * Calculate SNR using the online ETC and download the CSV table of results
 * Place the results in `calcs/`
 * Edit the observing parameters, using `params/planet_params.yaml` as a guide.
 * Run on the command line with `python bin_and_sim.py`

## I checked this against Pandexo
In one example (XO-3 b eclipse tested on 2019-02-28), Pandexo had about sqrt(2) times higher noise.
