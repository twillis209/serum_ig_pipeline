import os
import re
import pandas as pd
from scipy.stats import chi2

def compile_sumher_files(input_files):
    d = []

    for x in input_files:
        trait_A, trait_B = x.split('/')[3].split('_and_')
        snp_set = x.split('/')[4]

        with open(x, 'r') as infile:
            line = infile.readline()
            line = infile.readline()

        # Category Trait1_Her SD Trait2_Her SD Both_Coher SD Correlation SD
        _, h2_A, h2_A_se, h2_B, h2_B_se, cov, cov_se, rg, rg_se = line.split()

        rg_z = float(rg)/float(rg_se)

        rg_p = chi2.sf(rg_z**2, df = 1, loc = 0, scale = 1)

        d.append(
            {
                'trait.A' : trait_A,
                'trait.B' : trait_B,
                'snp.set' : snp_set,
                'h2.A.obs.sr' : float(h2_A),
                'h2.A.obs.se.sr' : float(h2_A_se),
                'h2.B.obs.sr' : float(h2_B),
                'h2.B.obs.se.sr' : float(h2_B_se),
                'gcov.obs.sr' : float(cov),
                'gcov.obs.se.sr' : float(cov_se),
                'rg.sr' : float(rg),
                'rg.se.sr' : float(rg_se),
                'rg.z.sr' : rg_z,
                'rg.p.sr' : rg_p
            }
        )

    return pd.DataFrame(d)

def tally_predictors_from_log_files(input_files):
    d = []

    for x in input_files:
        trait_A, trait_B = x.split('/')[3].split('_and_')

        with open(x, 'r') as fh:
            line = fh.readline()

            while not re.match('Have found valid summary statistics for ?all (\\d+) predictors', line):
                line = fh.readline()

            valid_predictors_A = int(re.match('Have found valid summary statistics for ?all (\\d+) predictors', line).group(1))

            line = fh.readline()

            while not re.match('Have found valid summary statistics for ?all (\\d+) predictors', line):
                line = fh.readline()

            valid_predictors_B = int(re.match('Have found valid summary statistics for ?all (\\d+) predictors', line).group(1))

        d.append(
            {
                'trait.A' : trait_A,
                'trait.B' : trait_B,
                'valid_predictors.A' : valid_predictors_A,
                'valid_predictors.B' : valid_predictors_B
            }
        )

    return pd.DataFrame(d)
