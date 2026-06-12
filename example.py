import sys
import numpy as np
from astropy.table import Table

sys.path.append("/Users/afaisst/Work/JWST_Projects/JWST_COSMOSWebCutout_Pipeline/github/JWST_CutoutTool/")
from getJWSTcutout_CW import cutout_jwst

fn = "/Users/afaisst/Work/Projects/JWST_PAH_z1/catalogs/final_selection_z0p92to0p98_mag24_lt18.csv"

srcs = Table.read(fn)
ids = srcs["ID_UNIQUE"].copy()
srcs["use_id"] = np.asarray(ids).astype(str)

results_dict , results_all, results_cons = cutout_jwst(srcs = srcs ,
                        hduexts = ["SCI","ERR","WHT","VAR_POISSON","VAR_RNOISE","VAR_FLAT"],
                        output_path = "./output2/",
                        cutout_size_arcsec = 5,
                        overlap_fraction_limit = 0.7,
                        keynames = ["use_id","ALPHA_J2000","DELTA_J2000"],
                        verbose = 0,
                        suppress_warnings = True,
                        MAKEPLOT = True
                       )


results_cons.write("results.csv", format="csv", overwrite=True)
#print(results_cons)
#print(results_all)
#print(results_dict)