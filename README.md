# JWST Image Cutout Tool

A simple script to create cutouts for a list of sources from multiple JWST images. This script was optimized for COSMOS-Web but should work for any list of JWST images as long as they have similar extensions (at least one `SCI` extension).

## How to set things up

First, the images directories and names have to be set up. This is very easy.

- Open `getJWSTcutout_CW.py` and look at the function `get_ListOfImages()`. As if May 2023, there are three different paths to images, including PRIMER-COSMOS (epoch 1), COSMOS-Web January 2023 data, and COSMOS-Web April 2023 data.
- Modify the `*_mainpath` to lead to these image products.
- In your terminal, do a `ls` and copy the results in the `*_rawlist` variable. The Python `.split()` function should be able to split these strings into lists of image names.
- That's it. The code is going to cycle through these images and assigns an *identifier* to each of the images. These can be used later to refer back to the images.


## How to use it

The code is run as follows:

```
sys.path.append("/Path/to/the/code/directory/")
from getJWSTcutout_CW import cutout_jwst

results_dict , results_all, results_cons = cutout_jwst(srcs = srcs ,
                        hduexts = ["SCI","ERR","WHT","VAR_POISSON","VAR_RNOISE","VAR_FLAT"],
                        output_path = "./cutouts/",
                        cutout_size_arcsec = 5,
                        overlap_fraction_limit = 0.7,
                        keynames = ["id","ra","dec"],
                        verbose = 0,
                        suppress_warnings = True
                       )
```

The different inputs are:
- srcs : `astropy table` <br>
        Path to the catalog containing the objects from which cutouts should be made
- hduexts: `list` [`int` or `str`]<br>
        List of HDU extension names (or numbers) for which cutouts are created.
        Note that they all will end up in a single FITS as extensions. The first
        HDU extension in this list must contain the WCS information (for WCS and pixel scale).
        For example: hduexts=["SCI","ERR","WHT","VAR_POISSON","VAR_RNOISE","VAR_FLAT"]
- output_path: `str`<br>
        Path to the directory where to save the cutouts
- cutout_size_arcsec: `float`<br>
        Size of cutouts in ARCSECONDS
- overlap_fraction_limit: `float`<br>
        Overlap fraction. For example 0.7 would mean 70% overlap is required. Note: do not
        set this too high (i.e., 1.0) as there could be random pixels set to 0. A value of 
        0.7 is probably best.
- keynames: `list` [`str`]
        List of key names in catalog for columns of object names, R.A., and Decl. For 
        example ["NUMBER","ALPHA_J2000","DELTA_J2000"]
- verbose: `int`<br>
        Set to -1 to be completely silent, set to 0 to be (almost) silent. Set to 1 to talk a bit. Set to 2 to talk a lot.
- suppress_warnings: `bool`<br>
        If set to True, suppress all warnings

Some words on the `overlap_fraction_limit`: When the code creates a cutout, it computes the fraction of number of zeros or NaNs over the total number of pixels in the cutout. The `overlap_fraction_limit` defines a limit on this fraction, below which no cutout is created. A value of 0.7 (70%) is probably fine in most cases. Be careful when setting this limit to 1.0 (100%). This would reject a cutout even if there is one Nan or zero (0.0) on the image!


In addition to the cutouts, there are three outputs created:
- *results_dict*: A dictionary (key is identifier of the image) for each image including name of the source, R.A., Decl., survey (PC=PRIMER-COSMOS, CW=COSMOS-Web), band (fXXXw), tile (AX, for CW January data and PC it is a placeholder Z0), and a flag (see below).
- *results_all*: This is the same as *results_dict*, however, all in one table. It for each source (i.e., row), it includes name, R.A., Decl. and flags for each image (called [identifier]_flag). 
- *results_cons*: This is the consolidated table, probably __the most useful output__. For each source, it contains name, R.A., Decl., survey, band, and tile. The latter three (survey, band, tile) contain one entry for each image in which the source is covered (delimited by a comma). For example, the survey = "CW,CW,PC", band = "f150w,f277w,f090w", etc. If a source is not covered, they value is "none".

As mentioned above, different flags are included:
- `FLAG = 0`: all is good
- `FLAG = 1`: A source is not contained in the image or there is another error (basically when `Cutout2D()` fails)
- `FLAG = 2`: A source is contained on the image, however, there is no coverage. This flag is based on the `overlap_fraction_limit` set by the user (see above).








