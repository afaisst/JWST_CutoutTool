import os
from webbrowser import get
import numpy as np

from astropy.io import fits
from astropy.table import Table
from astropy.nddata import Cutout2D
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS

from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from shapely import MultiPolygon

import warnings


def create_CW_tile_dictionary():

    tiles = '''A1	149.8703317	2.0856512	149.7198796	2.1403395	149.7908786	2.3354095	149.9413496	2.2807163
A2	150.0058959	2.0363591	149.8554506	2.0910612	149.9264667	2.2861269	150.07693	2.2314186
A3	150.1414523	1.9870553	149.9910155	2.0417704	150.0620479	2.2368306	150.2125019	2.1821081
A4	150.2769995	1.9377408	150.1265729	1.9924679	150.1976208	2.1875215	150.3480637	2.1327859
A5	150.4125359	1.8884166	150.2621212	1.9431545	150.3331838	2.1382005	150.4836139	2.0834528
A6	149.8045087	1.9048087	149.6540746	1.9594923	149.7250552	2.1545612	149.8755087	2.0998725
A7	149.9400575	1.8555218	149.7896293	1.9102182	149.8606274	2.1052826	150.011074	2.05058
A8	150.0755992	1.8062243	149.9251788	1.8609325	149.9961935	2.0559913	150.1466316	2.0012757
A9	150.2111325	1.7569171	150.0607214	1.8116361	150.131752	2.0066883	150.2821799	1.9519607
A10	150.3466557	1.7076011	150.1962556	1.7623299	150.2673014	1.9573744	150.4177173	1.9026358'''

    # Create tiles dictionary
    tiles_mod = tiles.replace("\t",",").split("\n")
    tiles_coords_dict = dict()
    for tt in tiles_mod:
        tid = tt.split(",")[0]
        tcoords = [(float(ttra),float(ttdec)) for ttra,ttdec in zip(tt.split(",")[1::2] , tt.split(",")[2::2]) ]
        tiles_coords_dict[tid] = tcoords

    return(tiles_coords_dict)

def get_ListOfImages():
    '''
    Returns the list of images and paths as well as identifiers. This has
    to be configures by hand. That's probably the easiest way. This has
    also to be updated from time to time as new data is coming in.
    
    '''

    primer_cosmos_mainpath = "/Volumes/MyBook_18TB/data/Work/COSMOS/COSMOS-Web/data/PRIMER-COSMOS/v0.2/"
    primer_cosmos_rawlist = '''
mosaic_nircam_f090w_PRIMER-COSMOS_epoch1_30mas_v0_2_i2d.fits	mosaic_nircam_f277w_PRIMER-COSMOS_epoch1_30mas_v0_2_i2d.fits
mosaic_nircam_f115w_PRIMER-COSMOS_epoch1_30mas_v0_2_i2d.fits	mosaic_nircam_f356w_PRIMER-COSMOS_epoch1_30mas_v0_2_i2d.fits
mosaic_nircam_f150w_PRIMER-COSMOS_epoch1_30mas_v0_2_i2d.fits	mosaic_nircam_f410m_PRIMER-COSMOS_epoch1_30mas_v0_2_i2d.fits
mosaic_nircam_f200w_PRIMER-COSMOS_epoch1_30mas_v0_2_i2d.fits	mosaic_nircam_f444w_PRIMER-COSMOS_epoch1_30mas_v0_2_i2d.fits
'''



    cosmoswebJan23_mainpath = "/Volumes/MyBook_18TB/data/Work/COSMOS/COSMOS-Web/data/images_Jan2023/v0.1/"
    cosmoswebJan23_rawlist = '''
    mosaic_nircam_f115w_COSMOS-Web_30mas_v0_1_i2d.fits	mosaic_nircam_f277w_COSMOS-Web_30mas_v0_1_i2d.fits
mosaic_nircam_f150w_COSMOS-Web_30mas_v0_1_i2d.fits	mosaic_nircam_f444w_COSMOS-Web_30mas_v0_1_i2d.fits
    '''

    cosmoswebApr23_mainpath = "/Volumes/MyBook_18TB/data/Work/COSMOS/COSMOS-Web/data/images_Apr2023/v0.2/"
    cosmoswebApr23_rawlist = '''
    mosaic_nircam_f115w_COSMOS-Web_30mas_A10_v0_2_i2d.fits	mosaic_nircam_f277w_COSMOS-Web_30mas_A10_v0_2_i2d.fits
mosaic_nircam_f115w_COSMOS-Web_30mas_A1_v0_2_i2d.fits	mosaic_nircam_f277w_COSMOS-Web_30mas_A1_v0_2_i2d.fits
mosaic_nircam_f115w_COSMOS-Web_30mas_A2_v0_2_i2d.fits	mosaic_nircam_f277w_COSMOS-Web_30mas_A2_v0_2_i2d.fits
mosaic_nircam_f115w_COSMOS-Web_30mas_A3_v0_2_i2d.fits	mosaic_nircam_f277w_COSMOS-Web_30mas_A3_v0_2_i2d.fits
mosaic_nircam_f115w_COSMOS-Web_30mas_A4_v0_2_i2d.fits	mosaic_nircam_f277w_COSMOS-Web_30mas_A4_v0_2_i2d.fits
mosaic_nircam_f115w_COSMOS-Web_30mas_A5_v0_2_i2d.fits	mosaic_nircam_f277w_COSMOS-Web_30mas_A5_v0_2_i2d.fits
mosaic_nircam_f115w_COSMOS-Web_30mas_A6_v0_2_i2d.fits	mosaic_nircam_f277w_COSMOS-Web_30mas_A6_v0_2_i2d.fits
mosaic_nircam_f115w_COSMOS-Web_30mas_A7_v0_2_i2d.fits	mosaic_nircam_f277w_COSMOS-Web_30mas_A7_v0_2_i2d.fits
mosaic_nircam_f115w_COSMOS-Web_30mas_A8_v0_2_i2d.fits	mosaic_nircam_f277w_COSMOS-Web_30mas_A8_v0_2_i2d.fits
mosaic_nircam_f115w_COSMOS-Web_30mas_A9_v0_2_i2d.fits	mosaic_nircam_f277w_COSMOS-Web_30mas_A9_v0_2_i2d.fits
mosaic_nircam_f150w_COSMOS-Web_30mas_A10_v0_2_i2d.fits	mosaic_nircam_f444w_COSMOS-Web_30mas_A10_v0_2_i2d.fits
mosaic_nircam_f150w_COSMOS-Web_30mas_A1_v0_2_i2d.fits	mosaic_nircam_f444w_COSMOS-Web_30mas_A1_v0_2_i2d.fits
mosaic_nircam_f150w_COSMOS-Web_30mas_A2_v0_2_i2d.fits	mosaic_nircam_f444w_COSMOS-Web_30mas_A2_v0_2_i2d.fits
mosaic_nircam_f150w_COSMOS-Web_30mas_A3_v0_2_i2d.fits	mosaic_nircam_f444w_COSMOS-Web_30mas_A3_v0_2_i2d.fits
mosaic_nircam_f150w_COSMOS-Web_30mas_A4_v0_2_i2d.fits	mosaic_nircam_f444w_COSMOS-Web_30mas_A4_v0_2_i2d.fits
mosaic_nircam_f150w_COSMOS-Web_30mas_A5_v0_2_i2d.fits	mosaic_nircam_f444w_COSMOS-Web_30mas_A5_v0_2_i2d.fits
mosaic_nircam_f150w_COSMOS-Web_30mas_A6_v0_2_i2d.fits	mosaic_nircam_f444w_COSMOS-Web_30mas_A6_v0_2_i2d.fits
mosaic_nircam_f150w_COSMOS-Web_30mas_A7_v0_2_i2d.fits	mosaic_nircam_f444w_COSMOS-Web_30mas_A7_v0_2_i2d.fits
mosaic_nircam_f150w_COSMOS-Web_30mas_A8_v0_2_i2d.fits	mosaic_nircam_f444w_COSMOS-Web_30mas_A8_v0_2_i2d.fits
mosaic_nircam_f150w_COSMOS-Web_30mas_A9_v0_2_i2d.fits	mosaic_nircam_f444w_COSMOS-Web_30mas_A9_v0_2_i2d.fits
    '''

    ## Start Table
    image_table = Table(names=["survey","band","tile","identifier","path"] , dtype=[str,str,str,str,str])

    ## Primer-COSMOS
    for f in primer_cosmos_rawlist.split():    
        identifier = "PC-{}-Z0-{}-{}-{}".format( f.split("_")[4],
                                            f.split("_")[5],
                                            f.split("_")[1],
                                            f.split("_")[2]
                                            )
        image_table.add_row(["PC",
                            f.split("_")[2],
                            "Z0",
                            identifier ,
                            os.path.join(primer_cosmos_mainpath , f)
                            ])
        
    ## COSMOS-Web January 2023 data
    for f in cosmoswebJan23_rawlist.split():    
        identifier = "CW-{}-Z0-{}-{}-{}".format( "Jan23",
                                            f.split("_")[4],
                                            f.split("_")[1],
                                            f.split("_")[2]
                                            )
        image_table.add_row(["CW",
                            f.split("_")[2],
                            "Z0",
                            identifier ,
                            os.path.join(cosmoswebJan23_mainpath , f)
                            ])
        

    ## COSMOS-Web April 2023 data
    for f in cosmoswebApr23_rawlist.split():    
        identifier = "CW-{}-{}-{}-{}-{}".format( "Apr23",
                                            f.split("_")[5],
                                            f.split("_")[4],
                                            f.split("_")[1],
                                            f.split("_")[2]
                                            )
        image_table.add_row(["CW",
                            f.split("_")[2],
                            f.split("_")[5],
                            identifier ,
                            os.path.join(cosmoswebApr23_mainpath , f)
                            ])
        
    return(image_table)


def cutout_jwst(srcs , hduexts , output_path , cutout_size_arcsec , overlap_fraction_limit , keynames , verbose , suppress_warnings):
    '''
    Creates cutouts for an JWST image for all sources in a catalog for which there is overlap.
    Stores the cutouts in a directory with user-defined naming.

    Note: Image paths are generated in `get_ListOfImages()`. Need to modify that if new images
    arrive! 

    Parameters
    ----------
    srcs : `astropy table`
        Path to the catalog containing the objects from which cutouts should be made
    hduexts: `list` [`int` or `str`]
        List of HDU extension names (or numbers) for which cutouts are created.
        Note that they all will end up in a single FITS as extensions. The first
        HDU extension in this list must contain the WCS information (for WCS and pixel scale).
        For example: hduexts=["SCI","ERR","WHT","VAR_POISSON","VAR_RNOISE","VAR_FLAT"]
    output_path: `str`
        Path to the directory where to save the cutouts
    cutout_size_arcsec: `float`
        Size of cutouts in ARCSECONDS
    overlap_fraction_limit: `float`
        Overlap fraction. For example 0.7 would mean 70% overlap is required. Note: do not
        set this too high (i.e., 1.0) as there could be random pixels set to 0. A value of 
        0.7 is probably best.
    keynames: `list` [`str`]
        List of key names in catalog for columns of object names, R.A., and Decl. For 
        example ["NUMBER","ALPHA_J2000","DELTA_J2000"]
    verbose: `int`
        Set to 0 to be silent. Set to 1 to talk a bit. Set to 2 to talk a lot.
    suppress_warnings: `bool`
        If set to True, suppress all warnings
    
    Returns
    --------
    Saves the image cutout in the output directory.
    Returns a catalog with statistics for each object. The catalog includes:
    - object ID
    - RA
    - DEC
    - Flag (0=good, 1=something wrong, 2=no overlap)
    
    '''
    #overlap_fraction_limit = 0.7
    nan_fraction_limit = overlap_fraction_limit*1.0
    
    
    if suppress_warnings: warnings.filterwarnings("ignore")


    ## Get list of images
    image_path_table = get_ListOfImages()

    ## Create a table for each image, key=identifier
    all_tab = dict()

    ## Create combined table (this is larger but contains all the info)
    tab_keys = [ "{}_flag".format(im["identifier"]) for im in image_path_table]
    tab_dtypes = list( np.repeat(float , len(image_path_table)) )
    all_tab_combined = Table(names=["ID","RA","DEC"] + tab_keys , dtype=[str , float , float] + tab_dtypes)

    ## Go through all images and create cutouts for all possible sources.
    for ii in range(len(image_path_table)):

        ## Start table for this image
        this_tab = Table(names=["ID","RA","DEC","identifier","flag"] , dtype=[str , float , float, str , int])

        ## This image path
        this_image_path = image_path_table["path"][ii]
        this_image_identfier = image_path_table["identifier"][ii]

        ## Open image
        with fits.open(this_image_path) as hdul:
            #hdul.info()

            if verbose >= 0: print("+++++++ Making cutouts for image {} ({}) +++++++".format(this_image_path.split("/")[-1],this_image_identfier))

            ## Cutout for each source on this image.
            for src in srcs:

                if verbose > 0: print("++ Creating Cutout for source: {} at R.A. = {:.3f} and Decl. = {:.3f}".format(src[keynames[0]],
                                                                                        src[keynames[1]],
                                                                                        src[keynames[2]]))


                ## Prepare the arrays to save the flags and hdus
                FLAGS = []
                hdus_new = []

                ## Go go through each HDU and create the cutouts
                for hh,hduext in enumerate(hduexts):
                    if verbose > 1: print(" -> Processing HDU extension {}".format(hduext))
                    
                    # get header: Note that we get the header
                    # that contains the WCS information separately. We
                    # assume that this is the first one in the list.
                    hdr = hdul[hduext].header
                    hdr0 = hdul[hduexts[0]].header
                    
                    # get WCS: same here, get the header that
                    # contains the WCS information separately. We 
                    # assume that is the first one in the list.
                    #hdr_wcs = WCS(hdul[hduext].header)
                    hdr0_wcs = WCS(hdul[1].header)

                    # Compute the pixel scale and the size of the cutouts in pixels.
                    pixscale = [hdr0["CDELT1"]*3600 , hdr0["CDELT2"]*3600]
                    size = u.Quantity((cutout_size_arcsec/pixscale[0],cutout_size_arcsec/pixscale[1]), u.pixel)

                    # create the cutouts. Assign a FLAG to it, such that
                    # FLAG = 0 if everything is OK
                    # FLAG = 1 if no cutout can be created because of any error (e.g., negative
                    # coordinate because the source is in the wrong fiels such as GOODS-N vs. COSMOS)
                    # FLAG = 2 if no cutout can be create because not enough overlap (user defined)
                    FLAG = 0
                    try:

                        # create cutout
                        position = SkyCoord(src[keynames[1]] , src[keynames[2]] , unit="degree" , frame='icrs')
                        tmp = Cutout2D(hdul[hduext].data, position, size, wcs=hdr0_wcs , copy=True, mode="partial")
                        cutout = tmp.data

                        # check overlap (if cutting is successful)
                        overlap_fraction = round(1 - len(np.where(cutout == 0)[0]) / (cutout.shape[0]*cutout.shape[1]),1)
                        fraction_nan = round(len(np.where(np.isnan(cutout))[0]) / (cutout.shape[0]*cutout.shape[1]),1)
                        if (overlap_fraction < overlap_fraction_limit) | (fraction_nan > nan_fraction_limit):
                            if verbose > 1: print("-> Not enough overlap ({}) or too many NaN ({}) to create cutout.".format(overlap_fraction,fraction_nan))
                            FLAG = 2

                    except:
                        if verbose > 1: print("-> No cutout can be made (wrong field or other error)")
                        FLAG = 1

                    
                    # Assemble this part of the HDU and append.
                    # Note that in some cases, the SCI HDU is created, but others are not
                    # because of slightly different overlap. In that case, we can probably
                    # still create the HDUL
                    if (FLAG == 0) | (FLAG == 2):
                        
                        # Make the first HDU extension the primary.
                        if hh == 0:
                            hdu_new = fits.PrimaryHDU(data = cutout.copy() , header=hdr.copy())
                        else:
                            hdu_new = fits.ImageHDU(data = cutout.copy() , header=hdr.copy())
                        
                        # update the header (just in case)
                        hdu_new.header.update(tmp.wcs.to_header())
                        hdu_new.verify('fix')

                        # Append to the HDUS list
                        hdus_new.append(hdu_new)

                        # Append flag
                        FLAGS.append(FLAG)

                        if verbose > 1: print("-> Cutout created.")

                    else:
                        
                        # no cutout can be created. Set placeholder to -99 and append flags.
                        hdus_new.append(-99)
                        FLAGS.append(FLAG)


                ## Went through all the extensions, now put everything together... if there
                # are no errors. Do not create cutout if the SCI flag is other than 0. However,
                # check that the user is requesting a SCI extension. If not, only create cutout
                # if FLAGS[0]==0 and none of the FLAGS is 2.
                if verbose > 2: print("-> Flags: ", FLAGS)
                if "SCI" in hduexts:
                    sel_hdu_test = np.where(np.asarray(hduexts) == "SCI")[0][0]
                else:
                    sel_hdu_test = 0
                
                if (FLAGS[sel_hdu_test] == 0) & (2 not in FLAGS):

                    if verbose >=0: print("-> Creating final cutout for {}".format(src[keynames[0]]))
                    
                    # cutout name
                    cutout_name = "{}-{}.fits".format( src[keynames[0]] , this_image_identfier )
                    if verbose > 0: print("-> Assembling final HDUL and save as {}".format(cutout_name))

                    # and save
                    hdul_new = fits.HDUList(hdus_new)
                    hdul_new.writeto(os.path.join( output_path , cutout_name ) , overwrite=True)
                    
                this_tab.add_row([src[keynames[0]],src[keynames[1]],src[keynames[2]],this_image_identfier,FLAGS[sel_hdu_test] ])


        ## Add this table to dictionary (containing all individual tables)
        all_tab[this_image_identfier] = this_tab.copy()



    ## Create large table
    for ii,src in enumerate(srcs):
        tmp = [ src[keynames[0]],src[keynames[1]],src[keynames[2]] ]
        for key in list(all_tab.keys()):
            tmp += [all_tab[key]["flag"][ii]]
        all_tab_combined.add_row(tmp)


    ## Create Consolidated Table
    all_tab_consolidated = Table(names=["ID","RA","DEC","survey","band","tile"] , dtype=[str,float,float,str,str,str])
    for ii,src in enumerate(srcs):
        this_row = all_tab_combined[tab_keys][ii]
        sel_good = np.where( np.asarray(list(this_row)) == 0)[0]
        if len(sel_good) == 1:
            this_survey = tab_keys[sel_good[0]].split("_flag")[0].split("-")[0]
            this_band = tab_keys[sel_good[0]].split("_flag")[0].split("-")[5]
            this_tile = tab_keys[sel_good[0]].split("_flag")[0].split("-")[2]
            all_tab_consolidated.add_row([src[keynames[0]],src[keynames[1]],src[keynames[2]],this_survey,this_band,this_tile])
        elif len(sel_good) > 1:
            this_surveys = []
            this_bands = []
            this_tiles = []
            for sel in sel_good:
                this_surveys.append( tab_keys[sel].split("_flag")[0].split("-")[0] )
                this_bands.append( tab_keys[sel].split("_flag")[0].split("-")[5] )
                this_tiles.append( tab_keys[sel].split("_flag")[0].split("-")[2] )
            all_tab_consolidated.add_row([src[keynames[0]],src[keynames[1]],src[keynames[2]],",".join(this_surveys),",".join(this_bands),",".join(this_tiles)])
        else:
            all_tab_consolidated.add_row([src[keynames[0]],src[keynames[1]],src[keynames[2]],"none","none","none"])



    if verbose > 0: print("++++++ ALL DONE FOR THIS IMAGE +++++++")


    if suppress_warnings: warnings.filterwarnings("default")

    return(all_tab , all_tab_combined, all_tab_consolidated)
