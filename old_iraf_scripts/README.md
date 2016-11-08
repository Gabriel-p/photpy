
# Prepare IRAF

A `loginuser.cl` should be created in the same IRAF folder where the `login.cl`
file and the `uparm` folder exist.

```
# LOGINUSER.CL -- Customization of LOGIN.CL
#
task    $casleored = /home/<user>/IRAF/custom/casleored.cl
task    ctio = /home/<user>/IRAF/custom/ctio.cl
task    campanas = /home/<user>/IRAF/custom/campanas.cl
task    align = /home/<user>/IRAF/custom/align.cl
task    align2 = /home/<user>/IRAF/custom/align2.cl
task    standards_phot = /home/<user>/IRAF/custom/standards_phot.cl
task    standards_cat = /home/<user>/IRAF/custom/standards_cat.cl
task    standards_fit = /home/<user>/IRAF/custom/standards_fit.cl
task    psfphot  = /home/<user>/IRAF/custom/psfphot.cl
task    $daomaster_auto  = $/home/<user>/IRAF/custom/daomaster_auto
task    $daomaster_semi  = $/home/<user>/IRAF/custom/daomaster_semi
task    $daomaster_manual  = $/home/<user>/IRAF/custom/daomaster_manual
task    daom  = /home/<user>/IRAF/custom/daom.cl
task    getdata  = /home/<user>/IRAF/custom/getdata.cl
task    clean  = /home/<user>/IRAF/custom/clean.cl
task    select  = /home/<user>/IRAF/custom/select.cl

keep
```

# Instructions

## Standard stars photometry

### getdata

#### Requirements
      
 1. Working Folder: folder containing the '.fits' files to process.
 1. At least one '.fits' file MUST exist in the working folder.
 1. Will fail if a star TOO bright is present in the frame. In this case, you
 can select a section of the frame(s) to be scanned for stars, leaving the
 section with the saturated star(s) out.
      
#### What it does
      
 1. Obtains sky mean value, sky standard deviation (STDDEV) and FWHM values
 from one or several '.fits' files, through an iterative process.
 End of iteration condition: FWHM, SKY MEAN and STDDEV values must ALL have
 a difference of less than abs(10%) with the previous calculated value; OR
 after 10 or more iterations have been executed. If the iteration fails to
 converge (that is: reaches the 10 iteration limit) then an average of all
 previous values is stored in the "imname_data" file. These previous
 values can be found in the "imname_iter" file. This is a time consuming
 process and tends to fail when the image frame has illuminated borders,
 foreign white lines or over-saturated stars (in this case you could try to
 set 'datamax' to a value well bellow the saturation limit for that image
 or scan only a SECTION of the image, one that excludes the over-saturated
 star(s)).

#### Files created, copied, moved or deleted
      
 1. "imname_iter" - Stores the values of the three parameters as the
 iterative process moves forward.
 1. "imname_data" - Stores the definitive values of the parameters.

### standards_phot

#### Requirements

 1. Working Folder: where all the standard stars are located.
 1. FWHM and sky standard deviation (STDDEV) values must be known (it is
 recommended to use one single value for ALL standard stars; this value
 should be an average of all values from all frames).        
 1. "You can save yourself a lot of trouble if you simply adopt a single
 radius for all the standards from all the nights for all the filters" -
 From 'A user's guide to Stellar CCD photometry with IRAF' by Massey &
 Davis. This radius is asked at the beginning of the script and it's
 recommended to be approx 4xFWHM.
 1. Only ONE frame per filter should be inside the working folder. This
 frame should show ALL the standard stars and NOT have saturated standard
 stars. In case one single frame can not show all standards and NOT have
 un-saturated standards at the same time, then two or more frames can be
 used; REMEMBERING in the next script (ie, while creating the "imsets" file
 in 'casleored3') to add an extra field, since there CAN  NOT be two frames
 with the same filter in one single field.
 1. "cluster_folders" file, created by the previous script, should be in the
 working folder. If this file is not present, a 'Warning' will be
 presented. Care should be taken then, when (if) running following scripts
 on CASLEO images (this is not a problem when processing CTIO images).

#### What it does

 1. Runs 'daofind' and an interactive 'tvmark' to mark ONLY the standard stars.
 1. Runs 'phot' task on '.coo.*' files containing ONLY standard stars.
 1. Asks if we want to keep the '.coo.*' files (not needed).

#### Files created, copied, moved or deleted

 1. "standards_phot_data" - File that stores ALL the data input.
 1. "apert_standard" - File that stores the aperture value used for standard
 stars (which will be used by the 'psfphot' script later on). This file is
 created, copied inside the folders where star clusters are, and then
 deleted.
 1. "standars2" - This file stores the names of the standard stars files,
 without the '.fit' extension (used by the next script).
 1. '.coo.*' files are deleted or kept according to the user's choice.
 1. One '*.mag.1' file for each image (used by the next script).
