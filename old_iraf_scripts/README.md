
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









6-) CASLEORED3 ; This script is executed by typing 'casleored3' in IRAF

      1- REQUIREMENTS TO RUN THIS SCRIPT

      - Working Folder: same as above (ie: the 'standard/' folder).
      
      * The file "standars2" created by the previous script must exist.
      * The filters identifiers in the images headers MUST be "U,B,V,I".
      Otherwise the script needs to be modified: when running the 'mknobsfile'
      task, replace "U,B,V,I" by the appropriate set of identifiers or update
      the images headers with this set of filters.
      * The "*.mag.1" files created by the previous script must exist in the
      working folder.
      
      2- WHAT IT DOES
      
      * Runs 'mkcatalog' task unless the user wants to use an existing
      'landolt.cat' file (if so, a 'flandolt.cat.dat' file must also be provided
      by the user). Both files are used by the next script.
      * Creation of the 'sets' file (used by the 'mknobsfile' task).
      * Runs 'mknobsfile' task.
      * IMPORTANT --> After this scripts ends, the file "noche.obs" MUST be
      edited so that the names of the standard stars matches those in the
      "landolt.cat" file.

      3- FILES CREATED, COPIED, MOVED OR DELETED
      
      * "landolt.cat" - This file is either created by the 'mkcatalog' task or
      by the user (file used by the next script).
      * "flandolt.cat.dat" - File that contains the format of "landolt.cat"          
      * "sets" - File used by the 'mknobsfile' task (not used by any other
      script, but still not deleted for control purposes).
      * "noche.obs" - File created by the 'mknobsfile' task (used by the next
      script).
      * "fnoche.obs.dat" - File that contains the format of "noche.obs"      


7-) CASLEORED4 ; This script is executed by typing 'casleored4' in IRAF

      1- REQUIREMENTS TO RUN THIS SCRIPT
      
      - Working Folder: same as above (ie: the 'standard/' folder).
            
      * The following files MUST exist in the working folder: noche.obs,
      fnoche.obs.dat, flandolt.cat.dat and landolt.dat. An existing version of
      the noche.cfg file can be used (otherwise one can be created by the
      'mkconfig' task).
      * The names of the stars in the 'noche.obs' file MUST match those names in
      the 'landolt.cat' file. 
      * "cluster_folders" file, created by the previous script, should be in the
      working folder. If this file is not present, a 'Warning' will be
      presented. Care should be taken then, when (if) running following scripts.
      * Check that the RMS is in the hundredth, standard deviation is 0.02-0.03
      and reduced chi is approx. 1.
      * You can store the data with ':vshow fitp_file' (this is not necessary,
      since a 'log' file with this name is created automatically)
      
      2- WHAT IT DOES
      
      * Allows the user to use an existing 'noche.cfg' file or else run
      'mkconfig' task to create one.
      * Runs 'fitparams' task. The 'standard deviation' should not be bigger
      than 0.03, the RMS should be in the hundredth and the 'reduced chi' value
      should be close to 1.
      * Stores a copy of the files "noche.cfg" and "noche.ans" in each of the
      cluster's folders (if they exist).
      
      3- FILES CREATED, COPIED, MOVED OR DELETED
      
      * "noche.cfg" - This file is either created by the 'mkconfig' task or by
      the user. A copy of this file is stored in each of the cluster's folders 
      (if the file 'clusters_folders' exists).
      * "noche.ans" - File created by the 'fitparams' task (file used by the
      last script). A copy of this file is stored in each of the cluster's
      folders (if the file 'clusters_folders' exists).
      * A 'log' file called: 'fitp_file' 



## PSF photometry
1. Run `psfphot`, one time per image frame.
1. Once *all* cluster's image frames have been processed by `psfphot`, run
`daom` script.

9-) PSFPHOT ; This script is executed by typing 'psfphot' in IRAF

      1- REQUIREMENTS TO RUN THIS SCRIPT
      
      - Working Folder: this script is to be run INSIDE each of the cluster's
      folders.
      
      * The script can be used under 'interactive' mode or completely automatic.
      The latter should be used with caution.
      
      * Image frames from different filters MUST be aligned. 
      
      * If used in 'automatic' mode, then the following will be used as default:
      
      datapars.noise = "poisson"
        datapars.ccdread = "RDNOISE"
        datapars.gain = "GAIN"
        datapars.exposur = "EXPTIME"
        datapars.filter = "FILTERS"
        
      If used under 'interactive' mode, these same values will be used, but the
    user is given the choice to accept them or exit the script and set them
    manually. If the user wishes to use different values, it is better to set
    them BEFORE running this script to avoid having to exit and re-run it
    after setting them.
      
      * This script can process a single image or a file containing a list of
      images. If the user wants to process a list of images in order, then a
      file must be created containing the names of that images one name by line.
      Then, when asked for the 'Image name' by the script simply introduce the
      name of that file with the '@' character at the beginning. Eg: if the file
      is named 'image_list', then it will look like this:
        
        image1.fits
        image2.fits
        image3.fits
        ...
        
        (the dots at the end indicate more images) and when asked by the script
        to "Input image name" the user will write: @image_list
        
      * The name of the image files should be introduced without the
      '.fits'/'.fit' extension (the same applies to image names written inside
      the input file if one is used) The script should work nonetheless if the
      extension is included.
        
      * 'datamax' value MUST be known since it will be asked at the beginning of
      the script.
        
      * GAIN and RDNOISE values can be supplied by a file called 'gain_rdnoise'
      that must exist in the same folder as the images. This file should contain
      only ONE line with the parameters GAIN and RDNOISE (in THAT order) present
      in the first line, separated by a single space. If the file 'gain_rdnoise'
      does not exist, the script will attempt to read these values from the
      image header, this means that GAIN and RDNOISE values MUST exist in the
      image header. * 'aperture' value used for standard stars should reside in
      a file called 'apert_standard'. If this file doesn't exist, the value will
      be asked and the file containing that value will be created in that folder
      (THIS MEANS THE APERTURE VALUE USED FOR STANDARD STARS MUST BEKNOWN). Only
      ONE aperture value for standard stars can be inputed. If MORE than one
      value of aperture was used for standard stars previously, then the code
      might need to be re-written.
        
        ^^ (If all scripts worked correctly (for CASLEO frames), both this files
        should already exist inside each of the cluster's folders)
        
        
      * A file 'imname_data' will be expected by the script (where 'imname' is
      the name of the image being processed; without any extension) This file
      should contain FWHM, SSTDEV and SKY MEAN values ONE PER LINE and in that
      order. If the file doesn't exist, the script will ask for this values.
        
        ^^ sky mean value , sky standard desviation value [STDDEV] and FWHM are
        calculated by the 'getdata' script. The radius of the bigger PSF star in
        the frame and aperture photometry radius to be used [pixels] are
        calculated by 'psfphot'.
           
        
      * 'psfrad' value (radius of bigger PSF star in the frame) will be
      calculated by the script as (4xFWHM + 1). The user can accept or reject
      this value (in the case the user reject the value, he will have to input
      that value manually). The same applies for the 'aperture' value (aperture
      photometry radius [pixels]), which will be calculated as (0+FWHM)
        
      * If 'smean*gain + rdnoise*rdnoise < 0', 'sigma =0' will be set. Otherwise
      'sigma = sqrt(smean*gain + rdnoise*rdnoise)/gain'
        
      * A 'WARNING' will be issued if 'sigma' and 'SSTDEV' values differ
      significantly --> a difference of less than abs(25%) is acceptable 
      (totally arbitrary criteria). Page 27 of 'A Reference Guide to the
      IRAF-DAOPHOT Package' by L. Davis, recommends to check for this
      difference.
        
      * According to a 'A user's guide to stellar CCD photometry with IRAF' by
      Massey & Davis (page 41): '...set datamin in datapars to some value like 
      (sky value) - 3*sigma...'. According to IRAF's datapar help: '...DAOPHOT
      users should either leave datamin set to INDEF or set it to a number
      between 5-7 sigma below the sky background value'. So we set this value to
      'smean-4.5*sigma', which the user can accept or reject (in which case a
      value will be asked for)
      
      2- WHAT IT DOES
      
      * Basically it follows the guide 'A user's guide to stellar CCD
      photometry with IRAF' by Massey & Davis (1992); although some things are
      changed. The steps are as follow:
        
        1. Use display and imexamine on a few frames to determine the typical
        full-width-half-max of stars (FWHM) and what would be a good value to
        use for the radius of the psf (i.e., what radius will contain the
        brightest star for which you wish to do photometry.) Also determine the
        rest of the necessary data (sky level, etc).
        2. Find stars.
           2.1. Determine 'delta' (delta = sqrt(smean*gain +
           rdnoise*rdnoise)/gain)
           2.2. Enter the sky value minus 3*delta as your value for datamin in
           datapars.
           2.3. Run daofind using as a threshold value 3.5*delta.
           2.4. Use tvmark to mark the stars found (imagename.coo.1). If you
           need to, rerun daofind with a larger or small threshold.
        3. Run aperture photometry using phot.
        4. Automatic search for PSF stars with 'pstselect' task (also, option
        for manual selection).
        5. Run 'psf' and add stars using the "a" key. Try to select bright,
        uncrowded stars.
           Then:
           5.1. Run 'nstar'. If there are neighbors, be sure to decrease the psf
           radius. Run 'substar' (also using the smaller sized psf radius) and
           display the resultant subtracted frame. Do the residuals of the PSF
           stars look consistent, or is one of them funny? If need be, start
           over.
           5.2. Remove any neighbor stars by editing the PSF stars out of the
           ".nst" file, and rerunning 'substar'. Run 'psf' on the subtracted
           file, using the normal psf radius again. Rerun 'nstar' on the
           original frame using the normal psf radius and the revised PSF. Run
           'substar' and display the results. Are the PSF stars nicely removed,
           and do the areas around the PSF stars look clean? If need be, run new
           iteration.
        6. Run 'allstar' on the original frame. Display the substracted frame
        and see if your stars have been nicely substracted off. If need be star
        over from PSF selection.
        7. Run daofind on the subtracted frame. Use tvmark to examine the
        results, and if need be add any extra stars. Add new stars, previously
        hidden (if any) and perform 'phot' and 'allstar' again.
        8. Compare .coo.1, .als.1 and .als.2 (if such file exists) files and
        either choose one of these files to keep or else rerun the script from
        the beginning ("2. Find stars").
        9. Calculation of aperture correction and use of the task 'pcalc' to
        correct the 'MAG' value in the '.als' or '.mag' file by adding the
        previously calculated aperture.
        10. Creation of the 'list' file (which stores the names of all the
        '.als' files), to be used by the 'daom' script next.
        11. Rerun the script for each frame in your "set" (e.g., all short and
        long exposures in each filter of a single field, say).
   
      3- FILES CREATED, COPIED, MOVED OR DELETED
      
      * "list" - List of the '*.als.*' files created. File used by the next
      script.
      * "aperture" - Control file that stores the aperture values for all the
      frames (not used by any script)
      * "imname_aperture" - Control file that stores the aperture values for
      every PSF star used.
      * "imname.parameters" - Control file that stores the relevant parameters
      used during the PSF photometry for the 'imname' image frame.
      * The most important files: '*.coo' files and '*.als' files  
      * A bunch of extra files '.sub', '.pst', etc.., are created (see
      "files_created_by_psfphot" for a larger description).


10-) DAOM ; This script is executed by typing 'daom' in IRAF

      1- REQUIREMENTS TO RUN THIS SCRIPT
      
      - Working Folder: this script is to be run INSIDE each of the cluster's
      folders, once ALL the cluster's frames in that folder have been processed
      by 'psfphot'.

      * The codes 'daomaster_auto', 'daomaster_semi' and 'daomaster_manual',
      MUST be compiled prior to running this script and they need to be located
      in the '/home/IRAF/custom/' folder. If they are NOT located in
      this folder then the 'loginuser.cl' (showed at the beginning of this
      manual) must be modified accordingly, if used.
          
      * It needs the files: 'list', the '*.als.1' (or '*.als.2') or '*.mag.1' 
      (if you only performed aperture photometry) files and 'noche.ans' file. It
      also needs the '*.fit' image frames (the same ones that are stored in the
      'list' file), because the script needs to extract from the headers of this
      files the information about 'FILTERS' and 'EXPTIME'. The file 'noche.cfg'
      can (and should) be created by this script (this is to avoid formatting
      errors, since this file hasn't the same format as the '.cfg' file created
      by the 'mkconfig' task).
      
      2- WHAT IT DOES
      
      * Converts the '.als' or '.mag' files to 'DAOMASTER' format ('.txt'
      files).
      * Identifies filter and exposure time information and creates one '.mch'
      file per filter, with the first frame in the '.mch' files being the one
      with the longest exposure.
      * Creates the file 'daom.mch'. This file contains the names 'vfilter.mag',
      'bfilter.mag', 'ufilter.mag' and 'ifilter.mag', IN THAT ORDER (ie: with
      the 'V' filter file in first place), with the correct file format.
      * Feeds this created files ('*.mch') to DAOMASTER (which can be used in an
      automatic, semi-automatic or manual way). The output files are:
      'ufilter.mag', 'vfilter.mag', 'bfilter.mag', 'ifilter.mag' and 'daom.raw'
      * Adds airmass information to 'daom.raw' file (the airmass for each filter
      is the one that corresponds to the frame with the longest exposure). The
      output is 'daom.obs' file.
      * Deletes all the '*.mch' files. Optionally deletes all the '*.txt' files.
      * Runs 'invertifit' task TWO times. The first time we get a
      'night_obs.out.BV' file, the second time the output is a
      'night_obs.out.VI' file. The difference with this second file is that, to
      obtain it, only the equations for 'V' and 'I' filters (in the 'noche.cfg'
      file) were used.
      3- FILES CREATED, COPIED, MOVED OR DELETED
      
      * "noche.cfg" - If this file doesn't exist already, it can be created by
      this script. This file is later renamed to 'noche.cfg.BV' and another file
      is created: 'noche.cfg.VI'
      * Creates one '.txt' file per '.als' file.
      * Creates one '.mch' file per filter and one 'daom.mch' file (deleted at
      the end of the script).
      * Creates one '.mag' file per filter and one 'daom.raw' file.
      * Creates one '.cor' file per filter and one 'daom.cor' file.
      * Creates one 'daom.obs' file.
      * Creates the FINAL OUTPUT FILES 'night_obs.out.BV' and 'night_obs.out.VI'
      (created by 'invertfit')
      * Creates the 'shift_entrada.BV' file to feed the 'shift' code
      * Creates one 'shift_entrada.BV' file (this is the the 'night_obs.out.BV'
      file converted to the format needed to feed the 'shift' code; stars with
      INDEFS in any magnitude are not present in this file)
      * Creates one 'exptime.found' file, which stores the longest exposures
      times found. ONLY IF 'Automatic mode' is selected.  


## Prepare images (CTIO/Campanas/Casleo)

1. Run the 'setairmass' task if the images are from CASLEO, the 'ctio' script
if the image frames come from the CTIO telescope (CTIO Yale 1-meter Telescope +
Y4KCam), or 'campanas' if they come from 'Las Campanas' observatory.

1-) CTIO ; This script is executed by typing 'ctio' in IRAF

      1- REQUIREMENTS TO RUN THIS SCRIPT
      
      - Working Folder: where the selected CTIO image frames are (*.fits)
      
      * The '*.fits' frames must ALL belong to the SAME FIELD. This is because
      of the last task performed by the script: the aligning. You can run this
      script on several fields, at the same time, but you'll have to skip the
      aligning part: since the fields are not all the same, the task used to
      perform the aligning won't work because it needs to identify the SAME set
      of stars in all the frames to align them. The other option is to run this
      script once for each stars frame set (recommended).
      * The 'trim section' must be known. The default is: [26:4052,26:4052], but
      this value could change.
      
      2- WHAT IT DOES

      * Renames '*.fits' files to '*.fit'
      * Runs 'setairmass' task on all files.
      * Adds GAIN (1.44) and RDNOISE (7) values to the image headers.
      * Exctracts 'FILTER' value (2,3,5 or 6) and writes correct 'FILTERS'
      values (B,V,I or U)
      * Exctracts 'TIME-OBS' value and writes 'UT' value.
      * Trims image's borders (default value: [26:4052,26:4052]) using 'imcopy'
      task.
      * Transposes images two times by [-*,*]
      * Aligns images.

      3- FILES CREATED, COPIED, MOVED OR DELETED
      
      * No files deleted or created, only the '*.fits' files are renamed '*.fit'
      and edited.


1-) CAMPANAS ; This script is executed by typing 'campanas' in IRAF

      1- REQUIREMENTS TO RUN THIS SCRIPT
      
      - Working Folder: where the selected 'Las Campanas' image frames are 
      (*.fits)
      
      * The '*.fits' frames must ALL belong to the SAME FIELD. This is because
      of the last task performed by the script: the aligning. You can run this
      script on several fields, at the same time, but you'll have to skip the
      aligning part: since the fields are not all the same, the task used to
      perform the aligning won't work because it needs to identify the SAME set
      of stars in all the frames to align them. The other option is to run this
      script once for each stars frame set (recommended).
      
      2- WHAT IT DOES
      
      * Renames '*.fits' files to '*.fit'
      * Adds RDNOISE (6.6) values to the image headers.
      * Exctracts 'FILTER' value and writes 'FILTERS' values (B,V,I or U)
      * Exctracts 'UTSTART' value and writes 'UT' value.
      * Transposes images two times by [*,-*] and [*,*]
      * Aligns images.

      3- FILES CREATED, COPIED, MOVED OR DELETED
      
      * No files deleted or created, only the '*.fits' files are renamed '*.fit'
      and edited.
      

2-) CASLEORED ; This script is executed by typing 'casleored' in IRAF

      1- REQUIREMENTS TO RUN THIS SCRIPT
      
      - Working Folder: where ALL image frames are (flats, bias and star frames)
      
      * The names of the files must have the following format:
        -Flat files: "flatwx.fit" or "skywx.fit", where 'w' stands for the
        correct FILTER identification and 'x' stands for the flat number.
        -Bias files: "biasx.fit" where 'x' stands for the bias number.
        -Star frames: "abcdwx.fit" where 'abcd' stands for the unique cluster
        identification, 'w' stands for the correct FILTER identification and 'x'
        stands for the time exposure.
      * GAIN, RDNOISE, TRIMSEC and BIASSEC values must be known since they will
      be asked for.
      * Either sky or dome flats must exist. If both types of files are found,
      script will not work.
      * The instrument MUST be 'direct' (i.e.: setinstrument = direct).
      * The task 'setairmass' must be used on ALL the *.fit files.
      
      2- WHAT IT DOES
      
      * Adds TRIMSEC, BIASSES, GAIN, RDNOISE and FILTERS values to the image
      headers of all the files.
      * Runs 'rfits' on ALL files to set datatype="ushort" (pixel value max:
      64000).
      * Runs 'zerocombine'.
      * Runs 'ccdproc' to apply 'biassec', 'trimsec' and 'zerocor'.
      * Runs 'flatcombine'.
      * Runs 'ccdproc' to apply 'flatcor'.
      * Runs 'imtranspose' to rotate the stars frames.

      3- FILES CREATED, COPIED, MOVED OR DELETED
      
      * "gain_rdnoise" - This file will be used by the 'psfphot' script later
      on.
      * Moves bias files to 'bias/' folder.
      * Moves flats files to 'flats/' folder.
      * Moves combined flats and bias to 'calib/' folder.
      * "imagenes" - This file holds the names of the files correctly modified
      when adding FILTER data to standards and clusters frame's image headers 
      (just for control, not used by any script).
  
        
3-) ALIGN ; This script is executed by typing 'align' in IRAF

      1- REQUIREMENTS TO RUN THIS SCRIPT
      
      - Working Folder: where star images are (both standards and cluster image
      frames); since Flats and Bias were moved to separate folders by the
      previous script.
            
      * The names of the files must have the same format as in CASLEORED.
      * ds9 must be opened.
      
      2- WHAT IT DOES
      
      * Identifies the files common to a single cluster or standard stars set.
      * Runs 'daofind' and interactive 'tvmark' on a user selected frame (the
      resultant ".coo.1" file will be used by the following script to perform
      the actual alignment of the frames).

      3- FILES CREATED, COPIED, MOVED OR DELETED
      
      * One ".coo.1" file per cluster or standard stars set.
      * One "cluster.abcd" file per cluster or standard stars set ('abcd' stands
      for the unique cluster identification).
      * One "frame.abcd" file per cluster or standard stars set.
      * "clusters" - This file stores the following names: "cluster_abcd".

      ALL THIS FILES WILL BE USED (AND THEN DELETED) BY THE FOLLOWING SCRIPT

      
4-) ALIGN2 ; This script is executed by typing 'align2' in IRAF

      1- REQUIREMENTS TO RUN THIS SCRIPT
      
      - Working Folder: same as above
            
      * The file 'gain_rdnoise' must exist (ie: 'casleored' must have been
      used). Else these values will be asked for and the file will be created.
      * All the files created by the previous script must exist.
      
      2- WHAT IT DOES
      
      * Uses the files created by the previous script ('align') to perform the
      actual frame alignment (for ALL the frames, both standard stars and
      clusters).
      * Image files are moved to individual unique folders (standard stars sets
      are grouped in a 'standard/' folder).

      3- FILES CREATED, COPIED, MOVED OR DELETED
      
      * All files created by the previous script are deleted.
      * Image files are moved to a specific folder with the unique name of the
      cluster (ie, the name of the folder is: 'abcd') or to a 'standard/' folder
      for standard stars sets. 
      * A copy of the file 'gain_rdnoise' is stored in each of the created
      folders (except the 'standard/' folder).
      * "cluster_folders" - This file holds the names of the cluster folders
      created (used by 'casleored2' and 'casleored4'). Is stored inside the
      'standard/' folder.



## Miscellaneous


11-) CLEAN ; This script is executed by typing 'clean' in IRAF (this script is
implemented inside 'daom' and will become obsolete)

      1- REQUIREMENTS TO RUN THIS SCRIPT
      
      - Working Folder: this script is to be run wherever the 'night_obs.out.BV'
      file is.
          
      * It needs the file: 'night_obs.out.BV'
      
      2- WHAT IT DOES
      
      * Converts the 'night_obs.out.BV' file to the format needed to feed the
      'shift' code (it also removes stars with INDEFS in any magnitude)
      
      3- FILES CREATED, COPIED, MOVED OR DELETED
      
      * The 'shift_entrada.BV' file is created.


12-) SELECT ; This script is executed by typing 'select' in IRAF (this script
has been implemented into 'psfphot')

      1- REQUIREMENTS TO RUN THIS SCRIPT
      
      - Working Folder: this script is to be run wherever the '*.coo.1' file is.
          
      * It needs the file: '*.coo.1'
      
      2- WHAT IT DOES
      
      * Removes all objects inside a rectangular section from the '*.coo.1'
      file. The limits of this section are:
        xmin (minimum value of x that defines the rectangle to be removed)
        xmax (maximum value of x that defines the rectangle to be removed)
        ymin (minimum value of y that defines the rectangle to be removed)
        ymax (maximum value of y that defines the rectangle to be removed)
        
        
             _______________________[xmax,ymax]
            |                       |
            |                       |
            |                       |
            |_______________________| 
        [xmin,ymin]
      
      
      3- FILES CREATED, COPIED, MOVED OR DELETED
      
      * The '*.coo.1' file is renamed: '*.coo.1_original'
      * A new file with the objects inside the rectangle removed is created with
      the name: '*.coo.1_1'
      * Subsequent uses of this script will append '_original' to the name of
      the input file and create a new one with a '_1' appended to the name of
      the input file.


13-) AVERAGE ; This script is executed by typing 'average' in IRAF

      1- REQUIREMENTS TO RUN THIS SCRIPT
      
      - Working Folder: this script is to be run wherever the '*_data' file are.
          
      * It needs the files: '*_data'
      
      2- WHAT IT DOES
      
      * Obtains an average of FWHM, Sky Mean and STDDEV values stored in the
      '*_data' files.
      * This script is to be used after the 'getdata' script was used over
      standard frames, to obtain the average values of FWHM and STDDEV used in
      the 'casleored2' script.
      
      
      3- FILES CREATED, COPIED, MOVED OR DELETED
      
      * One 'average' file is created.

