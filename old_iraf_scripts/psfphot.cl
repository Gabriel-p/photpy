################################################################################
#
#        ===========================================================
#        PROCEDURE TO PERFORM POINT SPREAD FUNCTION (PSF) PHOTOMETRY
#        ===========================================================
# 
#     To run properly this procedure must be executed inside the package:
#                           'noao.digiphot.daophot'
#
#          This script follows the 'A user's guide to stellar CCD photometry
#                     with IRAF' by Massey & Davis (1992).
#
#
#                            by Gabriel Perren 2009-2012
################################################################################


procedure psfphot (image, datamax)

file image {prompt = "Input image name"}
#real sky_mean {prompt = "Sky mean value [MEAN]"}
#real sky_std {prompt = "Sky standard desviation value [STDDEV]"}
#real fwhm {prompt = "Typical FWHM for stellar image"}
#real psf_rad {prompt = "(4xFWHM+1 aprox) Radius of bigger PSF star in the frame"}
#real aperture1 {prompt = "(FWHM aprox) Aperture photometry radius [pixels]"}
real datamax {prompt = "Maximum good data value?"}
bool interactive {prompt = "Run under interactive mode?"}
struct *file_var {mode="h", prompt = "Internal file name variable"}
struct *file_var2 {mode="h", prompt = "Internal file name variable"}
struct *flist {mode="h", prompt = "Internal file name variable"}


begin
      string imname, imname2
      bool filelist_check, check, dopsf, psf_done, no_neighbours, right_sustraction
      bool first_psf, first_nstarsubstar, use_psf, star_search, auxi, newstars, doapert
      bool sigma_warning, apert_warning, nopsfstar, pstselect, thresh_change, tvmark, pstauto
      bool core_psf, rem_psf, better_psf, better_psf2, do_pstselect, apert_null, dao_perf
      real dmax, fitrad, psfrad, nstarpsfrad, sigma, gain, rdnoise, sstd, value
      real fannulus, sannulu, fdannulu, datamin, sigma_old
      real apstd, aperture
      real smean, thresh
      real min, max, xmin, xmax, ymin, ymax
      real avg, mag1, mag2, delta_m, maxnpsf, psfnumber
      string m1, m2, null, expression
      real diff, oldfitrad, oldthresh, oldpsfrad, olddmax
      int k, i, l, n, m, o, p, var1, var2
      struct line, line2

      imname = image
      dmax = datamax
      inter=interactive

# ------------------------------------------------------------------------------------
# Package loaded control

	if (! defpac ("daophot")) {
	  print ('\n This script must be loaded inside the package noao/digiphot/daophot')
	  bye()
	}
      else {
      }
# ------------------------------------------------------------------------------------   

#     delete ((imname//'.parameters'),verify=no,>>&"/dev/null")

# ------------------------------------------------------------------------------------   
# Check image name (if input file has an '@', treat as a list of files iterating through the script)
# ------------------------------------------------------------------------------------ 

      filelist_check = no
      if (substr (imname, 1, 1) == "@") {
        filelist_check = yes
        k = strlen(imname)
	  imname2 = substr (imname, 2, k)  
        cp (imname2, 'fitfileslist')
        file_var2 = 'fitfileslist'
	  while (fscan (file_var2,line2) != EOF) {
	    imname = line2
	    goto runscript                                 # THESE GOTO STATEMENTS ACCOUNT FOR THE FACT THAT THE INPUT
	    print ('')                                     # FILE MIGHT BE A LIST OF FILES, ie: @filename
	    rerun:                                         # IF SO, THEN RUN THE SCRIPT FOR EACH FILE AND THEN EXIT
	    print ('\n New iteration with file: '//imname) # ELSE, RUN THE SCRIPT FOT THE SINGLE FILE INPUT.           
	    print ('')
	  }      
      }
      else {
        print ('\n Not a list of files: single file')
        print (' Running script')
      } 
         
	if (filelist_check) {
	  delete ('fitfileslist')
	  if (sigma_warning == no && apert_warning == no) {
	    print ('\n ----------------------------------------------------- ')
	    print (' Script "psfphot" finished correctly.                ')
	  }
	else if (sigma_warning == no && apert_warning == yes) {
	  print ('\n ----------------------------------------------------- ')
	  print (' Script "psfphot" finished with a WARNING:           ')
	  print ('\n WARNING: Not enough stars to compute an aperture correction')
	}
	else if (sigma_warning == yes && apert_warning == no) {
	  print ('\n ----------------------------------------------------- ')
	  print (' Script "psfphot" finished with a WARNING:           ')
	  print ('\n WARNING: sigma and sstd [SSTDEV] values differed significantly')
	  print (' sigma = '//sigma_old//' ; sstd = '//sstd)
	}
	else {
	  print ('\n ----------------------------------------------------- ')
	  print (' Script "psfphot" finished with two WARNINGS:        ')
	  print ('\n WARNING: sigma and sstd [SSTDEV] values differed significantly')
	  print (' sigma = '//sigma_old//' ; sstd = '//sstd)
	  print ('\n WARNING: Not enough stars to compute an aperture correction')		          
	}

	print ('\n After processing ALL THE CLUSTERs FRAMES:')
	print ('\n Move on to the next script: "daom"')
	print ('\n Remember this last script must be executed inside the package:')
	print ('\n             noao/digiphot/photcal                     ')
	print (' ----------------------------------------------------- ') 
	bye      
	}

	runscript:

      print ('\n###########################################', >> (imname//'.parameters'))	
      date >> (imname//'.parameters')
      print ('\nPSF Phot script has started', >> (imname//'.parameters'))	

      k = strlen(imname)
      if (substr (imname, k-4, k) == ".fits") {
        print ('\n Image file name: '//imname)
        print (' Is this the image name PLUS the \'.fits\' extension? (y/n)')
        print (' default = yes')
        check = yes
        scan (check)
        if (check) {
          imname = substr (imname, 1, k-5)
          print ('\n Renaming file from '//imname//'.fits to '//imname//'.fit')
          copy (imname//'.fits', imname//'.fit')
          delete (imname//'.fits')
          print ('\n\'.fits\' file renamed to \'.fit\' file', >> (imname//'.parameters'))		                  
          print ('\n Image name: '//imname)
          print ('Image name: ', imname, >> (imname//'.parameters'))
        }
        else {
          print ('\n Unknown \'image name\' error.')
          print (' Closing.')
          bye()
        }
      }
      else {
        if (substr (imname, k-3, k) == ".fit") {
          print ('\n Image file name: '//imname)
          print (' Is this the image name PLUS the \'.fit\' extension? (y/n)')
          print (' default = yes')
          check = yes
          scan (check)
          if (check) {
            imname = substr (imname, 1, k-4)
            print ('\n Image name: '//imname)
            print ('nImage name: ', imname, >> (imname//'.parameters'))
          }
          else {
            print ('\n Unknown \'image name\' error.')
            print (' Closing.')
            bye()
          } 
        }
        else {
          print ('\n Image name: '//imname)
          print ('\nImage name: ', imname, >> (imname//'.parameters'))
        }
      }
      
# ------------------------------------------------------------------------------------   
# End of 'Check image name'
# ------------------------------------------------------------------------------------ 
   

# ------------------------------------------------------------------------------------ 
# Input 'GAIN', 'RDNOISE' and 'apert_standard' values
# ------------------------------------------------------------------------------------ 

      # Acquire the values of "gain" and "rdnoise" from the file 'gain_rdnoise'
      # previously created by "casleored" script. If such file does not exists, input the values.
      #
      files ('*gain_rdnoise', > 'temp.grd')
      file_var = 'temp.grd'
      n=0
      while (fscan (file_var,line) != EOF) {
        n = n + 1
      }
      if (n == 0) {
        print ('\n No "gain_rdnoise" file.')
          
	  hselect.mode = "hl"
	  hselect.images = imname
	  hselect.fields = "GAIN"
	  hselect.expr = yes
	  hselect > "tempget"  # Get value from header
	  file_var = 'tempget'
	  while (fscan (file_var,gain) != EOF)
	  del ("tempget")

	  hselect.mode = "hl"
	  hselect.images = imname
	  hselect.fields = "RDNOISE"
	  hselect.expr = yes
	  hselect > "tempget"  # Get value from header
	  file_var = 'tempget'
	  while (fscan (file_var,rdnoise) != EOF)
	  del ("tempget")
          
        print ('\n GAIN = '//gain)
        print (' RDNOISE = '//rdnoise)
        print ('\n Accept GAIN and RDNOISE values read from header? (y/n)')
        print (' default = yes')
        check=yes
        scan (check)
        if (check) { # Do nothing
        }
        else {
	    print ('\n Input GAIN value')
	    scan (gain)          
	    print (' Input RDNOISE value')
	    scan (rdnoise)
        }
          
        print ('Gain = ', gain, >> (imname//'.parameters'))
        print ('Rdnoise = ', rdnoise, >> (imname//'.parameters'))
        print (gain, rdnoise, >> 'gain_rdnoise')          
          
      }
      else {
         file_var = 'gain_rdnoise'
         while (fscan (file_var, gain, rdnoise) != EOF)
         print ('Gain = ', gain, >> (imname//'.parameters'))
         print ('Rdnoise = ', rdnoise, >> (imname//'.parameters'))      
      }
      del ('temp.grd')

      # Acquire of the "aperture" value used for the standards stars from the file 'apert_standard'
      # previously created by "casleored2" script. If such file does not exists, input the value.
      #
      files ('*apert_standard', > 'temp.apert')
      file_var = 'temp.apert'
      m=0
      while (fscan (file_var,line) != EOF) {
          m = m + 1
      }
      if (m == 0) {
        print ('\n No "apert_standard" file.')
        print ('\n Input aperture used for standard stars [pixels]')
        scan (apstd)
        print ('Inputed aperture used for standard stars [pixels] = ', apstd, >> (imname//'.parameters')) 
        print (apstd, >> 'apert_standard')
      }
      else {      
        file_var = 'apert_standard'
        while (fscan (file_var, apstd) != EOF)
        print ('Aperture used for standard stars [pixels] (from file) = ', apstd, >> (imname//'.parameters')) 
      }
      del ('temp.apert')
      
# ------------------------------------------------------------------------------------ 
# End of 'Input 'GAIN', 'RDNOISE' and 'apert_standard' values'
# ------------------------------------------------------------------------------------


# ------------------------------------------------------------------------------------ 
# Get mean FWHM, SKY MEAN and SKY STANDARD DEVIATION value from several stars
# ------------------------------------------------------------------------------------

      print ('\n Read values from file obtained with \'getdata\' script')
      print (' (else input values manually)? (y/n)')
      print (' default = yes')
      check = yes
      scan (check)
      if (check) {
	  # ------------------------------------------------------------------------------------
	  # Search for 'imname//'_data'' file
	  #
	  auxi=no
	  while (auxi==no) { 
	    o=0
	    # Search for "imname//'_data'"
	    # ------------------------------------------------------------------------------------
	    files ('*'//imname//'_data', > 'tempfile')
	    flist = 'tempfile'
	    while (fscan (flist,line) != EOF) {
	      o = o + 1
	    }
	    del ("tempfile")
	    # ------------------------------------------------------------------------------------

	    if (o!=1) {
	      print('')
	      print (' ****************************')
	      print ('   Missing "'//imname//'_data" file.')
	      print (' ****************************')
	      print ('\n Continue and input values manually (else exit)? (y/n)')
	      scan (check)
	      if (check) {
	        print ('\n Input FWHM, SSTDEV and SKY MEAN values')
	        print ('\n Input FWHM value')
	        scan (fitrad)
	        print ('\n Input SSTDEV value')
	        scan (sstd)
	        print ('\n Input SKY MEAN value')
	        scan (smean)
	        auxi = yes
	      }
	      else {
	        bye()
	      }
	    } # closes the 'if o!=1'
	    else {
	      if (o==1) {
	        print ('')
	        print ('"'//imname//'_data" file found.')
	        auxi=yes
	        i = 1
	        file_var = imname//'_data'
	        while (fscan (file_var, value, null) != EOF) {
	          if (i == 1) {
	            fitrad = value
	            print ('\n FWHM = '//fitrad)
	          }
	          else {
	            if (i == 2) {
	              sstd = value
	              print ('\n SSTDEV = '//sstd)
	            }
	            else {
	              smean = value
	              print ('\n Sky Mean = '//smean)
	            }
	          } 
	          i = i+1   
	        }
	        print ('\n Accept values read from file (else input values manually)? (y/n)')
	        print (' default = yes')
	        check = yes
	        scan (check)
	        if (check) {
	        }
	        else {
	          print ('\n Input FWHM, SSTDEV and SKY MEAN values')
	          print ('\n Input FWHM value')
	          scan (fitrad)
	          print ('\n Input SSTDEV value')
	          scan (sstd)
	          print ('\n Input SKY MEAN value')
	          scan (smean)         
	        }
	      } # closes the 'if o=1'
	      else {
	        print (' Unknown error. Check code: line 360')
	        bye()
	      }
	    } # closes the 'else' before 'if o=1'

	  } # This bracket closes the 'auxi' 'while'

	} # This bracket closes the 'if' that checks for the existence of the "'//imname//'_data" file
      else {
        print ('\n Input FWHM, SSTDEV and SKY MEAN values')
        print ('\n Input FWHM value')
        scan (fitrad)
        print ('\n Input SSTDEV value')
        scan (sstd)
        print ('\n Input SKY MEAN value')
        scan (smean)          
      }  
     
# ------------------------------------------------------------------------------------ 
# End of 'Get mean FWHM, SKY MEAN and SKY STANDARD DEVIATION value from several stars'
# ------------------------------------------------------------------------------------


# ------------------------------------------------------------------------------------ 
# Input 'psfrad' and 'aperture' values
# ------------------------------------------------------------------------------------

      psfrad = 4*fitrad+1  # According to 'A Reference Guide to the IRAF-DAOPHOT Package'
      aperture = fitrad    # by L. Davis (page 31)
      if (inter) {
	  print ('\n Accept radius of bigger PSF star in the frame')
	  print (' (\'psfrad\') value as: '//psfrad//' (4xFWHM + 1)')
	  print (' (else input another \'psfrad\' value)? (y/n)')
	  print (' default = yes')
	  check = yes
	  scan (check)
	  if (check) { # Do nothing
	  }
	  else {
	    print ('\n Input new \'psfrad\' value:')
	    scan (psfrad)
	  }
	  print ('\n Accept aperture photometry radius [pixels]')
	  print (' (\'aperture\') value as: '//aperture//' (0+FWHM)')
	  print (' (else input another \'aperture\' value)? (y/n)')
	  print (' default = yes')
	  check = yes     
	  scan (check)
	  if (check) { # Do nothing
	  }
	  else {
	    print ('\n Input new \'aperture\' value:')
	    scan (aperture)
	    print (' New aperture value: aperture = ', aperture, >> (imname//'.parameters'))
	  }
      }

# ------------------------------------------------------------------------------------ 
# End of 'Input 'psfrad' and 'aperture' values'
# ------------------------------------------------------------------------------------

 
# ------------------------------------------------------------------------------------------- 
# Input parameters
# ------------------------------------------------------------------------------------------- 

      print ('\nParameters used during the PSF procedure', >> (imname//'.parameters'))
      print ('', >> (imname//'.parameters'))
      print ('Image name = ', imname, >> (imname//'.parameters'))
      print ('Sky mean value = ', smean, >> (imname//'.parameters'))
      print ('Sky standard desviation value [SSTDEV] = ', sstd, >> (imname//'.parameters'))
      print ('Maximum good data value = ', dmax, >> (imname//'.parameters'))
      print ('Typical FWHM for stellar images = ', fitrad, >> (imname//'.parameters'))
      print ('Radius of bigger PSF star in the frame = ', psfrad, >> (imname//'.parameters'))
      print ('Aperture photometry radius [pixels] = ', aperture, >> (imname//'.parameters'))

# ------------------------------------------------------------------------------------------- 
# End of 'Input parameters'
# -------------------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Parameters input, Daofind and Phot tasks
# ------------------------------------------------------------------------------    

      if (smean*gain + rdnoise*rdnoise < 0) {
        print ('\n smean*gain + rdnoise*rdnoise < 0')
        print (' Using sigma = 0.')
        sigma_old = 0
        sigma = 0
        print ('\nWARNING: smean*gain + rdnoise*rdnoise < 0', >> (imname//'.parameters'))
        print ('Using sigma=0', >> (imname//'.parameters'))
        print ('\nWARNING: smean*gain + rdnoise*rdnoise')
        print ('Using sigma=0')
      }
      else {
        sigma = sqrt(smean*gain + rdnoise*rdnoise)/gain
      }

      sigma_warning = no
      diff = sigma - sstd
      if ((diff <= -sstd*0.25) || (diff >= sstd*0.25)) { # A difference of less than abs(25%) is acceptable (totally arbitrary criteria)
                                                         # Page 27 of 'A Reference Guide to the IRAF-DAOPHOT Package'
                                                         # by L. Davis, recommends to check for this difference.
        print ('\n WARNING: sigma and sstd [SSTDEV] values differ significantly')
        print (' sigma = '//sigma//' ; sstd = '//sstd)
        print (' difference = '//diff)
        print (' ', >> (imname//'.parameters'))
        print (' WARNING: sigma and sstd [SSTDEV] values differ significantly = ', diff, >> (imname//'.parameters'))
        print (' sigma = '//sigma//' ; sstd = '//sstd, >> (imname//'.parameters'))
        print (' ', >> (imname//'.parameters'))
        sigma_warning = yes
        if (inter) {
	    print ('\n Continue (else exit)? (y/n)')
	    print (' default = yes')
	    check = yes
	    scan (check)
	    if (check) { # Do nothing
	    }
	    else {
	      bye
	    }
        }
      }

      unlearn datapars
      unlearn daopars
      print (' ----------------------------------------------------- ')
      print (' Datapars and daopars tasks                            ')
      print (' Basic image parameters                                ')

      if (sigma == 0.) {
        sigma = sstd
        print ('Using SSTDEV to calculate the Threshold', >> (imname//'.parameters'))
        print ('\n Using SSTDEV to calculate the Threshold (since sigma = 0)')
      }
      else {
	  print ('\n Use \'sigma\'= '//sigma//' value to calculate the threshold as:')
	  print (' threshold = 3.5*sigma = '//3.5*sigma//' or else use')
	  print (' \'SSTDEV\'= '//sstd//' value:')
	  print (' threshold = 3.5*STDDEV = '//3.5*sstd)
	  print ('\n Input y/n (yes = sigma ; no = STDDEV)')
	  check=yes
	  scan (check)
	  if (check) { # Do nothing
	  }
	  else {
	    sigma = sstd
	    print ('Using STDDEV to calculate the Threshold', >> (imname//'.parameters'))          
	  }
      }

      findpars.threshold = 3.5*sigma
      print (' Threshold = '//3.5*sigma)
      print ('Threshold = ', (3.5*sigma), >> (imname//'.parameters')) 
      datapars.emission = yes
      datapars.fwhmpsf = fitrad
      datapars.sigma = sstd
      # According to a 'A user's guide to stellar CCD photometry with IRAF' by Massey & Davis
      # (page 41): '...set datamin in datapars to some value like (sky value) - 3*sigma...'
      # According to IRAF's datapar help: '...DAOPHOT users should either leave datamin set to
      # INDEF or set it to a number between 5-7 sigma below the sky background value'
      #
      # So we set this value to smean-4.5*sigma
      print ('\n Datamin = smean-4.5*sigma = '//smean-4.5*sigma)
      print ('\n Accept datamin value? (y/n)')
      print (' default = yes')
      check=yes
      scan (check)
      if (check) {
	   print ('Datamin = ', (smean-4.5*sigma), >> (imname//'.parameters')) 
	   datamin = smean-4.5*sigma
	   datapars.datamin = datamin
      }
      else {
        print ('\n Input new datamin value:')
        scan (datamin)
        print ('Datamin = ', datamin, >> (imname//'.parameters')) 
        datapars.datamin = datamin
      }

      datapars.datamax = dmax
      daopars.psfrad = psfrad
      daopars.fitrad = fitrad # = FWHM
      print (' Minimum data value = ' // datamin                )
      print (' Maximum data value = ' // dmax                         )
      print (' Fitrad = ' // fitrad // ' pixels                      ')
      print (' Psfrad = ' // psfrad // ' pixels                      ')      
      print (' ----------------------------------------------------- ')
      
      if (inter) {
	  print ('\n Accept the following data to be inputed in \'datapars\' ')
	  print (' (else don\'t overwrite these values and use values already')
	  print (' stored in \'datapars\' task)? (y/n)')
	  print (' default = yes')
	  print ('\n datapars.noise = "poisson"')
	  print (' datapars.ccdread = "RDNOISE"')
	  print (' datapars.gain = "GAIN"')
	  print (' datapars.exposur = "EXPTIME"')
	  print (' datapars.filter = "FILTERS"')
	  check = yes
	  scan (check)
	  if (check) {
	    datapars.noise = "poisson"
	    datapars.ccdread = "RDNOISE"
	    datapars.gain = "GAIN"
	    datapars.exposur = "EXPTIME"
	    datapars.filter = "FILTERS"		      
	  }
	  else {
	    print ('\n Use stored values in \'datapars\' task (else exit, change')
	    print (' these values, and re-run the script)? (y/n)')
	    print (' default = yes')
	    check = yes
	    scan (check)
	    if (check) { # Do nothing
	    }
	    else {
	      bye()
	    }
	  }
      }
      else {
	  datapars.noise = "poisson"
	  datapars.ccdread = "RDNOISE"
	  datapars.gain = "GAIN"
	  datapars.exposur = "EXPTIME"
	  datapars.filter = "FILTERS"	      
      }

	psf_done=no
	while (psf_done == no) { # THIS WHILE CLOSES AFTER THE CHECKING OF THE
	                         # .coo.1, .als.1 and .als.2 (if it exists) FILES
      star_search = no
      thresh = 3.5*sigma
      thresh_change = no
      while (star_search == no) {

      dao_perf=yes
      while (dao_perf) {
          
	  if (thresh_change == no) {
	    print ('\n Perform NEW Daofind (else use EXISTING '//imname//'.coo.1 file)? (y/n)')
	    print (' (if answer = yes, this will DELETE any \'.coo.1\' file)')
	    print (' default = no')
	    check = no
	    scan (check)
	  }
	  else {
	    check = yes
	  }    
	  if (check) {
	    o=0
	    # Search for "imname.coo.1"
	    # ------------------------------------------------------------------------------------
	    files ('*'//imname//'.coo.1', > 'tempfile')
	    flist = 'tempfile'
	    while (fscan (flist,line) != EOF) {
	      o = o + 1
	    }
	    del ("tempfile")
	    # ------------------------------------------------------------------------------------

	    if (o!=1) {
	      dao_perf=no # File doesn't exist. Move on.
	    }
	    else {
	      if (o==1) {
	        print ('\n DELETE '//imname//'.coo.1 file? (y/n)')
	        scan (check)
	        if (check) {
	          delete (imname//'.coo.1',verify=no,>>&"/dev/null")
	          dao_perf=no
	        }
	        else {
	          dao_perf=yes
	        }
	      }
	      else {
	        print (' Unknown error. Try again.')
	        dao_perf=yes
	      }
	    }
	  }
	  else {
	    dao_perf=no
	    # ------------------------------------------------------------------------------------
	    # Search for 'imname.coo.1' file
	    # ------------------------------------------------------------------------------------
	    auxi=no
	    while (auxi==no) { 
	      o=0
	      # Search for "imname.coo.1"
	      # ------------------------------------------------------------------------------------
	      files ('*'//imname//'.coo.1', > 'tempfile')
	      flist = 'tempfile'
	      while (fscan (flist,line) != EOF) {
	        o = o + 1
	      }
	      del ("tempfile")
	      # ------------------------------------------------------------------------------------

	      if (o!=1) {
	        print ('\n ***********************************')
	        print ('   Missing '//imname//'.coo.1 file.')
	        print (' ***********************************')
	        print ('\n Search for '//imname//'.coo.1 file again (else run')
	        print (' \'Daofind\' task to create one)? (y/n)')
	        print (' default = yes')
	        check = yes
	        scan (check)
	        if (check) {
	          auxi=no
	        }
	        else {
	          auxi=yes
	          check = yes # Move on to 'Daofind' task
	        }
	      }
	      else {
	        if (o==1) {
	          print ('\n '//imname//'.coo.1 file found.')
	          auxi=yes
	          check = no # Do not perform 'Daofind' task
	        }
	        else {
	          print (' Unknown error. Try again.')
	          dao_perf=yes
	        }
	      }
	    }
	    # ------------------------------------------------------------------------------------
	    # End of 'Search for 'imname.coo.1' file'
	    # ------------------------------------------------------------------------------------ 		          
	  }
	} # This bracket closes the 'dao_perf' 'while'

	if (check) { 
        print ('------------', >> (imname//'.parameters')) 
        print ('Daofind task', >> (imname//'.parameters'))
        print ('------------', >> (imname//'.parameters'))
	  print ('\n ----------------------------------------------------- ')
	  print (' Daofind task:                                         ')
	  print (' Search for stars                                      ')
	  print (' Threshold value: ' // thresh)
	  print (' ----------------------------------------------------- ')
	  daofind.verif = no
	  daofind.verb = yes
	  daofind.interactive = no
	  daofind ((imname), (imname//'.coo.1'))		
	}                  

	file_var = (imname//'.coo.1')
	m=0
	while (fscan (file_var,line) != EOF) {
	  m = m + 1
	}
	print ('\n Number of stars found = '//(m-41))
      print ('Number of stars found = ', (m-41), >> (imname//'.parameters'))	      
	print (' Threshold = '//thresh)
	      
      print ('\n Star search done.')
      print ('\n Press \'y\' key to display stars found.')
      print ('\n Else press \'n\' key to DELETE the '//imname//'.coo.1 file,')
      print (' select a new \'threshold\' value and perform a new \'Daofind\' run.')
      print (' default = yes')
      check = yes
      scan (check)
      if (check) {
          
        tvmark=yes
        while (tvmark==yes) {
	    print ('\n ----------------------------------------------------- ')
	    print (' Display task:                                         ')
	    print (' Original Image: ' // imname                            )
	    print (' ----------------------------------------------------- ')
	    display ((imname), 1)
	      
	    print ('\n ----------------------------------------------------- ')
	    print (' Tvmark                                                ')
	    print (' -red point = detected star                            ')
	    print (' Add or remove stars in necessary                      ')
	    print ('\n -the "a" key to add an undetected star              ')
	    print (' -the "d" key to remove a detected star                ')
	    print (' -the "q" key to exit interactive tvmark               ')
	    print (' ----------------------------------------------------- ')
	    tvmark.interactive = no
	    tvmark.mark = 'point'
	    tvmark.font = "raster"
	    tvmark.txsize = 2
	    tvmark.color = 204
	    tvmark.number = no
	    tvmark.label = no
	    tvmark (1, (imname // '.coo.1'))
	      
	          # The first tvmark marks the stars in the image, found by 'daofind'
	          # The second one (below this), performs the 'interactive' tvmark.
	          # If I try to make the first tvmark 'interactive', then
	          # it doesn't mark the stars, this way the found star are marked

          tvmark.interactive = yes 
	    tvmark.number = yes
	    tvmark.toleran = 4
	    tvmark (1, (imname // '.coo.1'))   

	    auxi=no
	    while (auxi == no) {

          print ('\n Now, do you want to...:')
          print ('\n     [1]: Move on to \'phot\' task')
          print ('\n     [2]: Repeat \'tvmark\' over this \'.coo.1\' file')
          print ('\n     [3]: Remove all objects inside a rectangular section and')
          print ('          then re-load \'tvmark\' ')
          print ('\n     [4]: Repeat \'daofind\' with a different \'threshold\'')
          print ('          (this means DELETING the actual \'.coo.1\' file)')
          print ('\n     [5]: Save actual \'.coo.1\' file and exit scritp.')

		                  var2=2
	              scan (var2)
	              if (var2==1) {
	                  auxi=yes
	                  tvmark=no
		                      star_search = yes
	                  print('\n Moving on to \'phot\' task.')
	              }
	              else {
	                  if (var2==2) {
		                          auxi=yes
		                          tvmark=yes
		                          print('\n Repeating \'tvmark\' over \'.coo.1\' file.')
	                  }
	                  else {
	                      if (var2==3) {
	                         
	                          print ('\n Input xmin, ymin, xmax and ymax values like shown below')
	                          print ('')
	        print ('      _______________________[xmax,ymax]')
	        print ('     |                       |')
	        print ('     |                       |')
	        print ('     |                       |')
	        print ('     |_______________________|')
	        print (' [xmin,ymin]') 
	                      
	                          print ('\n Input minimum x value')
	                          scan (xmin)
	                          print ('\n Input minimum y value')
	                          scan (ymin)
	                          print ('\n Input maximum x value')
	                          scan (xmax)
	                          print ('\n Input maximum y value')
	                          scan (ymax)

                                  expression = "XCENTER < "//xmin
		      expression = expression // " "// '|| XCENTER > '//xmax
		      expression = expression // " "// '|| YCENTER < '//ymin
		      expression = expression // " "// '|| YCENTER > '//ymax

		      print ('\n Expression: '//expression)

		      rename.field = 'all'
		      rename (imname//'.coo.1', imname//'_original')

		      pselect.infiles = imname//'_original'
		      pselect.outfiles = imname//'.coo.1'
		      pselect.expr = expression
		      pselect.mode = "hl"
		      pselect
		      delete (imname//'_original')
		      
	                          auxi=yes
	                          tvmark=yes
	                          print('\n Repeating \'tvmark\' over \'.coo.1\' file.')		      
	                      }
	                      else {
	                          if (var2==4) {
              print ('\n DELETE '//imname//'.coo.1 file, change \'threshold\' value')
		          print (' and repeat \'Daofind\'? (y/n)')
		          print (' default = yes')
		          scan (check)
		          if (check) {
		              auxi=yes
		              tvmark=no
		              thresh_change = yes
		              print ('\n Actual threshold = ' //thresh)
		              print (' Input new threshold value')
		              print ('\n (A higher threshold means less stars)')
		              scan (thresh)
		              findpars.threshold = thresh
		              print ('New Findpars.threshold = ', thresh, >> (imname//'.parameters'))
		              del (imname// '.coo.1')
              }
              else {
	                                  auxi=no
              }	                          
	                          }
                      else {
                          if (var2==5) {
                              print ('\n Halt script')
                              bye()
                          }
		                      else {
		                          print (' Invalid choice. Try again.')
		                          auxi=no
		                      }
                      }
	                      }
	                  }
	              } # This bracket closes the 'if (var2==1)' 'else'
	          } # This bracket closes the 'while (auxi == no)' 'while'

              } # This bracket closes the 'tvmark' 'while'
          } # This bracket closes the 'if (check)' 'if' after 'Star search done'
          else {
              star_search = no
              print ('\n DELETE '//imname//'.coo.1 file, change \'threshold\' value')
		          print (' and repeat \'Daofind\' (else exit script)? (y/n)')
		          print (' default = yes')
		          check=yes
		          scan (check)
		          if (check) {
		              print ('\n Actual threshold = ' //thresh)
		              print (' Input new threshold value')
		              print ('\n (A higher threshold means less stars)')
		              thresh_change = yes
		              scan (thresh)
		              findpars.threshold = thresh
		              print ('New Findpars.threshold = ', thresh, >> (imname//'.parameters'))
		              del (imname// '.coo.1')		              
              }
              else {
                  print ('****************************', >> (imname//'.parameters'))
                  print ('Exited script after Daofind', >> (imname//'.parameters'))
                  bye()
              }
          }		          

      } # This bracket closes the 'star_search' 'while'
      
      print ('\n ----------------------------------------------------- ')
      print (' Phot task (first run)                                 ')
      print (' Aperture = '//aperture                                 )
      print (' ----------------------------------------------------- ')

#      print ('\n Recenter stars marked by tvmark?')
#      print (' (Useful when using .coo. files from other frames)')
#      print (' default = no')
#      check = no
#      scan (check)
#      if (check) {
#      
#          apphot
#          
#		      unlearn centerpars
#		      unlearn fitskypars
#		      unlearn photpars
#		      unlearn psf
#		      fitskypars.salgorithm = "mode" # From Massey-Davis guide to stellar CCD photometry
#		      centerpars.calgorithm = "centroid"
#		      centerpars.cbox = 2.5*fitrad
#		      fannulus = 4*fitrad
#		      fdannulu = 3.25*fitrad      
#		      fitskypars.annulus = fannulus
#		      fitskypars.dannulus = fdannulu      
#		      
#		      phot.interactive = no
#		      phot.radplots = no
#		      phot.update = yes
#		      phot.verbose = yes
#		      phot.verify = no
#		      photpars.apertures = aperture
#		      phot ((imname), "default", "default")
#		      
#		      daophot
#      
#      }
#      else {

      unlearn centerpars
      unlearn fitskypars
      unlearn photpars
      unlearn psf
      fitskypars.salgorithm = "mode" # From Massey-Davis guide to stellar CCD photometry
      centerpars.calgorithm = "none" 
      
      # According to 'A Reference Guide to the IRAF-DAOPHOT Package'
      # by L. Davis (page 31): cbox = 2xFWHM (or 5, wichever is greater)
      #                        annulus = 4xFWHM
      #                        dannulu = 2.5-4.0xFWHM
     
      # According to IRAF help: a reasonable value for 'cbox' is 2.5-4.0 * FWHM
      
      # According to 'A User's Guide to Stellar CCD Photometry with IRAF'
      # by Massey-Davis (page 47): cbox = 5 (approx 2.0-3.0xFWHM)
      #                            annulus = 10 (approx 3.0-4.0xFWHM)
      #                            dannulu = 10 (approx 3.0-4.0xFWHM)
       
      centerpars.cbox = 2.5*fitrad
      fannulus = 4*fitrad
      fdannulu = 3.25*fitrad      
      fitskypars.annulus = fannulus
      fitskypars.dannulus = fdannulu
      print ('\n Fitskypars.annulus = '//fannulus//' (4*fitrad)')
      print (' Fitskypars.dannulus = '//fdannulu//' (3.25*fitrad)')
      if (inter) {
		      print ('\n Accept \'annulus\' and \'dannulu\' values (else')
		      print (' input new values)? (y/n)')
		      print (' default = yes')
		      check = yes
		      scan (check)
		      if (check) {
		          print ('Fitskypars.annulus = ', fannulus, >> (imname//'.parameters')) 
		          print ('Fitskypars.dannulu = ', fdannulu, >> (imname//'.parameters')) 		          
		      }
		      else {
		          print (' Input new \'annulus\' value')
		          scan (fannulus)
		          print (' Input new \'dannulus\' value')
		          scan (fdannulu)		          
		          fitskypars.annulus = fannulus
		          fitskypars.dannulus = fdannulu
		          print ('Fitskypars.annulus = ', fannulus, >> (imname//'.parameters')) 
		          print ('Fitskypars.dannulus = ', fdannulu, >> (imname//'.parameters')) 		          
		      }
      }
      else {
	        print ('Fitskypars.annulus = ', fannulus, >> (imname//'.parameters')) 
	        print ('Fitskypars.dannulus = ', fdannulu, >> (imname//'.parameters'))       
      }   
      phot.interactive = no
      phot.radplots = no
      phot.update = yes
      phot.verbose = yes
      phot.verify = no
      photpars.apertures = aperture
      phot ((imname), "default", "default")
      
#      }
      
      txdump.mode = 'hl' 
      txdump.textfile = (imname//'.mag.1')
      txdump.headers = no
      txdump.fields = 'MSKY,STDEV'
      txdump.expr = 'MAG[1]!=INDEF'
      txdump > auxiliar

      file_var = ('auxiliar')
      m=0
      while (fscan (file_var,line) != EOF) {
          m = m + 1
      }
      
      delete ('auxiliar')      

      print ('')
      print (' Number of stars found (with MAG != INDEF) = '//m)
      maxnpsf = int(m*0.02) # The 'maxnpsf' number is the number of stars the 'pstselect'
      if (maxnpsf < 20) {    # task will look for. We set this number to 2% of the stars
          maxnpsf = 20       # found with MAG!=INDEF (but not less than 20).
      }
		            
# ------------------------------------------------------------------------------
# End of 'Parameters input, Daofind and Phot tasks'
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Perform PSF?
# ------------------------------------------------------------------------------

      if (inter) {
		      print ('\n Perform PSF photometry (else keep Aperture photometry)? (y/n)')
		      print (' default = yes')
		      check = yes
		      scan (check)
		      if (check) {
		          dopsf = yes
		      }
		      else {
		          dopsf = no
		          l = 3 # This parameter will be used when the 'list' file is created. The number 3 indicates
		      }         # that there is no '.als' file but rather a '.mag' file.
      }
      else {
          dopsf = yes
      }
      
# ------------------------------------------------------------------------------
# End of 'Perform PSF?'
# ------------------------------------------------------------------------------


      right_sustraction = no
	while (right_sustraction == no) { # This 'while' encloses the whole 'Core PSF photometry'
                                  # (or 'Aperture photometry') structure
      if (dopsf==yes) {
		      print ('\n ----------------------------------------------------- ')
		      print ('             Psf photometry...')
		      print ('\n ----------------------------------------------------- ') 
		      print ('', >> (imname//'.parameters'))
		      print ('************************************************', >> (imname//'.parameters'))
		      print ('Core PSF photometry', >> (imname//'.parameters')) 
      }
      else {
		      print ('\n ----------------------------------------------------- ')
		      print ('             Aperture photometry...')
		      print ('\n ----------------------------------------------------- ') 
		      print ('-------------------', >> (imname//'.parameters'))
		      print ('Aperture photometry', >> (imname//'.parameters')) 
		      print ('-------------------', >> (imname//'.parameters'))      
      }


# ------------------------------------------------------------------------------
# Search for up to 'maxnpsf' PSF stars
# ------------------------------------------------------------------------------
#
#     Selection of stars with no INDEF magnitude value
#

      better_psf=no
      first_psf = no
      while (first_psf == no) {
    
	    oldpsfrad = psfrad
	    oldfitrad = fitrad
	    olddmax = dmax    
    
	if (better_psf==no) {
        print ('\n Perform \'pstselect\' (else use existing '//imname//'.pstselect file)? (y/n)')
        print (' default = no')
        do_pstselect = no
        scan (do_pstselect)
	}
	else {
        do_pstselect=yes
	}

        if (do_pstselect==yes) { # Do nothing (move over to 'pstselect' task)
        }
        else {
	          # ------------------------------------------------------------------------------------
	          # Search for 'imname.pstselect' file
	          # ------------------------------------------------------------------------------------
	          auxi=no
	          while (auxi==no) { 

	              o=0

	              # Search for "imname.pstselect"
	              # ------------------------------------------------------------------------------------
	              files ('*'//imname//'.pstselect', > 'tempfile')
	              flist = 'tempfile'
	              while (fscan (flist,line) != EOF) {
	              o = o + 1
	              }
	              del ("tempfile")
	              # ------------------------------------------------------------------------------------

	              if (o!=1) {
		                print ('\n ***********************************')
		                print ('   Missing '//imname//'.pstselect file.')
		                print (' ***********************************')
		                print ('\n Search for '//imname//'.pstselect file again (else run')
		                print (' \'pstselect\' task to create one)? (y/n)')
		                print (' default = yes')
		                check = yes
		                scan (check)
		                if (check) {
		                    auxi=no
		                }
		                else {
		                    auxi=yes
		                    do_pstselect = yes # Move on to 'pstselect' task
		                }
	              }
	              else {
		                if (o==1) {
	                      print ('\n '//imname//'.pstselect file found.')
      delete (imname//'.pstselectold',verify=no,>>&"/dev/null")
      copy ((imname) // '.pstselect', (imname) // '.pstselectold')
	                      auxi=yes
	                      do_pstselect = no # Do not perform 'pstselect' task
		                }
		                else {
	                      print (' Unknown error. Check code.')
	                      bye()
		                }
	              }
	          }
	          # ------------------------------------------------------------------------------------
	          # End of 'Search for 'imname.pstselect' file'
	          # ------------------------------------------------------------------------------------
        }
    
		    if (do_pstselect) {

            # Search for "imname.pstselect"
            # ------------------------------------------------------------------------------------
            o=0
            files ('*'//imname//'.pstselect', > 'tempfile')
            flist = 'tempfile'
            while (fscan (flist,line) != EOF) {
            o = o + 1
            }
            del ("tempfile")
            # ------------------------------------------------------------------------------------
            if (o!=1) { # Do nothing
            }
            else {
                if (o==1) {
                    print ('\n '//imname//'.pstselect file found.')
                    print (' Delete it? (y/n)')
                    print (' default = yes')
                    check=yes
                    scan (check)
                    if (check) {
                        delete (imname//'.pstselect')
                    }
                    else {
                        print ('\n Then move the file to another folder, because')
                        print (' if this file exists it WILL be deleted.')
                        print (' Continue (make sure the the file was moved) Press y/n key')
                        scan (check)
                        delete (imname//'.pstselect',verify=no,>>&"/dev/null")
                    }
                }
            }            

	    print ('\n ---------------------------------------------------------------')
	    print (' Pstselect task (automatically select candidates for PSFs stars)')   
	    print ('---------------------------------------------------------------')
	    print ('\n---------------', >> (imname//'.parameters'))    
	    print ('Pstselect task', >> (imname//'.parameters'))    
	    print ('\n Current \'psfrad\' , \'fitrad\' and \'datamax\' values:')
	    print (' Psfrad = '//psfrad)
	    print (' Fitrad = '//fitrad)
	    print (' Datamax = '//dmax)
	    print ('\n Change one or more of this values before running \'pstselect\' task? (y/n)')
	    print (' default = yes')
	    check = yes
	    scan (check)
	    if (check) {
	        print ('\n Input new \'psfrad\'')
	        scan (psfrad)
	        print ('\nInputed Psfrad = ', psfrad, >> (imname//'.parameters'))
	        daopars.psfrad = psfrad
	        print ('\n Input new \'fitrad\'')
	        scan (fitrad)
	        print ('Inputed Fitrad = ', fitrad, >> (imname//'.parameters'))
	        daopars.fitrad = fitrad
	        print ('\n Input new \'dmax\'')
	        scan (dmax)
	        print ('Inputed Datamax = ', dmax, >> (imname//'.parameters'))
	        datapars.datamax = dmax      
	    }
	    else {
	        print ('\nPsfrad = ', psfrad, >> (imname//'.parameters'))
	        print ('Fitrad = ', fitrad, >> (imname//'.parameters'))
	        print ('Datamax = ', dmax, >> (imname//'.parameters'))
	    }

	    pstselect = yes
	    while (pstselect) {
	    
	      delete ((imname//'.pstselect'),verify=no,>>&"/dev/null")    
	      print ('\n ----------------------------------------------------- ')
	      print (' Pstselect task. Number of stars to find is set to 2% of ')
	      print (' the number of stars found with MAG!=INDEF; but not less')
	      print (' than 20. For this frame, maxnpsf = '//maxnpsf)
	      print (' Image: '// imname // '.fit')
	      print ('\n ----------------------------------------------------- ')      
	      
	      print ('\n Let \'pstselect\' select '//maxnpsf//' stars automatically and then')
	      print (' go through this stars to see which one to keep? (y/n)')
	      print (' Else: select the PSF stars interactively (ie: manually, one by one)')
	      print ('\n Note: you can change the number of PSF stars to be ')
	      print (' automatically selected on the next step')
	      print ('\n default = yes')
	      pstauto=yes
	      scan (pstauto)
	      if (pstauto==no) {
	          pstselect.interactive = yes
	          print ('\n         Keystroke Commands')
	          print ('\n p       Print photometry for star nearest the cursor')
	          print (' l       List the current psf stars')
	          print (' n       Select the next good candidate psf star from the list')
	          print (' a       Add star nearest cursor to psf star list')
	          print (' d       Delete psf star nearest cursor from psf star list')
	          print ('\n         Colon Commands')
	          print ('\n :p [n]  Print photometry for star n'
	          print (' :a [n]  Add star n to psf star list')
	          print (' :d [n]  Delete star n from psf star list')
	          print ('\n When the selection process is over:')
      print (' On the image:                                         ')
      print (' - key "q" to quit from the image                      ')
      print (' On the terminal:                                      ')
      print (' - key "q" to follow with the photometry               ')
      print (' ----------------------------------------------------- ')
	      }
	      else {
	          pstselect.interactive = no      
	      }
	      
	      pstselect.mode = "ql"      
	      pstselect.maxnpsf = maxnpsf
	      pstselect.mkstars = yes
	      pstselect.datapars = ""
	      pstselect.daopars = ""
	      pstselect.plottype = "mesh"
	      pstselect.icommands = ""
	      pstselect.gcommands = ""
	      pstselect.verify = no
	      pstselect.verbose = yes
	      pstselect ((imname), (imname//'.mag.1'), (imname//'.pstselect')) # Output file is 'imname.pstselect'
	      
	      print ('\n ----------------------------------------------------- ')
	      print (' Display task                                          ')
	      print (' Image: '// imname // '.fit')
	      print (' ----------------------------------------------------- ')
	      display ((imname), 1)
	  
	      print ('\n ----------------------------------------------------- ')
	      print (' Tvmark                                                ')
	      print ('   - red point = selected PSF star                   ')
	      print (' ----------------------------------------------------- ')
	      txdump.mode = 'hl' 
	      txdump.textfile = (imname // '.pstselect')
	      txdump.headers = no
	      txdump.fields = 'xcenter, ycenter, id'
	      txdump.expr = 'yes'
	      txdump > auxiliar
	      tvmark.interactive = no
	      tvmark.mark = 'point'
	      tvmark.radii = 10
	      tvmark.font = "raster"
	      tvmark.color = 204
	      tvmark.label = yes
	      tvmark.toleran = 3.5
	      tvmark (1, 'auxiliar')
	      del ('auxiliar')

	      if (pstauto==no) { # Do nothing
	      }
	      else {
      print ('\n You will be presented with the stars \'pstselect\' selected.')
      print (' You can only remove stars from this list. When you are finished,')
      print (' you can decide if you want to change some parameters and repeat')
      print (' the task, or move on keeping the selected PSF stars.')
      print ('\n ----------------------------------------------------- ')
      print (' Go through the PSF stars selected by \'pstselect\':   ')
      print ('\n On the graphical terminal:                          ')
	          print ('\n         Keystroke Commands')
	          print (' ?       Print help')
	          print (' p       Print photometry for current star')
	          print (' a       Accept star and proceed')
	          print (' d       Reject star and proceed')
	          print ('\n When the selection is over:')		      
      print (' On the image:                                         ')
      print (' - key "w" to write the selected PSF stars, and        ')
      print (' - key "q" to quit from the image                      ')
      print (' On the terminal:                                      ')
      print (' - key "q" to follow with the photometry               ')
      print (' ----------------------------------------------------- ')
      psf.interactive = yes
      psf.datapars = ""
      psf.daopars = ""
      psf.matchbyid = yes
      psf.mkstars = yes
      psf.showplots = yes
      psf.plottype = "radial"
      psf.icommands = ""
      psf.gcommands = ""
      psf.verif = no
      
      # Re-order 'imname.pstselect' file by 'YCENTER'
      txsort.mode = 'hl'
		        txsort.ascend = yes              # Sorts stars in ascending YCENTER order
		        txsort ((imname//'.pstselect'), 'YCENTER')		      
      
      psf ((imname), "default", (imname//'.pstselect'), "default", "default", "default")
      delete (imname//'.psg.1')
      delete (imname//'.psf.1.fits')
      delete (imname//'.pstselect')		      		      
      copy ((imname) // '.pst.1', (imname) // '.pstselect') # This file will feed the 'psf' task with
      delete (imname//'.pst.1')                             # the PSF stars
      delete (imname//'.pstselectold',verify=no,>>&"/dev/null")
      copy ((imname) // '.pstselect', (imname) // '.pstselectold')
	      }
                
	      print ('\n Original \'psfrad\' , \'fitrad\' and \'datamax\' values:')      
	      print (' Psfrad = '//oldpsfrad)
	      print (' Fitrad = '//oldfitrad)
	      print (' Datamax = '//olddmax)
	      print ('\n Current \'psfrad\' , \'fitrad\' and \'datamax\' values:')
	      print (' Psfrad = '//psfrad)
	      print (' Fitrad = '//fitrad)
	      print (' Datamax = '//dmax) 
	      print ('\n Keep the selected PSF stars and move on to the PSF photometry? (y/n)')
	      print (' Else: repeat \'pstselect\' with different \'psfrad\', \'fitrad\' and \'datamax\' ')
	      print (' default = yes')
	      check = yes
	      scan (check)
	      if (check) {
	          pstselect = no # This exits the 'pstselect' 'while'
	      }
	      else {
	          delete (imname//'.pstselect')                              
	          print ('\n Input new \'psfrad\'')
	          scan (psfrad)
	          print ('New Psfrad = ', psfrad, >> (imname//'.parameters'))
	          daopars.psfrad = psfrad
	          print ('\n Input new \'fitrad\'')
	          scan (fitrad)
	          print ('New Fitrad = ', fitrad, >> (imname//'.parameters'))
	          daopars.fitrad = fitrad
	          print ('\n Input new \'dmax\'')
	          scan (dmax)
	          print ('New Datamax = ', dmax, >> (imname//'.parameters'))
	          datapars.datamax = dmax
	      }
	      delete (imname//'.pst.1',verify=no,>>&"/dev/null")
	      delete (imname//'.psg.1',verify=no,>>&"/dev/null")
	      delete (imname//'.psf.1.fits',verify=no,>>&"/dev/null")     

	    }  # This bracket closes the 'pstselect' 'while'

	    
	    print ('\n End of Pstselect task')
	    print (' ---------------------')
	    print ('\nEnd of Pstselect task', >> (imname//'.parameters'))
	    print ('--------------', >> (imname//'.parameters'))
	    psfrad = oldpsfrad
	    fitrad = oldfitrad
	    dmax = olddmax
	    datapars.datamax = dmax #
	    daopars.psfrad = psfrad # We set this parameters to their original values
	    daopars.fitrad = fitrad #
		    
		    } # This bracket closes the 'do_pstselect' 'if' 
		    else {
		        print ('\n Existing \'.pstselect\' file used', >> (imname//'.parameters'))
		    }
    
# ------------------------------------------------------------------------------
# End of 'Search for up to 'maxnpsf' PSF stars'
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
#  Core PSF photometry
# ------------------------------------------------------------------------------

       better_psf=no
       better_psf2=no
       core_psf=yes
       while (core_psf==yes) {

          print ('', >> (imname//'.parameters'))
          print ('------------------------------------------------', >> (imname//'.parameters'))
          print ('PSF and neighbors stars substraction', >> (imname//'.parameters')) 
          print ('', >> (imname//'.parameters'))

          if (better_psf) { # Do nothing
          }
          else {
		          print (' ----------------------------------------------------- ')
		          print (' Display task                                          ')
		          print (' Image: '// imname // '.fit')
		          print (' ----------------------------------------------------- ')
		          display ((imname), 1)
		      
		          print ('\n ----------------------------------------------------- ')
		          print (' Tvmark                                                ')
		          print ('   - blue circle = selected PSF star                   ')
		          print (' ----------------------------------------------------- ')
		          txdump.mode = 'hl' 
		          txdump.textfile = (imname // '.pstselect')
		          txdump.headers = no
		          txdump.fields = 'xcenter, ycenter, id'
		          txdump.expr = 'yes'
		          txdump > auxiliar
		          tvmark.interactive = no
		          tvmark.mark = 'circle'
		          tvmark.radii = 10
		          tvmark.font = "raster"
		          tvmark.color = 206
		          tvmark.label = yes
		          tvmark.toleran = 3.5
		          tvmark (1, 'auxiliar')
		          del ('auxiliar')
          }

          print ('\n Fitrad = '//fitrad)
          print (' PSF radius = '//oldpsfrad//' (4xFWHM+1)')
          if (inter) {
		          print ('\n Accept "psfrad" (else input new one)? (y/n)')
		          print (' default = yes')
		          check = yes
		          scan (check)
		          if (check) {
		              psfrad = oldpsfrad
		              daopars.psfrad = psfrad
		          }
		          else {
		              print ('\n Input new "psfrad" value')
		              scan (psfrad)
		              daopars.psfrad = psfrad
		              print ('New Daopars.psfrad (First PSF run) = ', psfrad, >> (imname//'.parameters'))
		          }
          }
          else {
              psfrad = oldpsfrad
              daopars.psfrad = psfrad          
          }

          print ('\n Set constant PSF model (varorder = 0): \'y\'')
          print (' Set variable PSF model (varorder = 2): \'n\'')
          check = no
          scan (check)
          if (check) {
              daopars.varorder = 0
              print ('Constant PSF model (varorder = 0)', >> (imname//'.parameters'))
          }
          else {
              daopars.varorder = 2
              print ('Variable PSF model (varorder = 2)', >> (imname//'.parameters'))
          }
   
          print ('\n -----------------------------------------------------')
          print (' Psf task                                              ')
          print (' PSF computation (First run)                           ')
          print (' Psfrad = '//psfrad)
          print (' Fitrad = '//fitrad)
          psf.interactive = no
          psf.datapars = ""
          psf.daopars = ""
          psf.matchbyid = yes
          psf.mkstars = yes
          psf.showplots = yes
          psf.plottype = "radial"
          psf.icommands = ""
          psf.gcommands = ""
          psf.verif = no
          psf ((imname), "default", (imname//'.pstselect'), "default", "default", "default")
          
		      files ('*'//imname//'.psg.1', > 'temp.grd')
		      file_var = 'temp.grd'
		      n=0
		      while (fscan (file_var,line) != EOF) {
		          n = n + 1
		      }
		      if (n == 0) {
              nopsfstar = yes # This means 0 stars were selected as PSF and so, the PSF selection must
		      }                   # be performed again
		      else {
              nopsfstar = no   
              copy ((imname) // '.psg.1', (imname) // '.grp.1')		              
		      }
		      del ('temp.grd')               
          
          # At this point a new .pst file will have been created: "imname.pst.1"
          # with the new PSF stars.          

          if (nopsfstar == no) {

		          first_nstarsubstar = no
		          while (first_nstarsubstar== no) {
		              if (dopsf==yes) {
	              print (' ----------------------------------------------------- ')
	              print (' Nstar task (first run)                                ')
	              print (' PSF photometry of PSF and their neighbour stars       ')
	              print (' ----------------------------------------------------- ')
		              }
		              else {
	              print (' ----------------------------------------------------- ')
	              print (' Nstar task (first run)                                ')
	              print (' Aperture photometry of PSF and their neighbour stars  ')
	              print (' ----------------------------------------------------- ')		              
		              }
		              if (better_psf && better_psf2) {
	              print ('\n Fitrad = '//fitrad)
	              print ('\n Original \'psfrad\' = '//psfrad)
          print (' Actual \'psfrad\' = '//nstarpsfrad)
          print (' Input new value:')
          scan (nstarpsfrad)
		                  daopars.psfrad = nstarpsfrad
		                  print ('New Daopars.psfrad (First Nstar run) = ', nstarpsfrad, >> (imname//'.parameters'))
		              }
		              else {
	              daopars.psfrad = psfrad - fitrad
	              print ('Daopars.psfrad = ', (psfrad - fitrad), >> (imname//'.parameters'))
	              print ('\n Fitrad = '//fitrad)
	              print (' PSF radius = '//psfrad - fitrad//' (previous psfrad - fitrad)')
	              if (inter) {
              print ('\n Accept "psfrad"; should be \'previous psfrad\' - \'fitrad\' ')
              print (' according to Massey-Davis guide (else input new one)? (y/n)')
              print (' default = yes')
              check = yes
              scan (check)
              if (check) {
                  nstarpsfrad = psfrad
              }
              else {
                  print ('\n Input new "psfrad" value')
                  scan (nstarpsfrad)
              }
	              }
	              else {
              nstarpsfrad = psfrad
	              }
		                  daopars.psfrad = nstarpsfrad
		                  print ('New Daopars.psfrad (First Nstar run) = ', nstarpsfrad, >> (imname//'.parameters'))
		              }

		              nstar.verif = no
		              nstar ((imname), "default", "default", "default", "default")
		        
		              print ('\n ----------------------------------------------------- ')
		              print (' Substar task (first run)                              ')
		              print (' Sustraction of PSFs and their neighbour stars         ')
		              print (' ----------------------------------------------------- ')
		              print (' PSF radius = '//nstarpsfrad)     
		              print (' Fitrad = '//fitrad)
		              substar.verif = no
		              #      
		              # If the third file in 'substar' is "default", it will exclude from the
		              # substraction the PSF stars and thus only substract their neighbours. If
		              # this file is "", it will substract all the stars in the '.nst' file
		              # (ie: PSFs and their neighbours)
		              #      
		              substar ((imname), "default", "", "default", "default")
		          
                  if (better_psf) {
	              print ('\n ----------------------------------------------------- ')
	              print (' Display task                                          ')
	              print (' Image without PSFs or their neighbours: '//imname//'.sub.1.fits ')
	              print ('\n IMPROVED VERSION (maybe) DISPLAYED IN FRAME 2')
	              print (' ----------------------------------------------------- ')
                      display ((imname // '.sub.1.fits'), 2)
                  }
                  else {
	              print ('\n ----------------------------------------------------- ')
	              print (' Display task                                          ')
	              print (' Image without PSFs or their neighbours: '//imname//'.sub.1.fits ')
	              print (' ----------------------------------------------------- ')
                      display ((imname // '.sub.1.fits'), 1)
                  }

		              print ('\n ----------------------------------------------------- ')
		              print (' Tvmark                                                ')
		              print ('   - red point = detected star                         ')
		              print ('   - blue circle = selected PSF star                   ')
		              print (' ----------------------------------------------------- ')
		              tvmark.interactive = no       
		              tvmark.mark = 'point'
		              tvmark.font = "raster"
		              tvmark.color = 204
		              tvmark.label = no
		              tvmark.number = no
		              tvmark.toleran = 3.5
		              if (better_psf) {
		                  tvmark (2, imname // '.coo.1')
		              }
		              else {
		                  tvmark (1, imname // '.coo.1')
		              }
		              txdump.textfile = (imname // '.pst.1')
		              txdump.headers = no
		              txdump.fields = 'xcenter, ycenter, id'
		              txdump.expr = 'yes'
		              txdump > auxiliar
		              tvmark.interactive = no      
		              tvmark.mark = 'circle'
		              tvmark.radii = 10
		              tvmark.font = "raster"
		              tvmark.color = 206
		              tvmark.label = yes
		              tvmark.toleran = 3.5
		              if (better_psf) {
		                  tvmark (2, 'auxiliar')
		              }
		              else {
		                  tvmark (1, 'auxiliar')
		              }
		              del ('auxiliar')

                  print ('\n ------------------------------------------------------------')
                  print ('\n Where the PSFs and their neighbors substracted cleanly? (y/n)')
                  print (' (else select among various options)')
                  print (' default = yes')
                  check=yes
                  scan (check)
                  
                  if (check) {
                      better_psf=no
                      first_nstarsubstar = yes  # Leaves the 'first_nstarsubstar' 'while'
                      core_psf = no             # Leaves the 'core_psf' 'while'
                      first_psf = yes           # Leaves the 'first_psf' 'while'
                  }
                  else {
                      auxi=no
                      while (auxi==no) {
	                  print ('\n Choose one of this options:')
	                  print ('\n [1]- Run \'nstar\' and \'substar\' again with different \'psfrad\'')
	                  print ('\n [2]- Remove one or several PSF stars from list and')
	                  print ('      repeat the substraction process')
	                  print ('\n [3]- Select PSF stars again (MANUALLY, ie: \'.pstselect\' file will be DELETED)')
              if (better_psf) {
                  print ('\n [4]- Return to first PSF stars list')
              }
              else { # Do nothing
              }
	                  print ('\n Input option number:')
	                  print (' default = [1]')
	                  var2=1
              scan (var2)
		                      if (better_psf) {
		                          if (var2 != 1 && var2 != 2 && var2 != 3 && var2 != 4) {
		                              auxi=no
		                              print ('\n Wrong input. Try again.')
		                          }
		                          else {
		                              auxi=yes
		                          }
		                      }
		                      else {
		                          if (var2 != 1 && var2 != 2 && var2 != 3) {
		                              auxi=no
		                              print ('\n Wrong input. Try again.')
		                          }
		                          else {
		                              auxi=yes
		                          }
		                      }
                      }
                      better_psf=yes
		              
	              if (var2==1) {
                          better_psf2=yes
                          first_nstarsubstar = no # Stays inside the 'first_nstarsubstar' 'while'
                          del (imname // '.nst.1')
                          del (imname // '.nrj.1')
                          del (imname // '.sub.1.fits') 
	              }
	              else {
	                  if (var2==2) {
	                      better_psf2=no
                             first_nstarsubstar = yes # Leaves the 'first_nstarsubstar' 'while', but stays
                                                       # inside the 'core_psf' 'while'                              
                              rem_psf=yes
                              while (rem_psf==yes) {
		                              print ('\n Input the number of the PSF star you wish to remove:')
		                              print (' (you can remove more after this)')
		                              scan (psfnumber)
		                              # Remotion process ----------------------
                                  expression = "ID != "//psfnumber

		      print ('\n Expression: '//expression)

		      rename.field = 'all'
		      rename (imname//'.pstselect', imname//'_pstoriginal')

		      pselect.infiles = imname//'_pstoriginal'
		      pselect.outfiles = imname//'.pstselect'
		      pselect.expr = expression
		      pselect.mode = "hl"
		      pselect
		      delete (imname//'_pstoriginal')
		                              # Remotion process ----------------------		                              
		                              print ('\n Remove another one (else continue)? (y/n)')
		                              print (' default = yes')
		                              check=yes
		                              scan (check)
		                              if (check) {
		                              }
		                              else {
                                      rem_psf=no
                          del (imname // '.psf.1.fits')
                          del (imname // '.grp.1')
                          del (imname // '.psg.1')
                          del (imname // '.pst.1')
                          del (imname // '.nst.1')
                          del (imname // '.nrj.1')
                          del (imname // '.sub.1.fits')                                      
		                              }
                              }
	                  }
	                  else {
	                      if (var2==3) {
	                          better_psf2=no
		                              print ('\n Selecting PSF stars again')
		                              print ('Selecting PSF stars again', >> (imname//'.parameters'))
	                          first_nstarsubstar = yes       # Leaves the 'first_nstarsubstar' 'while',
	                          core_psf = no                  # and leaves the 'core_psf' 'while' but stays
	                          del (imname // '.psf.1.fits')  # inside the 'first_psf' 'while'
	                          del (imname // '.grp.1')
	                          del (imname // '.psg.1')
	                          del (imname // '.pst.1')
	                          del (imname // '.nst.1')
	                          del (imname // '.nrj.1')
	                          del (imname // '.sub.1.fits')
	                          del (imname // '.pstselect')
	                      }
	                      else {
	                          if (var2==4) {
	                              print ('\n Selecting first PSF stars list again')
	                              print ('Selecting first PSF stars list again', >> (imname//'.parameters'))
                          first_nstarsubstar = yes       # Leaves the 'first_nstarsubstar' 'while'
                          del (imname // '.psf.1.fits')
                          del (imname // '.grp.1')
                          del (imname // '.psg.1')
                          del (imname // '.pst.1')
                          del (imname // '.nst.1')
                          del (imname // '.nrj.1')
                          del (imname // '.sub.1.fits')
                          del (imname // '.pstselect')
                          copy ((imname) // '.pstselectold', (imname) // '.pstselect')
                          better_psf=no
	                              better_psf2=no
	                          }
	                      }
	                  }
	              }
	          }    


		          } # This bracket closes the 'first_nstarsubstar' while 

		          print ('', >> (imname//'.parameters'))
		          print ('End of PSF and neighbors stars substraction', >> (imname//'.parameters')) 
		          print ('------------------------------------------------', >> (imname//'.parameters'))
		          print ('', >> (imname//'.parameters'))
          
          } # This bracket closes the 'nopsfstar' 'if'
          else {
              first_nstarsubstar = yes       # Leaves the 'first_nstarsubstar' 'while', but stays
                                             # inside the 'first_psf' 'while'
          }
          
        } # This bracket closes the 'core_psf' while 

      } # This bracket closes the 'first_psf' while 
      
      delete (imname//'.pstselectold',verify=no,>>&"/dev/null")
# ------------------------------------------------------------------------------           

      file_var = (imname//'.pst.1')
      n = 0
      m = 0
      while (fscan (file_var,line) != EOF) {
          n = n + 1
      }
      file_var = (imname//'.nst.1')	      
      while (fscan (file_var,line) != EOF) {
          m = m + 1
      }	      

      if ((m-42)/2 == (n-24)) {  # If this is true, that means there are NO PSF neighbour stars.
          no_neighbours = yes
          print ('\n No PSF neighbour stars. Skipping second \'substar\' run.')
          print ('\n (\'Substar\' will still run just so the numbering')
          print ('  on the files is right. This does\'nt affect the outcome)')
          print (' No PSF neighbour stars. Skipping second \'substar\' run.', >> (imname//'.parameters'))          
          copy ((imname) // '.fit', (imname) // '.sub.2.fits')
      }
      else {
          no_neighbours = no
		      print ('\n ----------------------------------------------------- ')
		      print (' Substar task  (second run)                            ')
		      print (' Sustraction of PSF neighbours stars                   ')
		      print (' ----------------------------------------------------- ')
		      print (' PSF radius = '//nstarpsfrad)     
		      print (' Fitrad = '//fitrad)
		      #      
		      # If the third file in 'substar' is "default", it will exclude from the
		      # substraction the PSF stars and thus only substract their neighbours. If
		      # this file is "", it will substract all the stars in the '.nst' file
		      # (ie: PSFs and their neighbours)
		      #         
		      substar.verif = no
		      substar ((imname), "default", "default", "default", "default")

		      #
		      # At this point a file "imname.sub.2.fits" should have been created
		      # This file contains all the stars, except the PSFs neighbours
		      #	              
      }

      print ('\n ----------------------------------------------------- ')
      print (' Psf task (second run, without neighbours)             ')
      print (' PSF computation                                       ')
      print (' ----------------------------------------------------- ')
#-----------------------------------------------------------------------------------------------      
      daopars.psfrad = psfrad # This goes back to the 'psfrad' value used when the first PSF was performed
      print ('\n Fitrad = '//fitrad)
      print (' PSF radius = '//psfrad)     
      if (no_neighbours==no) {
          if (inter) {
	      print ('\n Accept PSF radius "psfrad"; should be the same one you used when')
	      print (' the PSF selection was made, according to Massey (else input new one)? (y/n)')
	      print (' default = yes')
	      check = yes
	      scan (check)
	      if (check) { # Do nothing
	      }
	      else {
	          print ('\n Input new "psfrad" value')
	          scan (psfrad)
	          daopars.psfrad = psfrad
	          print ('New Daopars.psfrad (Second PSF run) = ', psfrad, >> (imname//'.parameters'))
	      }
		      }
      }
#-----------------------------------------------------------------------------------------------       
     
      if (no_neighbours==no) {
          print ('\n datamin = '//datamin)
          print (' 2*sigma = '//2*sigma)
          if (inter) {
		          print ('\n Diminish \'datamin\' value by 2 SIGMAS? (y/n)')
		          print ('\n (New datamin value will be = '//datamin -2*sigma//' )')
		          print (' default = yes')
		          check = yes
		          scan (check)
		          if (check) {
      datapars.datamin = datamin -2*sigma # Due to added noise by subtraction, we reset 'datamin'
                                                # to 2 sigma below the previous value.
      print ('New Datapars.datamin (Second PSF run) = ', (datamin -2*sigma), >> (imname//'.parameters'))
	}
	else { # Leave 'datamin' as it is
	}
	  }
	  else {
	      datapars.datamin = datamin -2*sigma # Due to added noise by subtraction, we reset 'datamin'
	                                                # to 2 sigma below the previous value.
	      print ('New Datapars.datamin (Second PSF run) = ', (datamin -2*sigma), >> (imname//'.parameters'))	  
	  }
	}
	else {
        datapars.datamin = datamin -2*sigma
        print ('New Datapars.datamin (Second PSF run) = ', (datamin -2*sigma), >> (imname//'.parameters'))
	}
      psf.interac = no
      psf ((imname // '.sub.2.fits'), (imname // '.mag.1'), (imname) // '.pst.1', "default", "default", "default")

      txdump.textfile = (imname // '.pst.1')
      txdump.headers = no
      txdump.fields = 'xcenter, ycenter'
      txdump.expr = 'yes'
      txdump > auxiliar
      copy ('auxiliar', (imname) // '.psf1')
      #
      # This file (imname.psf1) which only contains the coordinates of the PSF selected stars,
      # will be used for the Aperture Correction further bellow.
      #
      del ('auxiliar')
      
#----------------------------------------------------------------------------------------------
# Check how good is the revised PSF

      copy ((imname // '.sub.2.fits.psg.1'), (imname // '.psg1'))           # Files to be used by the
      copy ((imname // '.sub.2.fits.psf.1.fits'), (imname // '.psf2.fits')) # next 'nstar' and 'substar' run (below)

      print ('\n ----------------------------------------------------- ')
      print (' Nstar task (second run)                               ')
      print (' PSF photometry of PSF and their neighbour stars       ')
      print (' ----------------------------------------------------- ')
      daopars.psfrad = psfrad
      print ('Daopars.psfrad = ', psfrad, >> (imname//'.parameters'))
      print ('\n Fitrad = '//fitrad)
      print (' PSF radius = '//psfrad)
      if (no_neighbours==no) {
          if (inter) {
		          print ('\n Accept PSF radius "psfrad"; should be the same one you used when')
	      print (' the PSF selection was made, according to Massey (else input new one)? (y/n)')
	      print (' default = yes')
	      check = yes
		          scan (check)
		          if (check) { # Do nothing
		          }
		          else {
		              print ('\n Input new "psfrad" value')
		              scan (psfrad)
		              daopars.psfrad = psfrad
		              print ('New Daopars.psfrad (Second Nstar run) = ', psfrad, >> (imname//'.parameters'))
		          }
          }
      }

      nstar.verif = no
      nstar ((imname), (imname // '.psg1'), (imname // '.psf2.fits'), "default", "default")      

      print ('\n ----------------------------------------------------- ')
      print (' Substar task (third run)                              ')
      print (' Sustraction of PSFs and their neighbour stars         ')
      print (' ----------------------------------------------------- ')
      print (' PSF radius = '//psfrad)     
      print (' Fitrad = '//fitrad)
      substar.verif = no
      #      
      # If the third file in 'substar' is "default", it will exclude from the
      # substraction the PSF stars and thus only substract their neighbours. If
      # this file is "", it will substract all the stars in the '.nst' file
      # (ie: PSFs and their neighbours)
      #      
      substar ((imname), "default", "", (imname // '.psf2.fits'), "default")

      if (no_neighbours==no) {
          if (inter) {
		          print ('\n ----------------------------------------------------- ')
		          print (' Display task                                          ')
		          print (' Frame 1: substraction with PREVIOUS PSF model         ')
		          print (' Frame 2: substraction with CURRENT PSF model          ')              
		          print (' ----------------------------------------------------- ')
		          display ((imname // '.sub.1.fits'), 1)
		          display ((imname // '.sub.3.fits'), 2)
		          
		          print ('\n Can you see an improvement in this new version of the substracted')
		          print (' frame (using the CURRENT PSF model) over the previous one (using')
		          print (' the PREVIOUS PSF model)?')
		          print ('\n Are any neighbors still visible?')
		          print ('\n Do you want to keep one of this two PSF models and move on (else')
		          print (' perform the \'neighbors substraction process\' and \'psf\' task one')
		          print (' more time and THEN move on)? (y/n)')
		          print (' default = yes')
		          check = yes
		          scan (check)
          }
          else {
              check = yes
          }
      }
      else {
          check = yes
      }
      if (check) {
          if (no_neighbours==no) {
              if (inter) {
		              print ('\n Keep CURRENT PSF model (else keep PREVIOUS one)? (y/n)')
		              print (' default = yes')
		              auxi=yes
		              scan (auxi)
              }
              else {
                  auxi = yes
              }
          }
          else {
              auxi = yes
          }
          if (auxi) {
              delete ((imname // '.sub.3.fits')) # This is just to keep the numbering right
              delete ((imname // '.nst.2'))
              delete ((imname // '.nrj.2'))
              delete ((imname // '.psg1'))
              delete ((imname // '.psf2.fits'))
          }
          else {
              delete ((imname // '.sub.3.fits'))
              delete ((imname // '.nst.2'))
              delete ((imname // '.nrj.2'))
              delete ((imname // '.psg1'))
              delete ((imname // '.psf2.fits'))
		    delete ((imname // '.sub.2.fits.psf.1.fits'))      
	      copy ((imname // '.psf.1.fits'), (imname // '.sub.2.fits.psf.1.fits'))
	      delete ((imname // '.psf.1.fits'))		      
          }
      }
      else {
          print ('\n ----------------------------------------------------- ')
          print (' Substar task (fourth run)                             ')
          print (' Sustraction of PSFs neighbour stars                   ')
          print (' ----------------------------------------------------- ')
          print (' PSF radius = '//psfrad)     
          print (' Fitrad = '//fitrad)
          substar.verif = no
          #      
          # If the third file in 'substar' is "default", it will exclude from the
          # substraction the PSF stars and thus only substract their neighbours. If
          # this file is "", it will substract all the stars in the '.nst' file
          # (ie: PSFs and their neighbours)
          #      
          substar ((imname), "default", (imname // '.pst.1'), (imname // '.psf2.fits'), "default")   
          
		      print ('\n ----------------------------------------------------- ')
		      print (' Psf task (third run, without neighbours)              ')
		      print (' PSF computation                                       ')
		      print (' ----------------------------------------------------- ')
		      print (' PSF radius = '//psfrad)     
		      print (' Fitrad = '//fitrad)

		      datapars.datamin = datamin -2*sigma # Due to added noise by subtraction, we reset 'datamin'
		                                                # to 2 sigma below the previous value.
		      print ('New Datapars.datamin (Third PSF run) = ', (datamin -2*sigma), >> (imname//'.parameters'))
		      psf.interac = no
		      psf ((imname // '.sub.4.fits'), (imname // '.mag.1'), (imname // '.pst.1'), "default", "default", "default")
		      
		      delete ((imname // '.psg1'))
		      delete ((imname // '.psf2.fits'))
          delete ((imname // '.nst.2'))
          delete ((imname // '.nrj.2'))	
          
		      delete ((imname // '.sub.1.fits'))
		      copy ((imname // '.sub.3.fits'), (imname // '.sub.1.fits'))      
		      delete ((imname // '.sub.3.fits'))
		      
		      delete ((imname // '.sub.2.fits'))
		      copy ((imname // '.sub.4.fits'), (imname // '.sub.2.fits'))      
		      delete ((imname // '.sub.4.fits'))          
		      
	delete ((imname // '.psf.1.fits'))      
		      copy ((imname // '.sub.2.fits.psf.1.fits'), (imname // '.psf.1.fits'))
		      delete ((imname // '.sub.2.fits.psf.1.fits'))		      
		      
		      delete ((imname // '.sub.2.fits.psg.1'))
		      copy ((imname // '.sub.4.fits.psg.1'), (imname // '.sub.2.fits.psg.1'))
		      delete ((imname // '.sub.4.fits.psg.1'))
		      
		      delete ((imname // '.sub.2.fits.psf.1.fits'))
		      copy ((imname // '.sub.4.fits.psf.1.fits'), (imname // '.sub.2.fits.psf.1.fits'))
		      delete ((imname // '.sub.4.fits.psf.1.fits'))
		      
		      delete ((imname // '.sub.2.fits.pst.1'))
		      copy ((imname // '.sub.4.fits.pst.1'), (imname // '.sub.2.fits.pst.1'))
		      delete ((imname // '.sub.4.fits.pst.1'))            
      }

#----------------------------------------------------------------------------------------------         

      if (dopsf==no) {
          right_sustraction = yes # This is so that after the aperture photometry is done
      }                           # the script won't start from the beginning

      if (dopsf==yes) { # This means the PSF photometry will be made (ie: 'allstar' task will be used).

		      sannulu = 2*fitrad
		      daopars.sannulu = sannulu
		      print ('\n --------------------------------')
		      print ('\n Next step: \'Allstar\' task')
		      print ('\n Daopars.sannulu = '//sannulu)
		      if (inter) {
	      print ('\n Accept daopars.sannulu value (else input new one)? (y/n)')
	      print (' default = yes')
	      check = yes
	      scan (check)
	      if (check) { # Do nothing
	      }
	      else {
	          print (' Input new daopars.sannulu value')
	          scan (sannulu)
	          daopars.sannulu = sannulu
	          print ('New Daopars.sannulu (Second PSF run) = ', psfrad, >> (imname//'.parameters'))
	      }
		      }
		      
		      print ('\n ----------------------------------------------------- ')
		      print (' Allstar task                                          ')
		      print (' PSF photometry                                        ')
		      print (' ----------------------------------------------------- ')
		      allstar.verif = no
		      allstar ((imname), "default", (imname // '.sub.2.fits.psf.1.fits'), "default", "default", "default")
		      l = 1 # This number is used to know whether only one or two 'allstar' runs were performed

		      print ('')
		      print ('')
		      print ('\n -----------------------------------------------------------')
		      print (' Display task                                               ')
		      print (' Compare the two frames to see if the stars are cleanly gone')
		      print (' -----------------------------------------------------------')
		      display ((imname // '.fit'), 1)
		      display ((imname // '.sub.3.fits'), 2)

		      print ('\n Is the substraction right (else star over from PSF selection)? (y/n) ')
		      print (' default = yes')
		      right_sustraction=yes
		      scan (right_sustraction)
    
		      if (right_sustraction == yes) { # Do nothing
		      }
		      else {
		          del (imname // '.psf.1.fits')
		          del (imname // '.grp.1')
		          del (imname // '.psg.1')
		          del (imname // '.pst.1')
		          del (imname // '.nst.1')
		          del (imname // '.nrj.1')
		          del (imname // '.sub.1.fits')
		          del (imname // '.sub.2.fits')
		          del (imname // '.sub.2.fits.psf.1.fits')
		          del (imname // '.sub.2.fits.psg.1')
		          del (imname // '.sub.2.fits.pst.1')
		          del (imname // '.psf1')                 
		          del ((imname) // '.als.1')
		          del ((imname) // '.arj.1')        
		          del ((imname) // '.sub.3.fits')
#		          del ((imname) // '.pstselect')
		          delete (imname//'.sub.4.fits.pst.1',verify=no,>>&"/dev/null")
		          delete (imname//'.sub.4.fits.psf.1.fits',verify=no,>>&"/dev/null")
		          delete (imname//'.sub.4.fits.psg.1',verify=no,>>&"/dev/null")                    
		      }

		      print ('', >> (imname//'.parameters'))
		      print ('End of Core PSF Reduction', >> (imname//'.parameters')) 
		      print ('************************************************', >> (imname//'.parameters'))
		      print ('', >> (imname//'.parameters'))
      }

  } # This last bracket closes the 'right_sustraction' while

# ------------------------------------------------------------------------------
# End of 'Core PSF Reduction'
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Second search of stars, previously hidden
# ------------------------------------------------------------------------------

      if (dopsf==yes) {
      
          if (inter) {
	      print ('\n ----------------------------------------------------- ')
	      print (' Second search (add previously hidden stars if any)    ')
	      print (' ----------------------------------------------------- ')
	      print ('\n Stars previously found are shown over the substracted frame')

	      display ((imname // '.sub.3.fits'), 1)
	      tvmark.interactive = no
	      tvmark.mark = 'point'
	      tvmark.font = "raster"
	      tvmark.color = 204
	      tvmark.number = no
	      tvmark.label = no
	      tvmark.toleran = 3.5
	      tvmark (1, (imname // '.coo.1'))   

	      o=0 # This 'int' indicates that the file 'imname.sub.3.fits.coo.1' DOES NOT exists
	      print ('\n Perform search for new (previously hidden) stars? (y/n)')
	      print (' default = no')
	      newstars=no
	      scan (newstars)
		      }
		      else {
		          newstars = no
		      }

		      if (newstars) {

		          o=1 # This variable indicates that the file 'imname.sub.3.fits.coo.1' DOES exists
		              # (it doesn't exist YET, but it'll be created in the next step)
		          star_search = no
		          thresh = 3.5*sigma + 2*sigma
		          print ('Second star search, threshold = ', thresh, >> (imname//'.parameters'))
		          while (star_search == no) {

		                print ('\n Perform daofind search (else mark the stars manually)? (y/n)')
		                print (' default = yes')
		                check = yes
		                scan (check)
		                if (check) {
		                    print ('\n Actual threshold = ' // thresh)
		                    print ('\n Input new threshold value? (y/n)')
		                    print (' default = yes')
		                    check = yes
		                    scan (check)
		                    if (check) {
		                        print ('\n (A higher threshold means less stars)')
		                        scan (thresh)
		                        findpars.threshold = thresh
		                    }
		                    print ('\n ----------------------------------------------------- ')
		                    print (' Daofind task                                          ')
		                    print (' Second search (add previously hidden stars if any)    ')
		                    print (' ----------------------------------------------------- ')
		                    datapars.datamin = datamin -2*sigma
		                    print ('Second star search, datapars.datamin = ', (datamin -2*sigma), >> (imname//'.parameters'))
		                    daofind ((imname // '.sub.3.fits'), "default")
		                    #
		                    # At this point the file 'imname.sub.3.fits.coo.1' should be created

		                    print ('\n ----------------------------------------------------- ')
		                    print (' Display task                                          ')
		                    print (' Daofind results after the first "allstar" run')
		                    print (' ----------------------------------------------------- ')
		                    display ((imname // '.sub.3.fits'), 2)

		                    print ('\n ----------------------------------------------------- ')
		                    print (' Tvmark                                                ')
		                    print ('   - red point = detected star previously hidden       ')
		                    print (' Add or remove stars if necessary                      ')
		                    print (' ----------------------------------------------------- ')
		                    txdump.textfile = (imname // '.sub.3.fits.coo.1')
		                    txdump.headers = no
		                    txdump.fields = 'xcenter, ycenter, id'
		                    txdump.expr = 'yes'
		                    txdump > auxiliar1 

		                    tvmark.interactive = no
		                    tvmark.mark = 'point'
		                    tvmark.font = "raster"
		                    tvmark.color = 204
		                    tvmark.number = no
		                    tvmark.label = yes
		                    tvmark.toleran = 3.5
		                    tvmark (2, 'auxiliar1')
		          
		                    tvmark.interac = yes 
		                    tvmark.mark = 'point'
		                    tvmark.font = "raster"
		                    tvmark.color = 204
		                    tvmark.number = no      
		                    tvmark.label = yes
		                    tvmark.toleran = 3.5
		                    tvmark (2, imname // '.sub.3.fits.coo.1')
		                    del ('auxiliar1')
		                }
		                else {
		                
		                    print ('\n ----------------------------------------------------- ')
		                    print (' Display task                                          ')
		                    print (' Daofind results after the first "allstar" run')
		                    print (' ----------------------------------------------------- ')
		                    display ((imname // '.sub.3.fits'), 2)

		                    print ('\n ----------------------------------------------------- ')
		                    print (' Add stars manually                                    ')
		                    print (' ----------------------------------------------------- ')
		                    findpars.threshold = 5000. # A high enough threshold so that no star will be marked
		                    daofind ((imname // '.sub.3.fits'), "default")
		                    #
		                    # At this point the file 'imname.sub.3.fits.coo.1' should be created

		                    tvmark.interac = yes 
		                    tvmark.mark = 'point'
		                    tvmark.font = "raster"
		                    tvmark.color = 204
		                    tvmark.number = no      
		                    tvmark.label = yes
		                    tvmark.toleran = 3.5
		                    tvmark (2, imname // '.sub.3.fits.coo.1')
		                    findpars.threshold = thresh
		                }

		                # At this point the file 'imname.sub.3.fits.coo.1' exists,
		                # created either by Daofind or manually
		                #
		                print ('\n Perform new search (ie: delete created \'.coo\' file')
		                print (' and repeat the second search for stars; else continue) ? (y/n)')
		                print (' default = yes')
		                check = yes
		                scan (check)
		                if (check) {
		                    star_search = no
		                    del (imname// '.sub.3.fits.coo.1')
		                }
		                else {
		                    star_search = yes
		                }
		          } # This bracket closes the 'star_search' while

		          # ------------------------------------------------------------------------------
		          # Add new stars, previously hidden (if any) and perform Phot and Allstar again
		          # ------------------------------------------------------------------------------
		                print ('\n Add new (previously hidden) stars found and marked')
		                print (' in this second search (else disregard new found stars)? (y/n)')
		                print (' default = yes')
		                check = yes
		                scan (check)
		                if (check) {
		                    print ('\n ----------------------------------------------------- ')
		                    print (' Phot task (second run)                                ')
		                    print (' ----------------------------------------------------- ')
		                    phot.interactive = no
		                    phot.radplots = no
		                    phot.update = yes
		                    phot.verbose = yes
		                    phot.verify = no
		                    photpars.apertures = aperture
		                    phot ((imname), (imname // '.sub.3.fits.coo.1'), "default")
		                    
		         #          tmerge ((imname//'.mag.1')//","//(imname//'.mag.2'),(imname//'.mag.3'),"append")
		         #          pconcat (infiles=imname//'.mag.1,'//imname//'.mag.2', outfile=imname//'.mag.3')
		                    txconcat ((imname // '.mag.1' // ',' // imname // '.mag.2'), (imname // '.mag.3'))          
		                    # This line replaces the obsolete 'pappend' command
		          
		                    print ('                                                       ')      
		                    print (' ----------------------------------------------------- ')
		                    print (' Allstar task                                          ')
		                    print (' PSF photometry (second run)                           ')
		                    print (' ----------------------------------------------------- ')
		                    allstar.verif = no
		                    allstar ((imname), (imname // '.mag.3'), (imname // '.sub.2.fits.psf.1.fits'), "default", "default", "default")
		                    l = 2 # This number is used to know whether only one or two 'allstar' runs were performed
		                 }
		                 else { # Do nothing
		                 }

		          # ---------------------------------------------------------------------------------------------------
		          # End of 'Add new stars, previously hidden (if any) and perform Phot and Allstar again'
		          # ---------------------------------------------------------------------------------------------------

		      } # This bracket closes the 'newstars' if
		      else { # Do nothing
		      }
      }

# ------------------------------------------------------------------------------
# End of 'Second search of stars, previously hidden'
# -----------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Check .coo.1, .als.1 and .als.2 (if such file exists) files
# -----------------------------------------------------------------------------

      if (dopsf==yes) {
      
		      if (l==1) { # l = 1 means NO 'Second stars search' was performed.
		          if (inter) {
	          print ('\n ----------------------------------------------------- ')
	          print (' Display task                                          ')
	          print (' Display .als.1 (Frame1) and .coo.1 (Frame2) files')
	          print (' ----------------------------------------------------- ')
	          display ((imname // '.fit'), 1)
	          display ((imname // '.fit'), 2)

	          print ('\n ----------------------------------------------------- ')
	          print (' Tvmark                                                ')
	          print ('   - Frame1: detected stars by "allstar"               ')
	          print ('   - Frame2: detected stars by "daofind"               ')
	          print (' ----------------------------------------------------- ')
	          txdump.textfile = (imname // '.als.1')
	          txdump.headers = no
	          txdump.fields = 'xcenter, ycenter'
	          txdump.expr = 'yes'
	          txdump > auxiliar1 

	          tvmark.interactive = no
	          tvmark.mark = 'point'
	          tvmark.font = "raster"
	          tvmark.color = 204
	          tvmark.number = no
	          tvmark.label = no
	          tvmark.toleran = 3.5
	          tvmark (1, 'auxiliar1')
	          
	          tvmark.interac = no 
	          tvmark.mark = 'point'
	          tvmark.font = "raster"
	          tvmark.color = 204
	          tvmark.number = no      
	          tvmark.label = no
	          tvmark.toleran = 3.5
	          tvmark (2, imname // '.coo.1')
	          del ('auxiliar1')
		          }

		          check=no
		          while (check == no) {

                  if (inter) {
          print ('\n What do you want to do now?:')
          print ('     [1]: keep \'.als.1\' file and continue,')
          print ('     [2]: Start over from \'Daofind\' search')
          print ('          (the very beginning!!)')

	              scan (var1)
		              }
		              else {
		                  var1 = 1
		              }
		              if (var1==1) {
		                  check=yes
		                  psf_done=yes
		              }
		              else {
		                  if (var1==2) { # DELETE ALL FILES AND START FROM THE BEGINNING
		                      print ('\n Are you sure you want to delete all files')
		                      print (' and start over from the beginning? (y/n)')
		                      print (' default = no')
		                      check=no
		                      scan (check)
		                      if (check) { 
	                      check=yes
	                      del (imname // '.psf.1.fits')
	                      del (imname // '.grp.1')
	                      del (imname // '.psg.1')
	                      del (imname // '.pst.2')
	                      del (imname // '.nst.1')
	                      del (imname // '.nrj.1')
	                      del (imname // '.psf1')
	                      del (imname // '.sub.1.fits')
	                      del (imname // '.sub.2.fits')
	                      del (imname // '.sub.2.fits.psf.1.fits')
	                      del (imname // '.sub.2.fits.psg.1')
	                      del (imname // '.sub.2.fits.pst.1')
	                      del (imname // '.sub.2.fits.pst.2')        
	                      del (imname // '.als.1')
	                      del (imname // '.arj.1')        
	                      del (imname // '.sub.3.fits')
	                      del (imname // '.coo.1')
	                      del (imname // '.mag.1')
	                      del (imname // '.pst.1')
	                      del (imname // '.pstselect')
      del (imname // '.sub.3.fits.coo.1',verify=no,>>&"/dev/null")
	                      print ('', >> (imname//'.parameters'))
	                      print ('Re-do Daofind', >> (imname//'.parameters')) 
	                      print ('', >> (imname//'.parameters'))
	                 }
	                 else {
	                     check = no
	                 }
		                  }
		                  else {
		                    print (' Invalid choice. Try again.')
		                    check=no
		                  }
		              }
		          }

		      } # This bracket closes the 'if (l==1)' if

		      else {
		          print ('\n ----------------------------------------------------- ')
		          print (' Display task                                          ')
		          print (' Display .coo.1 (Frame1), .als.1 (Frame2) and .als.2 (Frame3) files')
		          print (' ----------------------------------------------------- ')
		          display ((imname // '.fit'), 1)
		          display ((imname // '.fit'), 2)
		          display ((imname // '.fit'), 3)

		          print ('\n ----------------------------------------------------- ')
		          print (' Tvmark                                                ')
		          print ('   - Frame1: detected stars by "daofind"               ')
		          print ('   - Frame2: detected stars by first "allstar" run     ')
		          print ('   - Frame3: detected stars by second "allstar" run    ')
		          print (' ----------------------------------------------------- ')
		          txdump.textfile = (imname // '.als.1')
		          txdump.headers = no
		          txdump.fields = 'xcenter, ycenter'
		          txdump.expr = 'yes'
		          txdump > auxiliar1 

		          txdump.textfile = (imname // '.als.2')
		          txdump.headers = no
		          txdump.fields = 'xcenter, ycenter'
		          txdump.expr = 'yes'
		          txdump > auxiliar2 
		          
		          tvmark.interactive = no
		          tvmark.mark = 'point'
		          tvmark.font = "raster"
		          tvmark.color = 204
		          tvmark.number = no
		          tvmark.label = no
		          tvmark.toleran = 3.5
		          tvmark (2, 'auxiliar1')
		        
		          tvmark.interactive = no
		          tvmark.mark = 'point'
		          tvmark.font = "raster"
		          tvmark.color = 204
		          tvmark.number = no
		          tvmark.label = no
		          tvmark.toleran = 3.5
		          tvmark (3, 'auxiliar2')

		          tvmark.interac = no 
		          tvmark.mark = 'point'
		          tvmark.font = "raster"
		          tvmark.color = 204
		          tvmark.number = no      
		          tvmark.label = no
		          tvmark.toleran = 3.5
		          tvmark (1, imname // '.coo.1')

		          del ('auxiliar1')
		          del ('auxiliar2')

		          print ('\n Which file to keep:')
		          print ('     [1]: keep .als.1,')
		          print ('     [2]: keep .als.2,')
		          print ('     [3]: Start over from \'Daofind\' search')
		          print ('          (the very beginning!!)')

		          check=no
		          while (check == no) {

                  var2=2
		              scan (var2)
		              if (var2==1) {
		                  check=yes
		                  del (imname // '.als.2')
		                  l=1
		                  print('\n Deleting '//imname//'.als.2 file')
		                  psf_done=yes
		              }
		              else {
		                  if (var2==2) {
		                  
	# ------------------------------------------------------------------------------------
	# Search for imname // '.als.2' file

		          o=0
		          files ('*'//imname //'.als.2', > 'tempfile')
		          flist = 'tempfile'
		          while (fscan (flist,line) != EOF) {
		              o = o + 1
		          }
		          del ("tempfile")
		          # ------------------------------------------------------------------------------------

		          if (o!=1) {
		              print('')
		              print (' ****************************')
		              print ('   Missing '//imname //'.als.2 file.')
		              print (' ****************************')
		              check=no
		          }
		          else {
		              if (o==1) {
		                  print ('')
		                  print (imname//'.als.2 file found.')
                      check=yes
                      psf_done=yes	                  
		              }
		              else {
		                  print (' Unknown error. Check code.')
		                  check=no
		              }
		          }
	# ------------------------------------------------------------------------------------
		                  }
		                  else {
		                      if (var2==3) { # DELETE ALL FILES AND START FROM THE BEGINNING
	                      print ('\n Are you sure you want to delete all files')
	                      print (' and start over from the beginning? (y/n)')
	                      print (' default = no')
	                      check=no
	                      scan (check)
	                      if (check) {                       
	                          check=yes
	                          del (imname // '.psf.1.fits')
	                          del (imname // '.grp.1')
	                          del (imname // '.psg.1')
	                          del (imname // '.pst.2',verify=no,>>&"/dev/null")
	                          del (imname // '.nst.1')
	                          del (imname // '.nrj.1')
	                          del (imname // '.psf1')
	                          del (imname // '.sub.1.fits')
	                          del (imname // '.sub.2.fits')
	                          del (imname // '.sub.2.fits.psf.1.fits')
	                          del (imname // '.sub.2.fits.psg.1')
	                          del (imname // '.sub.2.fits.pst.1')
	                          del (imname // '.sub.2.fits.pst.2',verify=no,>>&"/dev/null")
	                          del (imname // '.pstselect')      
	                          del (imname // '.als.1')
	                          del (imname // '.arj.1')        
	                          del (imname // '.sub.3.fits',verify=no,>>&"/dev/null")
	                          del (imname // '.sub.3.fits.coo.1')
	                          del (imname // '.coo.1')
	                          del (imname // '.mag.1')
	                          del (imname // '.pst.1')
	                          del (imname // '.mag.2')
	                          del (imname // '.mag.3')
	                          del (imname // '.als.2')
	                          del (imname // '.arj.2')        
	                          del (imname // '.sub.3.fits') 
	                          del (imname // '.sub.4.fits')
	                          print ('', >> (imname//'.parameters'))
	                          print ('Re-do Daofind', >> (imname//'.parameters')) 
	                          print ('', >> (imname//'.parameters'))
	                      }
	                      else {
	                          check=no
	                      }
		                      }
		                      else {
		                          print (' Invalid choice. Try again.')
		                          check=no
		                      }
		                  }
		              }
		          } # This bracket closes the 'while (check == no)' while
		      }  # This bracket closes the 'if (l==1)' else

          doapert = yes
      } # THIS BRACKET CLOSES THE 'dopsf' 'if' THAT OPENS ABOVE
      else {
          doapert = no   # This parameter indicates that no PSF photometry has been performed
                         # and thus no '.als' file exists, rather a '.mag' file exists.
          psf_done = yes # This is so the script will exit the 'psf_done' 'while' that opens before the
      }                  # firts 'daofind'
      
# ------------------------------------------------------------------------------
# End of 'Check .coo.1, .als.1 and .als.2 (if such file exists) files'
# -----------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Closing bracket
# -----------------------------------------------------------------------------
  } # THIS BRACKET CLOSES THE 'psf_done' 'while' THAT OPENS
    # BEFORE THE FIRST DAOFIND
# ------------------------------------------------------------------------------
# End of 'Closing bracket'
# -----------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Calculation of Aperture Correction
# ------------------------------------------------------------------------------

	apert_warning = no
	print ('\n Aperture correction')
	print ('\n ----------------------------------------------------- ')
	print (' Phot task (third run; second if no new stars added after PSF fitting)')
	print (' Computation of the aperture correction                ')
	print (' ----------------------------------------------------- ')
	photpars.aperture = (aperture // ',' // apstd)
	fitskypars.annulus = 2*fannulus  # In Massey-Davis guide, twice the value used first (ie: 2*fannulus)
	fitskypars.dannulus = fdannulu/2  # In Massey-Davis guide, half the value used first (ie: fdannulu/2)

      apert_null=yes
      while (apert_null==yes) {

          print ('')
		print ('------------------------------------------------', >> (imname//'.parameters')) 
		print ('Aperture correction', >> (imname//'.parameters'))
		print ('', >> (imname//'.parameters'))
		print ('Third Phot task (second if no new stars added after PSF fitting), Fitskypars.annulus = 2*fannulus = ', 2*fannulus, >> (imname//'.parameters'))
		print ('Third Phot task (second if no new stars added after PSF fitting), Fitskypars.dannulus = fdannulu/2 = ', fdannulu/2, >> (imname//'.parameters'))
		phot.verif = no
		phot ((imname // '.sub.2.fits'), (imname // '.psf1'), (imname // '.res'))
		txdump.textfile = (imname // '.res')
		txdump.headers = no
		txdump.fields = 'mag'
		txdump.expr = 'yes'
		txdump > auxiliar
		print ('\nMagnitudes and aperture corrections', >> (imname//'.parameters'))
		print ('Aperture Potometry: r1 = ' // aperture // '; Standar stars aperture: r2 = ' // apstd, >> (imname//'.parameters'))
		print ('', >> (imname//'.parameters'))
		file_var = 'auxiliar'
		avg = 0
		max = -1000
		min = 1000
		i = 0
		while (fscan (file_var, m1, m2) != EOF) {
		   if (m1 != 'INDEF' && m2 != 'INDEF') {
		      mag1 = real (m1)
		      mag2 = real (m2)
		      delta_m = mag2 - mag1
		      avg += delta_m
		         if (delta_m > max) {
		            max = delta_m
		         }
          if (delta_m < min) {
             min = delta_m
          }
		      i += 1
		      print (mag1, mag2, delta_m, >> (imname//'_aperture'))
		   }
		}
		delete ('auxiliar')
		print ('------------------------------------------------', >> (imname//'_aperture'))
		      apert_warning = no
          print ('Stars used to obtain aperture correction: ', i, >> (imname//'.parameters')) 
		if (i > 5) {
	apert_null=no
		   avg = (avg - max - min) / (i - 2) # To the total average 'avg' I substract the maximum and the minimum and 
		                                     # then I perform the average (ie, divide by the number of 'delta_m' used)
             print ('Aperture correction maxmin average:', avg, >> (imname//'.parameters'))
		}
		else if (i != 0) {
		    avg = avg / i
		    print ('\n Aperture correction was obtained with '//i//' PSF stars')
		    print ('\n Do you want to move on with an aperture correction of '//avg//'')
		    print (' or do you want to decrease \'datamin\' value by 10% and repeat \'phot\' task?')
		    print ('\n Repeat \'phot\' task = yes')
		    print (' Move on = no')
		    print ('\n default = yes')
		    check = yes
		    scan (check)
		    if (check) {
		        print ('Aperture correction: ', avg, >> (imname//'.parameters'))
    		    print ('Decreasing \'datamin\' value by 10% and repeating \'phot\' task', >> (imname//'.parameters'))
	    delete (imname//'.res',verify=no,>>&"/dev/null")
	    if (datamin>0.) {
	datamin = datamin - (0.1*datamin)
	    }
	    else {
	datamin = datamin + (0.1*datamin)
	    }
	    datapars.datamin = datamin
	    print ('\n Datamin = '//datamin)
	    print ('\n New datamin value = ', datamin, >> (imname//'.parameters'))
		    }
		    else {
		        apert_null=no # Exits the 'apert_null' 'while'
		        print ('Using less than 5 stars to perform aperture correction', >> (imname//'.parameters'))
		        print ('Aperture correction average:', avg, >> (imname//'.parameters'))
		    }
		}
		else {
		   avg = 0
		   print ('Not enough stars to compute an aperture correction, aperture correction = 0.', >> (imname//'.parameters'))
		   print ('Not enough stars to compute an aperture correction, aperture correction = 0.', >> (imname//'_aperture'))
		   print ('\n WARNING: Not enough stars to compute an aperture correction')
		   print ('\n Aperture correction set to 0.')
		   print ('')
		   apert_warning = yes
		}
		
		if (avg==0) {
		    print ('\n Aperture correction was not obtained since no suitable PSF stars')
		    print (' could be found (they all had MAG = INDEF)')
		    print ('\n Do you want to move on with an aperture correction of 0 or do you')
		    print (' want to decrease \'datamin\' value by 10% and repeat \'phot\' task?')
		    print ('\n Repeat \'phot\' task = yes')
		    print (' Move on = no')
		    print ('\n default = yes')
		    check = yes
		    scan (check)
		    if (check) {
		        print ('Aperture correction: ', avg, >> (imname//'.parameters'))
    		    print ('Decreasing \'datamin\' value by 10% and repeating \'phot\' task', >> (imname//'.parameters'))		    
	    delete (imname//'.res',verify=no,>>&"/dev/null")
	    if (datamin>0.) {
	datamin = datamin - (0.1*datamin)
	    }
	    else {
	datamin = datamin + (0.1*datamin)
	    }
	    datapars.datamin = datamin
	    print ('\n Datamin = '//datamin)
	    print ('\n New datamin value = ', datamin, >> (imname//'.parameters'))
		    }
		    else {
    		    print ('Moving on with null aperture correction', >> (imname//'.parameters'))		    
		        print ('Aperture correction: ', avg, >> (imname//'.parameters'))
		        apert_null=no # Exits the 'apert_null' 'while'
		    }
		}
	} # This bracket closes the 'apert_null' 'while'

	print ('', >> (imname//'_aperture'))
	print ('Aperture correction: ', avg, >> (imname//'_aperture'))
#       printf ("%-10.3f", avg, >> (imname//'.parameters'))
	print (imname//' ', avg, >> "aperture") # This file will contain the 'aperture' values of all the frames.
	print (' ', >> (imname//'.parameters'))
	print ('End of aperture correction', >> (imname//'.parameters'))
	print ('End of aperture correction', >> (imname//'_aperture'))

      if (doapert == yes) {
		      print ('\n Performing correction...')
		      if (l == 1) {
		          pcalc.mode = "hl"
		          pcalc.infile = imname//'.als.1'
		          pcalc.field = "MAG"
		          pcalc.value = 'MAG+'//avg
		          pcalc
		      }                            # This task corrects the 'MAG' value in the imname.als.*
		      else {                       # file by adding the previously calculated aperture (avg)
		          pcalc.mode = "hl"
		          pcalc.infile = imname//'.als.2'
		          pcalc.field = "MAG"
		          pcalc.value = 'MAG+'//avg
		          pcalc
		      }
		      print ('\n Aperture correction performed.')
		      print ('Aperture correction performed.', >> (imname//'_aperture'))
     print (' ', >> (imname//'.parameters'))
		      print ('Aperture correction performed.', >> (imname//'.parameters'))
     print ('------------------------------------------------', >> (imname//'.parameters'))
      }
      else {
		      print ('\n Performing correction...')
          pcalc.mode = "hl"
          pcalc.infile = imname//'.mag.1'
          pcalc.field = "MAG"
          pcalc.value = 'MAG+'//avg                                # This task corrects the 'MAG' value in the 'imname.mag.1'
          pcalc                                                    # file by adding the previously calculated aperture (avg)
		      print ('\n Aperture correction performed.')
		      print ('Aperture correction performed.', >> (imname//'_aperture'))          
      }		      

# ------------------------------------------------------------------------------
# End of 'Aperture Correction'
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
      # This file ("list") will be used by the 'daom.cl' script later on
      #  
      if (l==1) {
          print (imname//'.als.1', >> "list")
      }
      else {
          if (l==2) {
		          print (imname//'.als.2', >> "list")
          }
          else {
              if (l==3) {
                  print (imname//'.mag.1', >> "list")
              }
              else {
                  print ('\n Unknown error 2800. Check code.')
              }
          }
      }
       
      # Delete following files to save disk space
      print ('\n Delete '//imname//'.sub.1.fits file (the one')
      print (' with the PSF stars removed)? (y/n)')
      print (' default = yes')
      check = yes
      scan (check)
      if (check) {
          del ((imname) // '.sub.1.fits')
      }
      print ('\n Delete '//imname//'.sub.3.fits file (the one')
      print (' with ALL the stars removed)? (y/n)')
      print (' default = yes')
      check = yes
      scan (check)
      if (check) {
          del ((imname) // '.sub.3.fits')
      }      
      del ((imname) // '.sub.2.fits')
      if (l==2) {
          del ((imname) // '.sub.4.fits')
      }

# ------------------------------------------------------------------------------

	if (filelist_check) {     # IF THE INPUT IS A LIST OF FILES, THEN RERUN THE SCRIPT FOR THE NEXT
    	    goto rerun            # IMAGE IN LINE.
	}

      if (sigma_warning == no && apert_warning == no) {
		      print ('\n ----------------------------------------------------- ')
		      print (' Script "psfphot" finished correctly.                ')
      }
      else if (sigma_warning == no && apert_warning == yes) {
		      print ('\n ----------------------------------------------------- ')
		      print (' Script "psfphot" finished with a WARNING:           ')
		      print ('\n WARNING: Not enough stars to compute an aperture correction')
      }
      else if (sigma_warning == yes && apert_warning == no) {
		      print ('\n ----------------------------------------------------- ')
		      print (' Script "psfphot" finished with a WARNING:           ')
          print ('\n WARNING: sigma and sstd [SSTDEV] values differed significantly')
          print (' sigma = '//sigma//' ; sstd = '//sstd)
      }
      else {
		      print ('\n ----------------------------------------------------- ')
		      print (' Script "psfphot" finished with two WARNINGS:        ')
          print ('\n WARNING: sigma and sstd [SSTDEV] values differed significantly')
          print (' sigma = '//sigma//' ; sstd = '//sstd)
		      print ('\n WARNING: Not enough stars to compute an aperture correction')		          
      } 
      
#      print ('\n After processing ALL THE CLUSTERs FRAMES:')
#      print ('\n Move on to the next script: "daom"')
#      print ('\n Remember this last script must be executed inside the package:')
#      print ('\n             noao/digiphot/photcal                     ')
#      print (' ----------------------------------------------------- ')    

end

