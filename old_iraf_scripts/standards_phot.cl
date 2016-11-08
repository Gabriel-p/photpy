################################################################################
#
#        ===========================================================
#            PERFORM DAOFIND AND PHOT TASK ON STANDARD FRAMES
#        ===========================================================
# 
#     To run properly this procedure must be executed inside the package:
#                           'noao.digiphot.daophot'
#
#                          Gabriel Perren 2009-2016
#
################################################################################ 

procedure standards_phot (fwhm, sky_std, aperture, datamax)

real fwhm  {prompt = "Typical FWHM for stellar images"}
real sky_std {prompt = "Sky standard desviation value [STDDEV]"}
real datamax {prompt = "Maximum good data value?"}
real aperture {prompt = "[aprox 4xFWHM] Aperture photometry radius [pixels]"}
struct *file_var {mode="h", prompt = "Internal file name variable"}
struct *flist {mode="h", prompt = "Internal file name variable"}
struct *flist2 {mode="h", prompt = "Internal file name variable"}


begin

      string imname, file_name, coofile
      struct line, line2, name
      real fitrad, sstd, cbox, radius, apert2, thresh, dmax
      bool check, check2
      int a, b, m

      fitrad = fwhm
      sstd = sky_std
      apert2 = aperture
      dmax = datamax
      
      print ('', >> "standards_phot_data")
      print ('Fitrad: ', fitrad, >> "standards_phot_data")
      print ('STDDEV: ', sstd, >> "standards_phot_data")
      print ('Aperture: ', apert2, >> "standards_phot_data")
      print ('Datamax: ', dmax, >> "standards_phot_data")

# ------------------------------------------------------------------------------
# Control 0
    
	    if (! defpac ("daophot")) {
	        print ('')
	        print (' This script must be loaded inside the package noao/digiphot/daophot')
	        bye()
	    }
	    else { # Do nothing
	    }
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Control 1

	    print ('')
	    print (' IRAF must be located inside the working folder')
	    print (' (where the standards frames are)')
	    print ('')
	    print (' ONLY THE STANDARS FRAMES should be inside this folder!')
	    print ('')
	    print (' Is this correct? (y/n)')
	    check = no
	    scan (check)
	    if (check == yes) { # Do nothing
	    }
	    else {
	        bye()
	    }
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Control 2

	    print ('')
	    print (' Only ONE frame per filter should be inside the')
	    print (' working folder. This frame should show ALL the standar')
	    print (' stars and NOT have saturated standar stars.')
	    print ('')
	    print (' In case one single frame can not show all standars and NOT')
	    print (' have saturated standars at the same time, then two or more')
	    print (' frames can be used, REMEMBERING in the next script (while creating')
	    print (' the "imsets" file) to add an extra field, since there CAN NOT be')
	    print (' two frames with the same filter in one single field.')
	    print ('')
	    print (' Is this correct? (y/n)')
	    check = no
	    scan (check)
	    if (check) { # Do nothing
	    }
	    else {
	        bye()
	    }
# ------------------------------------------------------------------------------------


    
      # ------------------------------------------------------------------------------
      # This value, 'apert2', will be used by the "casleopsf" script later on
      #
      delete ('apert_standard',verify=no,>>&"/dev/null")
      print (apert2, > "apert_standard")
      # ------------------------------------------------------------------------------

      files ('*cluster_folders', > 'temp')
      flist = 'temp'
      b=0
      while (fscan (flist,line) != EOF) {
          b = b + 1
      }
      delete ('temp')
      if (b>0) {
          flist = 'cluster_folders'
		      while (fscan(flist, line) != EOF) {              # A copy of the 'apert_standard' file is stored
		          cp ('apert_standard', ('../'//line//'/'))    # in each of the folders where star clusters are
		          
      delete ('../apert_standard',verify=no,>>&"/dev/null") # This statement delets the file created by this 'while'
                                                            # in the upper folder (this happens because the file 'cluster_folders'
                                                            # has an empty line at the end)		          
		      }
      }
      else {
          print ('')
          print (' WARNING: no \'cluster_folders\' file found')
          print (' File \'apert_standard\' has not been copied')
          print (' to stars frames folders, but to actual folder.')
          print ('')
          print (' Continue (else exit)? (y/n)')
          scan (check)
          if (check) {
          }
          else {
              bye()
          }    
      }

 
# ------------------------------------------------------------------------------------
# Daofind task 
# ------------------------------------------------------------------------------------

      # ------------------------------------------------------------------------------
      # We store the names of the standards frames in a file ("standars")
      # and store the name of the file inside the *struct variable 'flist'
      #
      files ("*.fit", > "standars")
      file_name = ("standars")
      flist = file_name
      # ------------------------------------------------------------------------------

      thresh = 50000
      # This while goes through ALL the standard stars frames
      #
      while (fscan (flist,line) != EOF) {

          k = strlen (line)
          if (substr (line, k-3, k) == ".fit") {
              imname = substr (line, 1, k-4)
              print ('')
              print (' Image file: '//imname//'.fit')
              sleep 2
          }
          else {
              print ('line = '//line)
              print ('substr = '//substr (line, k-3, k))
              print ('') 
              print (' Image file: '//imname//'.fit')
              print (' Incorrect image file')
              print ('Closing...')
              sleep 2
              bye()
          }

          # This while goes through a single frame all the times it is needed
          #
          check2 = no
#          a = 1 # This value increases to aviod an 'ERROR: Operation would overwrite existing file'
          while (check2 == no) {

              print ('\n Datapars and daopars tasks                            ')
              print (' Basic image parameters                                ')
              datapars.fwhmpsf = fitrad
              datapars.sigma = sstd
              datapars.datamin = INDEF
              datapars.datamax = dmax
              datapars.noise = 'poisson'
              datapars.ccdread = 'RDNOISE'
              datapars.gain = 'GAIN'
              datapars.exposur = 'EXPTIME'
              datapars.airmass = 'AIRMASS'
              datapars.filter = 'FILTERS'
              datapars.obstime = 'UT'
              findpars.threshold = thresh
              print (' FWHM = ' // fitrad // ' pixels                        ')
              print (' Threshold = '// thresh)
              print (' ----------------------------------------------------- ')
    
              print ('\n ----------------------------------------------------- ')
              print (' Display task:                                         ')
              print (' Original Image: ' // imname                            )
              print (' ----------------------------------------------------- ')
              display ((imname), 1)

              print (' Marking '//imname//' file')
              print (' Use existing .coo.* file? (y/n)')
              print (' default = yes')
              check = yes
              scan (check)
              
              if (check) {
                  print (' Input name of existing .coo.* file to use')
#                 print (' (WHITHOUT the .coo.* extension)')

                  files ('*.coo.1', > 'list')
		              flist2 = 'list'
		              m=0
		              while  (fscan (flist2,name) != EOF) {
		                  print (name)
		                  m=m+1
		              }    
                  del ('list')

                  if (m==0) {
                      print ('\n NO .coo.1 FILES FOUND. ANSWER \'no\' TO THE NEXT QUESTION')
                  }
                  else {
		                  scan (coofile)
		#                 coofile = substr (coofile, 1, k-6)
		                  copy (coofile, imname//'.coo.1')
				              tvmark.interactive = no
				              tvmark.outimage = ""
				              tvmark.mark = 'circle'
				              tvmark.font = "raster"
				              tvmark.txsize = 2
				              tvmark.radii = 10
				              tvmark.color = 204
				              tvmark.number = yes
				              tvmark.label = no
				#             tvmark (1, (imname // '.coo.'// a))
				              tvmark (1, (imname // '.coo.1'))
                  }

              }
              else {
              
		              print ('                                                       ')
		              print (' ----------------------------------------------------- ')
		              print (' Daofind task:                                         ')
		              print (' Search for stars                                      ')
		              print (' ----------------------------------------------------- ')
		              daofind.verif = no
		              daofind.verb = no
		              daofind ((imname), "default")

		              print ('                                                       ')
		              print (' ----------------------------------------------------- ')
		              print (' Tvmark                                                ')
		              print ('   - circle = detected star                            ')
		              print (' ----------------------------------------------------- ')
		              print ('')

		#				      file_var = (imname // '.coo.'// a)
						      file_var = (imname // '.coo.1')
						      m=0
						      while (fscan (file_var,line2) != EOF) {
						          m = m + 1
						      }              
		              print ('')
		              print (' Number of stars marked: '//(m-41))
		              print ('')
		              print (' We set a really high threshold (50000) so that NO stars are')
		              print (' marked and the user selects the standards by hand. Remember')
		              print (' that ALL standard stars should be marked (with the \'a\' key).')
		              print ('')
		              print (' If you want \'Daofind\' to search for stars automatically,')
		              print (' then quit the interactive \'tvmark\' with the \'q\' key, answer')
		              print (' \'no\' to the next question, decrease the "threshold" value and')
		              print (' the task will run again with this new threshold.')
		              print ('')
		              print (' After this, you will have to edit the \'.coo\' file so that ONLY')
		              print (' the standard stars remain in it, so if extra stars are being')
		              print (' marked, you can edit them out them later.') 
		              print ('')
		              print (' Continue: press \'y\' or \'n\' key (you will enter the interactive \'tvmark\')')
		              scan (check)
		              
		              tvmark.interactive = no
		              tvmark.outimage = ""
		              tvmark.mark = 'circle'
		              tvmark.font = "raster"
		              tvmark.txsize = 2
		              tvmark.radii = 10
		              tvmark.color = 204
		              tvmark.number = yes
		              tvmark.label = no
		#             tvmark (1, (imname // '.coo.'// a))
		              tvmark (1, (imname // '.coo.1'))
		      
		              # The first tvmark marks the stars in the image, found by 'daofind'
		              # The second one (below this), performs the 'interactive' tvmark.
		              # If I try to make the first tvmark 'interactive', then
		              # it doesn't mark the stars, this way the found star are marked

		              tvmark.interactive = yes 
		              tvmark.outimage = ""
		              tvmark.mark = 'circle'
		              tvmark.font = "raster"
		              tvmark.txsize = 2
		              tvmark.radii = 10
		              tvmark.color = 204
		              tvmark.number = yes
		              tvmark.label = no
		#             tvmark (1, (imname // '.coo.'// a))
		              tvmark (1, (imname // '.coo.1'))
              
              }

              print ('\n Are ALL the standard stars marked (else repeat) (y/n)')
#              print (' "threshold" value and repeat)? (y/n)')
              check2=yes
              print (' default = yes')
#              print ('')
              scan (check2)

              if (check2 == no) {
                  delete (imname // '.coo.1')
#                  print (' Renaming '//(imname // '.coo.'// a)//' file to '//(imname // '.coo.'// a+1))
#                  rename.files = (imname // '.coo.'// a)
#                  rename.newname = (imname // '.coo.'// a+1)
#                  rename.field = all
#                  print ('')
#                  print (' Actual threshold value: '//thresh)
#                  print ('')
#                  print (' Input new threshold value')
#                  print ('(a higher value means less stars found)')
#                  scan (thresh)
#                  a = a+1
              }
              else {
                  # We edit the .coo.* file so that only the standard stars remain
                  print ('\n Edit the following file so that ONLY the standard stars remain')
                  print ('\n Continue: press \'y\' or \'n\' key (you will enter \'edit\' mode)')
                  scan (check)
#                  vi (imname // '.coo.'// a)
                  vi (imname // '.coo.1')
              }
              
          } # This bracket closes the 'check2' while

          print (imname, >> 'standars2')

      } # This bracket closes the first while

# ------------------------------------------------------------------------------------
# End of 'Daofind task' 
# ------------------------------------------------------------------------------------


# ------------------------------------------------------------------------------------
# Phot task 
# ------------------------------------------------------------------------------------

      print ('                                                       ')
      print (' ----------------------------------------------------- ')
      print (' Phot task                                             ')
      print (' Aperture photometry (radius = '// apert2// ' pixels)  ')
      print (' ----------------------------------------------------- ')

      phot.verif = no
      phot.verb = yes
      phot.interac = no
      phot.radplots = no
      photpars.aperture = apert2
      centerpars.calgorithm = "centroid"      

      # Page 9 of 'A user's guide to stellar CCD photometry...' by Massey & Davis
      print (' If the FWHM of the frames is unusually large (>=4.), the size')
      print (' of "cbox" must be something like twice the FWHM ('//2*fitrad//' in this case)')
      print (' Otherwise a "cbox" value of 5 is fine')
      print (' Since the FWHM is '//fitrad//' it is recommended to use a "cbox" value of: ')
      if (fitrad >= 4.) {
          cbox = 2*fitrad
          print (' cbox = '//cbox)
      }
      else {
          cbox = 5.
          print (' cbox = '//cbox)
      }
      print ('')
      print (' Accept this "cbox" value (else input new one)? (y/n)')
      print (' default = yes')
      check = yes
      scan (check)
      if (check) {
      }      
      else {
          print (' Input new "cbox" value')
          scan (cbox)
      }    
      centerpars.cbox = cbox
      radius = apert2 + 5
      fitskypars.annulus = radius
      fitskypars.dannulu = 5
      phot ("@standars2", "default", "default")

# ------------------------------------------------------------------------------------
# End of 'Phot task' 
# ------------------------------------------------------------------------------------
      
      print ('')
      print (' Keep *.coo.* files? (y/n)')
      print (' default = yes')
      check = yes
      scan (check)
      if (check) {
      }
      else {
          del ('*.coo.*')
      }

      delete ("standars")

      print ('                                                       ')
      print (' ----------------------------------------------------- ')
      print (' Script "standards_phot" finished correctly.               ')
      print (' Move on to the "standards_cat" script                   ')
      print ('') 
      print (' Remember this last script must be executed inside the package:')
      print ('')
      print ('             noao/digiphot/photcal                     ')
      print (' ----------------------------------------------------- ')


end
