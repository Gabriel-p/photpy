

procedure ctio ()

struct *flist {mode="h"}
struct *flist2 {mode="h"}
struct *file_var {mode="h"}


begin

      struct name, line3
      real line, boxsize, bigbox, niterate, tolerance
      real line2
      real ut
      real filter
      bool check, imalign, doalign
      string reference
      int k, m

      if (! defpac ("daophot")) {
          print ('')
          print (' This script must be loaded inside the package \'noao.digiphot.daophot\'')
          bye()
      }

      astutil  

      rename.files = '*.fits'
      rename.newname = 'fit'
      rename.field = 'extn'
      rename

      files ("*.fit", >> 'list')
      
      setairmass.equinox = 'equinox'
      setairmass.st = 'LST'
      setairmass.ut = 'time-obs'
      setairmass.date = 'date-obs'
      setairmass.update = yes
      setairmass ('*.fit')

      hedit (images="*.fit", fields="GAIN", value=1.44, add=no, addonly=yes, delete=no, verify=no, update=yes)
      hedit (images="*.fit", fields="RDNOISE", value=7., add=no, addonly=yes, delete=no, verify=no, update=yes)
      
      flist = 'list'
      while  (fscan (flist,name) != EOF) {
          print ('')
          print (' Editing '//name//' file...')
          print ('')
		      hselect.mode = "hl"
		      hselect.images = name
		      hselect.fields = "FILTER"
		      hselect.expr = yes
		      hselect > "tempctio"
			    flist2 = 'tempctio'
			    while (fscan (flist2,line) != EOF)
			    del ("tempctio")
			    filter = line
			    
          if (filter == 2) {
              hedit (images=(name), fields="FILTERS", value="B", add=yes, addonly=no, delete=no, verify=no, update=yes)
          }
          else {
              if (filter == 3) {
						      hedit (images=(name), fields="FILTERS", value="V", add=no, addonly=yes, delete=no, verify=no, update=yes)
              }
              else {
                  if (filter == 5) {
								      hedit (images=(name), fields="FILTERS", value="I", add=no, addonly=yes, delete=no, verify=no, update=yes)
                  }
                  else {
                      if (filter == 6) {
										      hedit (images=(name), fields="FILTERS", value="U", add=no, addonly=yes, delete=no, verify=no, update=yes)
                      }
                      else {
                          print ('ERROR IN NAME OF FILE: '//name)
                          print ('Filter value found: '//filter)
                          print ('Closing...')
                          bye()
                      }
                  }
              }
          }    			    
			    
		      
		      hselect.mode = "hl"
		      hselect.images = name
		      hselect.fields = "TIME-OBS"
		      hselect.expr = yes
		      hselect > "tempctio"
			    flist2 = 'tempctio'
			    while (fscan (flist2,line2) != EOF)
			    del ("tempctio")
			    ut = line2
		      hedit (images=(name), fields="UT", value=(ut), add=no, addonly=yes, delete=no, verify=no, update=yes)		      
      }


      print ('')
      print (' Trim images (this will trim ALL frames to')
      print (' the section: [26:4052,26:4052])? (y/n)')
      print (' default = no')
      check = no
      scan (check)
      if (check) {
          imcopy.input = '*.fit[26:4052,26:4052]'
          imcopy.output = '*.fit'
          imcopy.verbose = yes
          imcopy
      }  
      print ('')
      print (' Trimming done')

      
      print ('')
      print (' Perform transposing? (y/n)')
      print (' default = no')
      check = no
      scan (check)
      if (check) {
      
		      print ('')
		      print (' First transposing (this could take a while)...')
		      imtranspose (input="*.fit[-*,*]", output="*.fit")
		      
		      print ('')
		      print (' Second transposing (this could take a while)...')      
		      imtranspose (input="*.fit[-*,*]", output="*.fit")
		      
          daophot
      }
      

      print ('')
      print (' Perform align? (y/n)')
      print (' default = yes')
      check = yes
      scan (check)
      if (check) {
      
          findpars.threshold = 500000
          imalign = yes
          while (imalign) {
          
              flist = 'list'
              while  (fscan (flist,name) != EOF) {
                  print (name)
                  k = strlen(name)
                  if (substr (name, k-3, k) == ".fit") {
                      name = substr (name, 1, k-4)
                      print (name//'.fits', >> 'listfits')
                  }
              }
              print ('')
              print (' Input the \'reference\' image')
              scan (reference)
          
              print ('                                                       ')
              print (' ----------------------------------------------------- ')
              print (' Display task:                                         ')
              print (' Original Image: ' // reference                         )
              print (' ----------------------------------------------------- ')
              display ((reference), 1)

              print ('                                                       ')
              print (' ----------------------------------------------------- ')
              print (' Daofind task:                                         ')
              print (' Search for stars                                      ')
              print (' ----------------------------------------------------- ')
              daofind.verif = no
              daofind.verb = no
              daofind ((reference), "default")

              print ('                                                       ')
              print (' ----------------------------------------------------- ')
              print (' Tvmark                                                ')
              print ('   - circle = detected star                            ')
              print (' ----------------------------------------------------- ')
              print ('')

				      file_var = (reference // '.coo.'// 1)
				      m=0
				      while (fscan (file_var,line3) != EOF) {
				          m = m + 1
				      }              
              print ('')
              print (' Number of stars marked: '//(m-41))
              print ('')
              print (' We set a really high threshold (500000) so that NO stars are')
              print (' marked and the user selects the stars by hand.')
              print ('')
              print (' Continue: press \'y\' or \'n\' key (you will enter the interactive \'tvmark\')')
					    check = yes
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
              tvmark (1, (reference // '.coo.'// 1))
      
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
              tvmark (1, (reference // '.coo.'// 1))

              print ('')
              print (' Were stars correctly marked (else delete')
              print (' "*.coo.1" file and repeat)? (y/n)')
              print ('')
				      print (' default = yes')
				      check = yes
              scan (check)

              if (check == no) {
                  delete (reference // '.coo.'// 1)
                  imalign = yes
              }
              else {
                  imalign = no
                  
                  unlearn phot
						      unlearn centerpars
						      unlearn fitskypars
						      unlearn photpars                  
						      phot.interactive = no
						      phot.radplots = no
						      phot.update = yes
						      phot.verbose = no
						      phot.mode = 'hl'
						      phot.verify = no
                  phot (reference, "default", "default")
                
						      txdump.mode = 'hl' 
						      txdump.textfile = (reference//'.mag.1')
						      txdump.headers = no
						      txdump.fields = 'XCENTER,YCENTER'
						      txdump.expr = 'MAG[1]!=INDEF'
						      txdump > referimalign
              }          
          
          }
          
          doalign = yes
          while (doalign) {
          
		          print ('')
		          print (' Input \'imalign\' parameter values:')
		          print ('')
		          boxsize = 7
		          print (' Actual boxsize value = '//boxsize//' Input new value:')
		          scan (boxsize)
		          bigbox = 11
		          print (' Actual bigbox value = '//bigbox//' Input new value:')
		          scan (bigbox)
		          niterate = 8
		          print (' Actual niterate value = '//niterate//' Input new value:')
		          scan (niterate)
		          tolerance = 10
		          print (' Actual tolerance value = '//tolerance//' Input new value:')
		          scan (tolerance) 
		          
		          print ('')
		          print (' Performing \'imalign\'...')         
		          
		          imalign.input = '@list'
		          imalign.referenc = reference
		          imalign.coords = 'referimalign'
		          imalign.output = '@listfits' 
		          imalign.boxsize = boxsize
		          imalign.bigbox = bigbox
		          imalign.negative = no
		          imalign.niterate = niterate
		          imalign.tolerance = tolerance
		          imalign.maxshift = INDEF
		          imalign.shiftimages = yes
		          imalign.interp_type = "nearest"
		          imalign.boundary_typ = "nearest"
		          imalign.trimimages = yes
		          imalign.verbose = yes
		          imalign.mode = "hl"   
		          imalign

              delete ('@list')

		          print ('')
		          print (' Displaying shifted images')
		          print ('')
		          m=1
		          flist = ('listfits')
		          while (fscan(flist, name) != EOF) {
		              display (name, m)
		              m=m+1
		          }
		          
		          print ('')
		          print (' Were the images shifted correctly (else input new parameters and re-do)? (y/n)')
				      print (' default = yes')
		          check=yes
		          scan (check)
		          if (check) {
		              doalign = no
		          }
				      rename.files = '*.fits'
				      rename.newname = 'fit'
				      rename.field = 'extn'
				      rename  	
          }
          
      delete ('listfits')
      delete ('referimalign')
      delete (reference//'.coo.1')
      delete (reference//'.mag.1') 
      }

      delete ('list')      
      print ('')
      print (' Done!')

end


