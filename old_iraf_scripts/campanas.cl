

procedure campanas ()

struct *flist {mode="h"}
struct *flist2 {mode="h"}
struct *file_var {mode="h"}


begin

      struct name, line3, line2, ut
      real boxsize, bigbox, niterate, tolerance
      string line, filter, reference
      bool check, imalign, doalign
      int k, m, o

      if (! defpac ("daophot")) {
          print ('')
          print (' This script must be loaded inside the package \'noao.digiphot.daophot\'')
          bye()
      }

      astutil  

      rename.files = '*.fits'
      rename.newname = 'fit'
      rename.field = 'extn'
      rename.mode = "hl"
      rename

      files ("*.fit", >> 'list')
      
      hedit (images="*.fit", fields="RDNOISE", value=6.6, add=no, addonly=yes, delete=no, verify=no, update=yes)
      
      flist = 'list'
      while  (fscan (flist,name) != EOF) {
          print ('')
          print (' Editing '//name//' file...')
          print ('')
		      hselect.mode = "hl"
		      hselect.images = name
		      hselect.fields = "FILTER"
		      hselect.expr = yes
		      hselect > "tempcamp"
			    flist2 = 'tempcamp'
			    while (fscan (flist2,line) != EOF)
			    del ("tempcamp")
			    filter = line
			    
          if (filter == "B") {
              hedit (images=(name), fields="FILTERS", value="B", add=yes, addonly=no, delete=no, verify=no, update=yes)
          }
          else {
              if (filter == "V") {
						      hedit (images=(name), fields="FILTERS", value="V", add=no, addonly=yes, delete=no, verify=no, update=yes)
              }
              else {
                  if (filter == "I") {
								      hedit (images=(name), fields="FILTERS", value="I", add=no, addonly=yes, delete=no, verify=no, update=yes)
                  }
                  else {
                      if (filter == "U") {
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
		      hselect.fields = "UTSTART"
		      hselect.expr = yes
		      hselect > "tempcamp"
			    flist2 = 'tempcamp'
			    while (fscan (flist2,line2) != EOF)
			    del ("tempcamp")
			    ut = line2
		      hedit (images=(name), fields="UT", value=(ut), add=no, addonly=yes, delete=no, verify=no, update=yes)		      
      }

      
      print ('')
      print (' Perform transposing? (y/n)')
      print (' default = yes')
      check = yes
      scan (check)
      if (check) {
      
		      print ('')
		      print (' First transposing (this could take a while)...')
		      imtranspose (input="*.fit[*,-*]", output="*.fit")
		      
		      print ('')
		      print (' Second transposing (this could take a while)...')      
		      imtranspose (input="*.fit[*,*]", output="*.fit")
		      
          daophot
      }
      

      print ('')
      print (' Perform align? (y/n)')
      print (' default = yes')
      check = yes
      scan (check)
      if (check) {
      
          datapars.sigma = 1
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

							# ------------------------------------------------------------------------------------
							# Search for 'referimalign' file
		          o=0
		          files ('*referimalign', > 'tempfile')
		          flist = 'tempfile'
		          while (fscan (flist,line2) != EOF) {
		              o = o + 1
		          }
		          del ("tempfile")

		          if (o!=1) {
						      file_var = (reference // '.coo.'// 1)
						      m=0
						      while (fscan (file_var,line3) != EOF) {
						          m = m + 1
						      }
						      m = m-41
		          }
		          else {
		              if (o==1) {
				              tvmark.interactive = no
				              tvmark.outimage = ""
				              tvmark.mark = 'circle'
				              tvmark.font = "raster"
				              tvmark.txsize = 2
				              tvmark.radii = 10
				              tvmark.color = 204
				              tvmark.number = yes
				              tvmark.label = no
				              tvmark (1, 'referimalign')
				              
								      file_var = ('referimalign')
								      while (fscan (file_var,line3) != EOF) {
								          m = m + 1
								      } 	
		              }
		              else {
		                  print (' Unknown error. Check code.')
		                  bye()
		              }
		          }
							# End of 'Search for 'referimalign' file'
							# ------------------------------------------------------------------------------------

              print ('')
              print (' Number of stars marked: '//(m))
              print ('')
              print (' We set a really high threshold (500000) so that NO stars are')
              print (' marked and the user selects the stars by hand.')
              print ('\n If this is NOT the first run, then '// (m) //' stars will be marked')
              print ('\n MARK MORE THAN 5 STARS')
              print ('\n Continue: press \'return\' key (you will enter the interactive \'tvmark\')')
              scan (check)
              
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
              print (' default = yes')
              check = yes
              scan (check)

              if (check == no) {
                  delete (reference // '.coo.'// 1)
                  imalign = yes
              }
              else {
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
						      txdump >> referimalign
						      
						      file_var = ('referimalign')
						      m=0
						      while (fscan (file_var,line3) != EOF) {
						          m = m + 1
						      }              
		              print ('\n Number of stars with MAG[1]!=INDEF  = '//(m))

		              if (m < 5) {
		                  imalign = yes
		                  print (' Too few stars. Select new ones.')
		                  print ('')
		                  delete (reference // '.coo.'// 1)
		                  delete (reference // '.mag.'// 1)
		              }
		              else {
		                  imalign = no
		              }

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
		          print (' Use \'shift\' file? (y/n)')
		          print (' (the file MUST be named \'shift\')')
		          print ('\n One image per line with the x and y shifts')
		          print (' in columns  one  and  two:')
		          print (' eg: imname x_shift y_shift')
		          print ('\n Shifts are: x_shift=Xref-Xin')
		          print (' default = no')
		          check=no
		          scan (check)
		          if (check) {
		              imalign.shifts = 'shift'
		          }
		          else {
		              imalign.shifts = ''
		          }    

		          
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

              print ('\n Display images? (y/n)')
              print (' default=yes')
              check=yes
              scan (check)
              if (check) {
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
				      }
				      else {
				          doalign = no
				      }
				      rename.files = '*.fits'
				      rename.newname = 'fit'
				      rename.field = 'extn'
				      rename.mode = "hl"
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


