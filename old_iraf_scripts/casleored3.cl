################################################################################
#
#        ===========================================================
#             PROCEDURE TO PERFORM MKCATALOG AND MKNOBSFILE TASK
#        ===========================================================
# 
#     To run properly this procedure must be executed inside the package:
#                           'noao.digiphot.photcal'
#
#                            by Gabriel Perren 2009
#
################################################################################      
      

procedure casleored3 () 

struct *flist

begin

     struct line
     bool check, auxi, auxi2, mknobs
     real toleran
     int o
   
# ------------------------------------------------------------------------------------
# Control 0
# ------------------------------------------------------------------------------------
    
    if (! defpac ("photcal")) {
        print ('')
        print (' This script must be loaded inside the package noao/digiphot/photcal')
        bye()
    }
    else { # Do nothing
    }
# ------------------------------------------------------------------------------------
# Control 0
# ------------------------------------------------------------------------------------

   
# ------------------------------------------------------------------------------------
# Control 1
# ------------------------------------------------------------------------------------

    print ('')
    print (' IRAF must be located inside the working folder')
    print (' (where the standards frames are)')
    print ('')
    print (' Is IRAF located inside the working folder? (y/n)')
    scan (check)
    if (check == yes) { # Do nothing
    }
    else {
        bye()
    }
# ------------------------------------------------------------------------------------
# Control 1
# ------------------------------------------------------------------------------------


   
# ------------------------------------------------------------------------------------
# Mkcatalog task 
# ------------------------------------------------------------------------------------
 
      print ('')
      print (' Create new "landolt.cat" file (else use existing one)? (y/n)')
      print (' (".cat" file is created through "mkcatalog" task)')
      print ('')
      scan (check)
      if (check) {
          print ('                                                       ')
          print (' ----------------------------------------------------- ')
          print (' Mkcatalog task                                        ')
          print (' Catalog name = "landolt"                              ')
          print ('')
          print ('Ctrl + D equals <EOF>')
          print (' ----------------------------------------------------- ')
          mkcatalog.review = no
          mkcatalog.verify = no
          mkcatalog.edit = yes
          mkcatalog ("landolt.cat")
      }
      else {
      
          # ------------------------------------------------------------------------------------
          # Search for 'landolt.cat' file
          # ------------------------------------------------------------------------------------
                auxi=no
                while (auxi==no) { 

                    o=0

                    # Search for "landolt.cat"
                    # ------------------------------------------------------------------------------------
                    files ('*landolt.cat', > 'tempfile')
                    flist = 'tempfile'
                    while (fscan (flist,line) != EOF) {
                        o = o + 1
                    }
                    del ("tempfile")
                    # ------------------------------------------------------------------------------------

                    if (o!=1) {
                        print('')
                        print (' ****************************')
                        print ('   Missing "landolt.cat" file.')
                        print (' ****************************')
                        print ('')
                        print (' If you choosed to use an existing file, then')
                        print (' copy this file to the "standard/" folder before')
                        print (' moving on (this file won\'t be used until the next')
                        print (' script: "casleored4").') 
                        auxi=no
                        print (' Continue (else exit)? (y/n)')
                        scan (check)
                        if (check) {
                        }
                        else {
                            bye()
                        }
                    }
                    else {
                        if (o==1) {
                            print ('')
                            print ('"landolt.cat" file found.')
                            auxi=yes
                        }
                        else {
                            print (' Unknown error. Check code.')
                            bye()
                        }
                    }
                }

          # ------------------------------------------------------------------------------------
          # End of 'Search for 'landolt.cat' file'
          # ------------------------------------------------------------------------------------
                      
      }    

# ------------------------------------------------------------------------------------
# End of 'Mkcatalog task'
# ------------------------------------------------------------------------------------



# ------------------------------------------------------------------------------------
# Creation of the 'imsets' file
# ------------------------------------------------------------------------------------

      print ('')
      print (' Use existing \'sets\' file (else create one now)? (y/n)')
      scan (check)
      if (check) {
          # ------------------------------------------------------------------------------------
          # Search for 'sets' file
          # ------------------------------------------------------------------------------------
                auxi=no
                while (auxi==no) { 

                    o=0

                    # Search for "sets"
                    # ------------------------------------------------------------------------------------
                    files ('*sets', > 'tempfile')
                    flist = 'tempfile'
                    while (fscan (flist,line) != EOF) {
                        o = o + 1
                    }
                    del ("tempfile")
                    # ------------------------------------------------------------------------------------

                    if (o!=1) {
                        print('')
                        print (' ****************************')
                        print ('   Missing "sets" file.')
                        print (' ****************************')
                        print ('')
                        print (' If you choosed to use an existing \'sets\' file, then')
                        print (' copy this file to the current working folder now.')
                        print ('')
                        print (' Try again (else move on to create a \'sets\' file)? (y/n)')
                        scan (check)
                        if (check) {
                            auxi=no
                        }
                        else {
                            auxi2 = yes
                            auxi = yes
                        }
                    }
                    else {
                        if (o==1) {
                            print ('')
                            print ('"sets" file found.')
                            auxi=yes
                            auxi2 = no
                        }
                        else {
                            print (' Unknown error. Check code.')
                            bye()
                        }
                    }
                }

          # ------------------------------------------------------------------------------------
          # End of 'Search for 'sets' file'
          # ------------------------------------------------------------------------------------
      }
      else {
          auxi2 = yes
      }
      
      if (auxi2) {
		      copy ('estandars2', 'sets')
		      print ('')
		      print (' Edit the following file, which will be used as the "imsets" file')
		      print ('')
		      print (' Remember there can only be ONE frame per filter (if you used two')
		      print (' frames with the same filter, then you must add a new field since')
		      print (' there CAN NOT be two frames with the same filter in one single field)')
		      print ('')
		      print (' Example:')
		      print ('')
		      print (' p13 : p13ab120 p13ai90 p13au300 p13av120')
		      print (' p15 : pg15b90 pg15i120 pg15u300 pg15v90')
		      print ('')
		      print (' Continue? (y/n)')
		      check = no
		      while (check == no) {
		          scan (check)
		      }
		      vi ('sets')
      }

# ------------------------------------------------------------------------------------
# End of 'Creation of the 'imsets' file'
# ------------------------------------------------------------------------------------


# ------------------------------------------------------------------------------------
# Mknobsfile task 
# ------------------------------------------------------------------------------------

      toleran = 5.
      mknobs = no
      while (mknobs == no) {
          print ('                                                       ')
          print (' ----------------------------------------------------- ')
          print (' Mknobsfile task                                       ')
          print (' Current tolerance value = '//toleran                   )
          print ('')
          print (' ----------------------------------------------------- ')
          print ('')
          print (' The \'mknobsfile\' task will now create the \'noche.obs\' file.')
          print (' After this, the file will be displayed. If you see "INDEF" values,')
          print (' you can change the \'tolerance\' value and repeat in the next step.')
          print ('')
          check = yes
          print (' Continue (else exit)? (y/n)')
          scan (check)
          if (check) {
          }
          else {
              bye()
          }
          
          print ('')
          print (' Set \'mknobsfile.allfilters\' parameter (choose y to set the')
          print (' parameter to \'yes\' and n to set it to \'no\' )')
          print ('')
          print (' (Output only objects which are successfully matched in all the filters')
          print (' specified by \'idfilters\', else output objects regardless of whether')
          print (' they were found for all the filters)')
          scan (check)
          if (check) {
              mknobsfile.allfilters = yes
          }
          else {
              mknobsfile.allfilters = no
          }
          print ('')
          print (' Use existing shifts file (MUST be named \'shifts\')? (y/n)')
          scan (check)
          if (check) {
		          # ------------------------------------------------------------------------------------
		          # Search for 'shifts' file
		          # ------------------------------------------------------------------------------------
		                auxi=no
		                while (auxi==no) { 

		                    o=0

		                    # Search for "shifts"
		                    # ------------------------------------------------------------------------------------
		                    files ('*shifts', > 'tempfile')
		                    flist = 'tempfile'
		                    while (fscan (flist,line) != EOF) {
		                        o = o + 1
		                    }
		                    del ("tempfile")
		                    # ------------------------------------------------------------------------------------

		                    if (o!=1) {
		                        print('')
		                        print (' ****************************')
		                        print ('   Missing "shifts" file.')
		                        print (' ****************************')
		                        print ('')
		                        print (' If you choosed to use an existing \'shifts\' file, then')
		                        print (' copy this file to the current working folder now.')
		                        print ('')
		                        print (' Try again (else move on without using a \'shifts\' file)? (y/n)')
		                        scan (check)
		                        if (check) {
		                            auxi=no
		                        }
		                        else {
		                            mknobsfile.shifts = ''
		                            auxi = yes
		                        }
		                    }
		                    else {
		                        if (o==1) {
		                            print ('')
		                            print ('"shifts" file found.')
		                            mknobsfile.shifts = 'shifts'                            
		                            auxi=yes
		                        }
		                        else {
		                            print (' Unknown error. Check code.')
		                            bye()
		                        }
		                    }
		                }

		          # ------------------------------------------------------------------------------------
		          # End of 'Search for 'shifts' file'
		          # ------------------------------------------------------------------------------------           
          }
          else {
              mknobsfile.shifts = ''
          } 
          
          mknobsfile.tolerance = toleran
          mknobsfile.verify = no
          mknobsfile.verbose = yes
          mknobsfile ("*mag*", "U,B,V,I", "sets", "noche.obs")

          vi ('noche.obs')
          
          print ('')
          print ('Is the output file correct (no INDEFs)? (y/n)')
          print ('')
          scan (check)

          if (check == yes) {
              mknobs = yes
          }
          else {
              mknobs = no
              print ('')
              print (' Tolerance value = '//toleran )
              print (' Input new tolerance value (a higher tolerance might get rid of the INDEFs)')
              scan (toleran)
              del ('noche.obs')
          }

      }

# ------------------------------------------------------------------------------------
# End of 'Mknobsfile task'
# ------------------------------------------------------------------------------------

      print ('                                                       ')
      print (' ----------------------------------------------------- ')
      print (' Script "casleored3" finished correctly.               ')
      print ('')
      print (' ******************************************************')
      print ('')
      print ('                     IMPORTANT')
      print ('')
      print (' BEFORE you move on to the next script:')
      print ('')
      print (' Edit the "noche.obs" file so that the names of the standard')
      print (' stars matches those in the "landolt.cat" file.')
      print ('') 
      print (' ******************************************************')
      print ('')
      print (' Then move on to the "casleored4" script               ')
      print ('')
      print (' Remember this last script must be executed inside the package:')
      print ('')
      print ('             noao/digiphot/photcal                     ')
      print (' ----------------------------------------------------- ')

end
