################################################################################
#
#        ===========================================================
#             PROCEDURE TO PERFORM MKCONFIG AND FITPARAMS TASK
#        ===========================================================
# 
#     To run properly this procedure must be executed inside the package:
#                           'noao.digiphot.photcal'
#
#                            by Gabriel Perren 2009
#
#
#    To run properly this procedure needs to be executed inside the folder where
#           the following files exist (image files are NOT necessary):
#
#             noche.obs, fnoche.obs, flandolt.cat.dat and landolt.dat
#
################################################################################         
      

procedure casleored4 ()      

struct *flist

begin

      struct line
      bool check, check2, auxi
      int o, b


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
    if (check) { # Do nothing
    }
    else {
        bye()
    }
# ------------------------------------------------------------------------------------
# Control 1
# ------------------------------------------------------------------------------------   
   
  

# ------------------------------------------------------------------------------------
# Mkconfig task 
# ------------------------------------------------------------------------------------  


check2=no
check=no
  
      print ('')
      print (' Do you want to use an existing "noche.cfg" file (else run \'mkconfig\' task to create one)? (y/n)')
      scan (check)
      if (check) {
      
          # ------------------------------------------------------------------------------------
          # Search for 'noche.cfg' file
          # ------------------------------------------------------------------------------------
          auxi=no
          while (auxi==no) { 

              o=0

              # Search for "noche.cfg"
              # ------------------------------------------------------------------------------------
              files ('*noche.cfg', > 'tempfile')
              flist = 'tempfile'
              while (fscan (flist,line) != EOF) {
              o = o + 1
              }
              del ("tempfile")
              # ------------------------------------------------------------------------------------

              if (o!=1) {
	                print('')
	                print (' ****************************')
	                print ('   Missing "noche.cfg" file.')
	                print (' ****************************')
	                print ('')
	                print (' Search for "noche.cfg" file again (else run')
	                print (' \'Mkconfig\' task to create one manually)? (y/n)')
	                scan (check)
	                if (check) {
	                    auxi=no
	                }
	                else {
	                    auxi=yes
	                    check2 = no # Move on to 'Mkconfig' task
	                }
              }
              else {
	                if (o==1) {
                      print ('')
                      print ('"noche.cfg" file found.')
                      auxi=yes
                      check2 = yes                      
	                }
	                else {
                      print (' Unknown error. Check code.')
                      bye()
	                }
              }
          }
          # ------------------------------------------------------------------------------------
          # End of 'Search for 'noche.cfg' file'
          # ------------------------------------------------------------------------------------ 
      }
       

      while (check2==no) {

            print (' ----------------------------------------------------- ')
            print ('                    Mkconfig task                      ')
            print ('')
            print (' When this task is over, edit the output file so the variables in')
            print (' the "transformation" section are equal to the ones listed in')
            print (' the "observations" section just above it.')
            print (' If the filters IDs were NOT changed, then this WONT be necessary.') 
            print ('')
            print (' Also, edit OUT the equations for the filters which were not used')
            print (' and edit in the right coefficients.')
            print ('')
            print (' CASLEO extinction coefficients:')
            print ('')
            print (' fit    u1=0.0, u3=0.000           fit    b1=0.0, b3=0.000')
            print (' const  u4=0.0, u2=0.49            const  b4=0.0, b2=0.27')
            print ('')
            print (' fit    v1=0.0, v3=0.000           fit    i1=0.0, i3=0.000')
            print (' const  v4=0.0, v2=0.12            const  i4=0.0, i2=0.02')
            print ('')
            print (' CTIO extinction coefficients:')
            print ('')
            print (' fit    u1=0.0, u3=0.000           fit    b1=0.0, b3=0.000')
            print (' const  u4=0.0, u2=0.4544678       const  b4=0.0, b2=0.2669824')
            print ('')
            print (' fit    v1=0.0, v3=0.000           fit    i1=0.0, i3=0.0')
            print (' const  v4=0.0, v2=0.1418861       const  i4=0.0, i2=0.1029978')
            print ('')
            print (' ----------------------------------------------------- ')
            print (' Continue, press \'y\' or \'n\' key...')
            scan (check)

            mkconfig.verbose = yes
            mkconfig ("noche.cfg", "flandolt.cat.dat", "noche.obs", "landolt")

            print ('')
            print (' If Warnings = 0 and Errors = 0, we are doing fine :)')
            print ('')
            print (' Continue with \'fitparams\' task (else, re-do Mkconfig task)? (y/n)')
            scan (check)
            if (check) {
                check2=yes
            }
            else {
                check2=no
                del ('noche.cfg') 
            }
      }
    
# ------------------------------------------------------------------------------------
# End of 'Mkconfig task'
# ------------------------------------------------------------------------------------ 



# ------------------------------------------------------------------------------------
# Fitparams task 
# ------------------------------------------------------------------------------------ 

      print ('                                                       ')
      print (' ----------------------------------------------------- ')
      print (' Fitparams task                                        ')
      print ('')
      print (' Check that the RMS is in the hundredth.')
      print (' The character "?" shows the help file')
      print (' ----------------------------------------------------- ')
      fitparams.weighting = "photometric"
      fitparams.addscatter = yes # "..Addscatter is recommended if weighting is "photometric"..." ; From 'fitparams' help
      fitparams ("noche.obs", "landolt.cat", "noche.cfg", "noche.ans")

# ------------------------------------------------------------------------------------
# End of 'Fitparams task'
# ------------------------------------------------------------------------------------


      files ('*cluster_folders', > 'temp')
      flist = 'temp'
      b=0
      while (fscan (flist,line) != EOF) {
          b = b + 1
      }
      delete ('temp')
      if (b>0) {
		      flist = ('cluster_folders')
		      while (fscan(flist, line) != EOF) {       # A copy of this files is stored
		          cp ('noche.ans', ('../'//line//'/'))  # in each of the folders where star clusters are
		          cp ('noche.cfg', ('../'//line//'/'))
		      }
		      delete ('../noche.ans',verify=no,>>&"/dev/null")
		      delete ('../noche.cfg',verify=no,>>&"/dev/null")
      }
      else {
          print ('')
          print (' WARNING: no \'cluster_folders\' file found')
          print (' Files \'noche.ans\' and \'noche.cfg\' were not copied')
          print (' to clusters folders, but to actual folder.')
          print ('')
      }



      print ('                                                       ')
      print (' ----------------------------------------------------- ')
      print (' Script "casleored4" ended.               ')
#      print ('')
#      print (' Before you move on to the "casleopsf" script')
#      print (' REMEMBER TO COPY THE FOLLOWING FILES TO THE FOLDER WHERE')
#      print (' THE STARS FRAMES ARE LOCATED:')
#      print ('')
#      print (' "noche.ans", "noche.cfg", "gain_rdnoise" and "apert_standard"')
      print ('')
      print (' Remember "casleopsf" must be executed inside the package:')
      print ('')
      print ('             noao/digiphot/daophot                     ')
      print (' ----------------------------------------------------- ')

end

