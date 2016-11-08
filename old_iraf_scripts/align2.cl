################################################################################
#
#        ===========================================================
#                      PROCEDURE TO PERFORM IMALIGN (2)
#        ===========================================================
# 
#     To run properly this procedure must be executed inside the package:
#                          'images.immatch'
#
#                       by Gabriel Perren 2009
#
################################################################################


procedure align2 ()

struct *list {mode="h"}
struct *list2 {mode="h"}
struct *list3 {mode="h"}
struct *list4 {mode="h"}


begin


      int i, j, k, m, q, n
      struct name, name, line
      string var[50]
      bool check, auxi
      real boxsize, bigbox, niterate, tolerance, gain, rdnoise


# ------------------------------------------------------------------------------------
# Control 0
# ------------------------------------------------------------------------------------

      print ('')
      print (' This script must be loaded inside the package images/immatch')
      print (' Is this package loaded? (y/n)')
      check=no
      scan (check)
      if (check) {
      }
      else {
         bye()
      }

# ------------------------------------------------------------------------------------
# Control 0
# ------------------------------------------------------------------------------------   

# ------------------------------------------------------------------------------------  
      # Check the existence of the file 'gain_rdnoise' previously created by the
      # "casleored" script. If such file does not exist, input the values and
      # the file will be created.
      #
      files ('*gain_rdnoise', > 'temp.grd')
      list3 = 'temp.grd'
      n=0
      while (fscan (list3,line) != EOF) {
          n = n + 1
      }
      if (n == 0) {
          print ('')
          print (' No "gain_rdnoise" file.')
          print ('')
          print (' Input GAIN value')
          scan (gain)
          print (' Input RDNOISE value')
          scan (rdnoise)
          print ('')
          print (gain,rdnoise, > "gain_rdnoise") # This file will be used by the 'casleopsf' script later on
      }
      else { # The file exists, then do nothing.
      }
      del ('temp.grd')
# ------------------------------------------------------------------------------------  

      mkdir ('standard') # In this folder the standard stars sets will be stored later on

      list = "clusters"
      j=0
      for (i=1; fscan(list, name) != EOF; i=i+1) {
         j=j+1
         k = strlen(name)
         name2 = substr (name, k-3, k)
         var[i] = name2
      }

      m=1
      while (m<=j) {
      print (var[m], >> 'cluster_folders') # This file will be used by 'casleored2' script.
      list = ('cluster.'//var[m])

          while (fscan(list, name) != EOF) {

              k = strlen(name)
              name = substr (name, 1, k-3)

              print ((name//'.out'), >> ('cluster.'//var[m]//'.out'))
          }

      m=m+1
      }

      mv ('cluster_folders', 'standard/')

      m=1
      while (m<=j) {

          print ('')
          print (' Input parameter values:')
          print ('')
          boxsize = 11
          print (' Actual boxsize value = '//boxsize//' Input new value:')
          scan (boxsize)
          bigbox = 51
          print (' Actual bigbox value = '//bigbox//' Input new value:')
          scan (bigbox)
          niterate = 8
          print (' Actual niterate value = '//niterate//' Input new value:')
          scan (niterate)
          tolerance = 10
          print (' Actual tolerance value = '//tolerance//' Input new value:')
          scan (tolerance)

          auxi=no
          while (auxi==no) {

              print ('')
              print (var[m]//' Frames')
              print ('')
              list = ('cluster.'//var[m])
              while (fscan(list, name) != EOF) {
              print (name)
              }

              list = ('frame.'//var[m])
              while (fscan(list, name) != EOF) {
                  print (name)
              }

              print ('')
              imalign.input = ('@cluster.'//var[m])
              imalign.reference = (name//'.fit')
              imalign.coords = (name//'.coo.1')
              imalign.output = ('@cluster.'//var[m]//'.out')
              imalign.boxsize = boxsize
              imalign.bigbox = bigbox
              imalign.negative = no
              imalign.niterate = niterate
              imalign.tolerance = tolerance
              imalign.maxshift = INDEF
              imalign.shiftimages = yes
              imalign.interp_type = "linear"
              imalign.boundary_typ = "nearest"
              imalign.trimimages = yes
              imalign.verbose = yes
              imalign.mode = "hl"
              imalign

              print ('')
              print (' Displaying shifted images')
              print ('')
              q=1
              list = ('cluster.'//var[m]//'.out')
              while (fscan(list, name) != EOF) {
                  display (name, q)
                  q=q+1
              }

              print ('')
              print (' Were the images shifted correctly (else input new parameters and re-do)? (y/n)')
              check=yes
              scan (check)

              if (check) {
                  auxi=yes

                  del ('cluster.'//var[m])
                  del ('frame.'//var[m])
                  del (var[m]//'*.coo.1')

                  print ('')
                  print (' Is this a standard star set: '//var[m]//' ? (y/n)')
                  scan (check)
                  if (check) {
                      list2 = ('cluster.'//var[m]//'.out')
                      while (fscan(list2, name2) != EOF) {

                          k = strlen(name2)
                          name2 = substr (name2, 1, k-4)

                          del (name2) # Delete original image frame
                          rename.mode = "hl"
                          rename.field = "all"
                          rename.mode = "hl"
                          rename (name2//'.out.fits', name2) # Rename .out images to original names
                          mv (name2, 'standard/')            # Move out shifted images
                      }

                      del ('cluster.'//var[m]//'.out')  
                  }
                  else {
                      mkdir (var[m])
                      list2 = ('cluster.'//var[m]//'.out')
                      while (fscan(list2, name2) != EOF) {

                          k = strlen(name2)
                          name2 = substr (name2, 1, k-4)

                          del (name2) # Delete original image frame
                          rename.mode = "hl"
                          rename.field = "all"
                          rename.mode = "hl"
                          rename (name2//'.out.fits', name2) # Rename .out images to original names
                          mv (name2, (var[m]//'/'))          # Move out shifted images
                          cp ('gain_rdnoise', (var[m]//'/'))
                      }

                      del ('cluster.'//var[m]//'.out')                  
                  }
             }
             else {
                 auxi=no
                 list = ('cluster.'//var[m]//'.out')
                 while (fscan(list, name) != EOF) {
                     del (name//'.fits')
                 }
                 print ('')
                 print (' Input new parameter values:')
                 print ('')
                 print (' Actual boxsize value = '//boxsize//' Input new value:')
                 scan (boxsize)
                 print (' Actual bigbox value = '//bigbox//' Input new value:')
                 scan (bigbox)
                 print (' Actual niterate value = '//niterate//' Input new value:')
                 scan (niterate)
                 print (' Actual tolerance value = '//tolerance//' Input new value:')
                 scan (tolerance)
             }

          }

      m=m+1
      }

del ('clusters')
del ('gain_rdnoise')


      print ('                                                       ')
      print (' ----------------------------------------------------- ')
      print (' Script "align2" finished correctly.                   ')
      print (' Move on to the "casleored2" script                    ')
      print ('') 
      print (' Remember this last script must be executed inside the package:')
      print ('')
      print ('             noao/digiphot/daophot                     ')
      print (' ----------------------------------------------------- ')


end
