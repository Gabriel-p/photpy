################################################################################
#
#        ===========================================================
#                      PROCEDURE TO PERFORM IMALIGN
#        ===========================================================
# 
#     To run properly this procedure must be executed inside the package:
#                          'noao.digiphot.daophot'
#
#                            by Gabriel Perren 2009
#
################################################################################


procedure align ()

struct *list {mode="h"}
struct *listm {mode="h"}
struct *listk {mode="h"}

begin

      struct name, name2, namem, namek
      string var[50], var2[50], frame
      int i, j, k, m, o, p, q, x
      real thresh
      bool auxi, auxi2


# ------------------------------------------------------------------------------------
# Control 0
# ------------------------------------------------------------------------------------
    
    if (! defpac ("daophot")) {
    print ('')
        print (' This script must be loaded inside the package noao/digiphot/daophot')
        bye()
    }
    else { # Do nothing
    }
# ------------------------------------------------------------------------------------
# Control 0
# ------------------------------------------------------------------------------------   



# ------------------------------------------------------------------------------------
# Structure to create a single file for each cluster, containing the name of the frame
# files that correspond to the same cluster
# ------------------------------------------------------------------------------------ 

      files ('*.fit', > 'frames')
      list = "frames"
      j=0
      for (i=1; fscan(list, name) != EOF; i=i+1) {
         j = j+1  # This int stores the amount of cluster files in the folder
         var[i] = name

         k = strlen(name)
         if (substr (name, k-3, k) == ".fit") {
             name2 = substr (name, 1, 4)
             print (' Name: '//name2)
         }
         else {
             print (' Wrong image name')
             bye()
         }

          var2[i] = name2 # This variables (var2[1], var2[2], etc.) hold the names of
      }                   # all the star frames, without the ".fit" extension and the 
                          # filter and time exposition information.

      m=1
      while (m<=j) {
          k=1
          while (k<=j) {
              print (var2[m], >> 'tempm')
              print (var2[k], >> 'tempk')
              listm = "tempm"
              listk = "tempk"
              while (fscan (listm, namem) != EOF)
              while (fscan (listk, namek) != EOF)
              if (namem == namek) {                         # This 'while' structure separates the frames of a single cluster
                  print (var[k], >> 'cluster_'//var2[m])    # into a unique file for each cluster 
              }
              del ('tempm')
              del ('tempk')
              k=k+1
          }
          m=m+1
      }

      files ('cluster_*', > 'clusters')
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
          o=1
          list = ('cluster_'//var[m])
          while (fscan(list, name) != EOF) {
              o=o+1
          }
                                            # This 'while' cleans the otput file from the previous 
          o=sqrt(o)                         # 'while' so that the names of the frames are not repited

          fields.mode = 'hl' 
          fields.files = ('cluster_'//var[m])
          fields.fields = "1-1"
          fields.lines = '1-'//o
          fields >> ('cluster.'//var[m])

          del ('cluster_'//var[m])

          m=m+1
      }

      del ('frames')

# ------------------------------------------------------------------------------------
# Structure to create a single file for each cluster, containing the name of the frame
# files that correspond to the same cluster
# ------------------------------------------------------------------------------------ 



# ------------------------------------------------------------------------------------
# Daofind on a cluster's frame in the cluster's file
# ------------------------------------------------------------------------------------ 

      list = "clusters"
      j=0
      for (i=1; fscan(list, name) != EOF; i=i+1) {
          j=j+1
          k = strlen(name)                  # This 'for' stores the unique names of each cluster
          name2 = substr (name, k-3, k)     # in the var[i] variables (this name is a FOUR lettered word)
          var[i] = name2

      m=1
      while (m<=j) {

          print ('')
          print (var[m]//' Frames')
          print ('')
          list = ('cluster.'//var[m])
          while (fscan(list, name) != EOF) {
              print (name)
          }

          print ('')
          print (' Input the name of the frame you want to display')
          print (' and run Daofind on (whitout the .fit extension)')
          print (' This will be the "reference" image')
          print ('')

          auxi=no
          while (auxi==no) {   # This 'while' checks that the input image name is right.
              scan (frame)
              print (frame, > 'temp')
              list = "temp"
              while (fscan (list, name) != EOF)
              del ('temp')

              q=0
              x=0
              list = ('cluster.'//var[m])
              while (fscan(list, name2) != EOF) {
                  x=x+1
                  k = strlen(name2)
                  name2 = substr (name2, 1, k-7)

                  if (name==name2) {
                      q=q+0
                  }
                  else  {
                      q=q+1
                  } 
              }

              if (q==x) {
                  print ('')
                  print (' Image name does not match any name stored in the cluster.'//var[m]//' file')
                  print (' Input image name again') 
                  print ('')
              }
              else {
                  auxi=yes
              }
          }

          print (frame, > 'frame.'//var[m]) # File used by the next script: 'align2'
      }

          print ('')
          print (' Input threshold value (a high threshold finds only the brighter stars)')
          scan (thresh)
          findpars.threshold = thresh

          auxi = no
          while (auxi == no) {
          
              daofind.verif = no
              daofind.verb = no
              daofind.interactive = no
              daofind ((frame), "default")

              display ((frame), 1)
    
              tvmark.interactive = no
              tvmark.mark = 'point'
              tvmark.font = "raster"
              tvmark.color = 204
              tvmark.number = no
              tvmark.label = no
              tvmark (1, (frame // '.coo.1'))
      
              # The first tvmark marks the stars in the image, found by 'daofind'
              # The second one (below this), performs the 'interactive' tvmark.
              # If I try to make the first tvmark 'interactive', then
              # it doesn't mark the stars, this way the found star are marked

              tvmark.interactive = yes 
              tvmark.mark = 'point'
              tvmark.font = "raster"
              tvmark.color = 204
              tvmark.number = yes
              tvmark.label = no
              tvmark (1, (frame // '.coo.1')) 


              print ('')
              print ( 'Add or delete stars with the cursor')
              print ('')
              print (' Perform new Daofind search with different "threshold"? (y/n)')
              print ('')
              scan (auxi2)
              if (auxi2==no) {
                  auxi = yes
                  centerpars.calgorithm = "centroid"
                  phot.mode = "hl" 
                  #phot.image = frame//'.fit'
                  #phot.coords = frame//'.coo.1'
                  #phot.output = "default"
                  phot.skyfile = ""
                  phot.centerpars = ""
                  phot.interactive = no
                  phot.radplots = no
                  phot.icommands = ""
                  phot.gcommands = ""
                  phot (frame//'.fit', frame//'.coo.1', "default")
                  del (frame//'.coo.1')
                  txdump.textfiles = frame//'.fit.mag.1'
                  txdump.fields = "XCENTER,YCENTER"
                  txdump.expr = 'MAG[1]!=INDEF'
                  txdump.headers = yes
                  txdump.parameters = yes
                  txdump.mode = "hl"
                  txdump > frame//'.coo.1'
                  del (frame//'.fit.mag.1')
              }
              else {
                  auxi = no
                  del (frame// '.coo.1')
                  print ('')
                  print (' Actual threshold = ' //thresh)
                  print (' Input new threshold value')
                  print ('')
                  print (' (A higher threshold means LESS stars)')
                  print ('')
                  scan (thresh)
                  findpars.threshold = thresh
              }
      }


      m=m+1
}



# ------------------------------------------------------------------------------------
# Daofind on a cluster's frame in the cluster's file
# ------------------------------------------------------------------------------------ 

      print ('                                                       ')
      print (' ----------------------------------------------------- ')
      print (' Script "align" finished correctly.               ')
      print (' Move on to the "align2" script                   ')
      print ('') 
      print (' Remember this last script must be executed inside the package:')
      print ('')
      print ('             images/immatch                     ')
      print (' ----------------------------------------------------- ')


end
