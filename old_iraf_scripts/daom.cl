################################################################################
#
#        ===========================================================
#                        PROCEDURE TO PERFORM DAOMASTER 
#        ===========================================================
# 
#     To run properly this procedure must be executed inside the package:
#                           'noao.digiphot.photcal''
#
#                            by Gabriel Perren 2009
#
################################################################################

procedure daomaster ()

struct *flist {mode="h"}
struct *list {mode="h"}
struct *list2 {mode="h"}

begin

      struct line
      real line2
      string infile, outfile, name, filter, name2, air
      bool auxi, check
      int u, b, v, i, o, k, modes
      real airmass, exptime, max
      real vairmass, uairmass, bairmass, iairmass
      real uvarexp[50], uvarexpmax, bvarexp[50], bvarexpmax, vvarexp[50], vvarexpmax, ivarexp[50], ivarexpmax
      string uvar[50], uvarname, bvar[50], bvarname, vvar[50], vvarname, ivar[50], ivarname

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

    auxi = yes
    while (auxi == yes) {
		    print ('\n Choose execution mode for script:')
		    print ('         [1]: Automatic mode')
		    print ('         [2]: Semi-Automatic mode')
		    print ('         [3]: Manual mode')  
		    print ('\n -Automatic mode starts with a critical match-up radius')
		    print (' of 30 and drops with a step of 1 until 5, which is repeated')
		    print (' 20 times before exit.')
		    print ('\n -Semi-Automatic mode let\'s you input the critical match-up')
		    print (' radius. The output files are chosen automatically.')
		    print ('\n -Manual mode let\'s you choose ALL \'daomaster\' parameters')
		    print (' (including output files)')
		    scan (modes)
		    if ((modes != 1) && (modes != 2) && (modes != 3)) {
		        print ('\n Wrong input. Try again')
		        auxi = yes
		     }
		     else {
		         auxi = no
		     }   
     }

# ------------------------------------------------------------------------------------
# Search for 'list' file
# ------------------------------------------------------------------------------------

      auxi=no
      while (auxi==no) { 

          o=0

          # Search for "list"
          # ------------------------------------------------------------------------------------
          files ('*list', > 'tempfile')
          flist = 'tempfile'
          while (fscan (flist,line) != EOF) {
              o = o + 1
          }
          del ("tempfile")
          # ------------------------------------------------------------------------------------

          if (o!=1) {
              print('')
              print (' ****************************')
              print ('   Missing "list" file.')
              print (' ****************************')
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
                  print ('"list" file found.')
                  auxi=yes
              }
              else {
                  print (' Unknown error. Check code.')
                  bye()
              }
          }
      }

# ------------------------------------------------------------------------------------
# End of 'Search for 'list' file'
# ------------------------------------------------------------------------------------



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
              if (modes != 1) {
		              auxi=no
		              print (' Create default "noche.cfg" file? (y/n)')
		              scan (check)
              }
              else {
                  check = yes
                  auxi = yes
              }    
              if (check) {
                  auxi=yes

                  print ('# Declare the catalog variables', >> 'noche.cfg')
                  print ('', >> 'noche.cfg')
                  print ('catalog', >> 'noche.cfg')
                  print ('', >> 'noche.cfg')
                  print ('V                2', >> 'noche.cfg')
                  print ('error(V)         3', >> 'noche.cfg')
                  print ('BV               4', >> 'noche.cfg')
                  print ('error(BV)        5', >> 'noche.cfg')
                  print ('UB               6', >> 'noche.cfg')
                  print ('error(UB)        7', >> 'noche.cfg')
                  print ('VI               8', >> 'noche.cfg')
                  print ('error(VI)        9', >> 'noche.cfg')
                  print ('', >> 'noche.cfg')
                  print ('# Declare the observations file variables', >> 'noche.cfg')
                  print ('', >> 'noche.cfg')
                  print ('observations', >> 'noche.cfg')
                  print ('', >> 'noche.cfg')
                  print ('x             2              # x coordinate', >> 'noche.cfg')
                  print ('y             3              # y coordinate', >> 'noche.cfg')
                  print ('', >> 'noche.cfg')
                  print ('mV            4              # instrumental magnitude in filter V', >> 'noche.cfg')
                  print ('error(mV)     5              # magnitude error in filter V', >> 'noche.cfg')
                  print ('', >> 'noche.cfg')
                  print ('mB            6              # instrumental magnitude in filter B', >> 'noche.cfg')
                  print ('error(mB)     7              # magnitude error in filter B', >> 'noche.cfg')
                  print ('', >> 'noche.cfg')
                  print ('mU            8              # instrumental magnitude in filter U', >> 'noche.cfg')
                  print ('error(mU)     9              # magnitude error in filter U', >> 'noche.cfg')
                  print ('', >> 'noche.cfg')
                  print ('mI            10             # instrumental magnitude in filter I', >> 'noche.cfg')
                  print ('error(mI)     11             # magnitude error en filter I', >> 'noche.cfg')
                  print ('', >> 'noche.cfg')
                  print ('XV            14             # airmass in filter V', >> 'noche.cfg')
                  print ('XB            15             # airmass in filter B', >> 'noche.cfg')
                  print ('XU            16             # airmass in filter U', >> 'noche.cfg')
                  print ('XI            17             # airmass in filter I', >> 'noche.cfg')
                  print ('', >> 'noche.cfg')
                  print ('', >> 'noche.cfg')
                  print ('# Sample transformation section for the Landolt UBVRI system', >> 'noche.cfg')
                  print ('', >> 'noche.cfg')
                  print ('transformation', >> 'noche.cfg')
                  print ('', >> 'noche.cfg')
                  print ('fit   u1=0.0, u3=0.000', >> 'noche.cfg')
                  print ('const u4=0.0, u2=0.49', >> 'noche.cfg')
                  print ('UFIT : mU = (UB + BV + V) + u1 + u2 * XU + u3 * UB + u4 * UB * XU', >> 'noche.cfg')
                  print ('', >> 'noche.cfg')
                  print ('fit   b1=0.0, b3=0.000', >> 'noche.cfg')
                  print ('const b4=0.0, b2=0.27', >> 'noche.cfg')
                  print ('BFIT : mB = (BV + V) + b1 + b2 * XB + b3 * BV + b4 * BV * XB', >> 'noche.cfg')
                  print ('', >> 'noche.cfg')
                  print ('fit   v1=0.0, v3=0.000', >> 'noche.cfg')
                  print ('const v4=0.0, v2=0.12', >> 'noche.cfg')
                  print ('VFIT : mV = V + v1 + v2 * XV + v3 * BV + v4 * BV * XV', >> 'noche.cfg')
                  print ('', >> 'noche.cfg')
                  print ('#fit   v21=0.0, v23=0.000', >> 'noche.cfg')
                  print ('#const v24=0.0, v22=0.12', >> 'noche.cfg')
                  print ('#ViFIT : mV = V + v21 + v22 * XV + v23 * VI + v24 * VI * XV', >> 'noche.cfg')
                  print ('', >> 'noche.cfg')                  
                  print ('#fit   r1=0.0, r2=0.08, r3=0.000', >> 'noche.cfg')
                  print ('#const r4=0.0', >> 'noche.cfg')
                  print ('#RFIT : mR = (V - VR)  + r1 + r2 * XR + r3 * VR + r4 * VR * XR', >> 'noche.cfg')
                  print ('', >> 'noche.cfg')
                  print ('fit   i1=0.0, i3=0.000', >> 'noche.cfg')
                  print ('const i4=0.0, i2=0.02', >> 'noche.cfg')
                  print ('IFIT : mI = (V - VI) + i1 + i2 * XI + i3 * VI + i4 * VI * XI', >> 'noche.cfg')
                  if (modes != 1) {
		                  print ('')
		                  print (' Edit \'noche.cfg\' file? (y/n)')
		                  scan (check)
		                  if (check) {
		                      vi ('noche.cfg')
		                  }
                  }
              }
              else {
                  auxi=no
              }
          }
          else {
              if (o==1) {
                  print ('')
                  print ('"noche.cfg" file found.')
                  auxi=yes
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



# ------------------------------------------------------------------------------------
# Search for 'noche.ans' file
# ------------------------------------------------------------------------------------

      auxi=no
      while (auxi==no) {

          o=0

          # Search for "noche.ans"
          # ------------------------------------------------------------------------------------
          files ('*noche.ans', > 'tempfile')
          flist = 'tempfile'
          while (fscan (flist,line) != EOF) {
              o = o + 1
          }
          del ("tempfile")
          # ------------------------------------------------------------------------------------

          if (o!=1) {
              print('')
              print (' ****************************')
              print ('   Missing "noche.ans" file.')
              print (' ****************************')
              auxi=no
              print (' Search again (else exit)? (y/n)')
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
                  print ('"noche.ans" file found.')
                  auxi=yes
              }
              else {
                  print (' Unknown error. Check code.')
                  bye()
              }
          }
      }

# ------------------------------------------------------------------------------------
# End of 'Search for 'noche.ans' file'
# ------------------------------------------------------------------------------------



# ---------------------------------------------------------------------------------------------
# Create .txt file to feed DAOMASTER's First run and identificate longest exposures and filters
# ---------------------------------------------------------------------------------------------

      u = 1
      b = 1
      v = 1
      i = 1
      list='list'
      while (fscan (list, infile) != EOF) {

          k = strlen(infile)
          if (substr (infile, k-5, k-1) == ".als.") {
              name = substr (infile, 1, k-6)
		          print ('')
		          print (' Image name: '//name)
          }
          else {
              if (substr (infile, k-5, k-1) == ".mag.") {
			            name = substr (infile, 1, k-6)
				          print ('')
				          print (' Image name: '//name)              
              }
              else {
		              print ('')
		              print (' Infile : '//infile)
		              print (' Image file extension found: '//substr (infile, k-4, k-1))
		              print (' Error in line 318. Check code.')
		              bye()
		          }    
				  }       
       
#          name=substr(infile, 1, stridx (".als.", infile))
          outfile=name//'.txt'
          
          hselect.mode = "hl"
          hselect.images = name
          hselect.fields = "EXPTIME"
          hselect.expr = yes
          hselect > "tempfile"
          flist = 'tempfile'
          while (fscan (flist,line2) != EOF)
          del ("tempfile")
          exptime = line2  # This variable, "exptime", holds the exposure time (duh!)
          print ('Exptime = '//exptime)

          hselect.mode = "hl"
          hselect.images = name
          hselect.fields = "FILTERS"
          hselect.expr = yes
          hselect > "tempfile"
          flist = 'tempfile'
          while (fscan (flist,line) != EOF)
          del ("tempfile")
          filter = line
          print ('Filter = '//filter)

        #
        # Check 'filter' input (shoud be 'u' or 'b' or 'v' or 'i' for CASLEO images)

          if (filter != "U" && filter != "B" && filter != "V" && filter != "I") {
              filter = "a"
              while (filter != "U" && filter != "B" && filter != "V" && filter != "I") {
                  print ('')
                  print (' The filter ID is incorrect. Please enter "U" or "B" or "V" or "I"')
                  scan (filter)
              }
          }
          else { # Do nothing
          }
        #
        # Check 'filter' input (shoud be 'u' or 'b' or 'v' or 'i' for CASLEO images)


          if (filter=="U") {
              print ("Converting ",infile," into ", outfile)
              txdump(textfile=infile,fields="ID,XC,YC,MAG,MERR,CHI,SHARP,CHI,SHARP",expr+,headers-, > outfile)
              uvar[u] = name
              uvarexp[u] = exptime
              u = u+1
          }
          else {
              if (filter=="B") {
                  print ("Converting ",infile," into ", outfile)
                  txdump(textfile=infile,fields="ID,XC,YC,MAG,MERR,CHI,SHARP,CHI,SHARP",expr+,headers-, > outfile)              
                  bvar[b] = name
                  bvarexp[b] = exptime
                  b = b+1
              }
              else {
                  if (filter=="V") {
                      print ("Converting ",infile," into ", outfile)
                      txdump(textfile=infile,fields="ID,XC,YC,MAG,MERR,CHI,SHARP,CHI,SHARP",expr+,headers-, > outfile)
                      vvar[v] = name 
                      vvarexp[v] = exptime    
                      v = v+1                 
                  }
                  else {
                      if (filter=="I") {
                          print ("Converting ",infile," into ", outfile)
                          txdump(textfile=infile,fields="ID,XC,YC,MAG,MERR,CHI,SHARP,CHI,SHARP",expr+,headers-, > outfile)
                          ivar[i] = name
                          ivarexp[i] = exptime
                          i = i+1
                      }
                  }
              }
          }

      }
      
      #
      # This sctructure identifies the frame with the longest exposure in each filter
      print ('')
      print ('')
      max = 30
      
      for (i=1; i<=max; i=i+1) { # 30 is a big enough number to contain all the variables
          if (uvar[2] == ""){
              uvarname = uvar[1] # If there's only one frame then this is the longest exposure
          }
          else {          
              if (uvar[i] == "" || uvar[i+1] == "") {
              }
              else {
                  if (uvarexp[i] >= uvarexp[i+1]) {
                      uvarexpmax = uvarexp[i]
                      uvarexp[i+1] = uvarexp[i]
                      uvarname = uvar[i]
                      uvar[i+1] = uvar[i]
                  }
                  else {
                      uvarexpmax = uvarexp[i+1]
                      uvarname = uvar[i+1]
                  }
              }
          }
      }
          
      for (i=1; i<=max; i=i+1) { # 30 is a big enough number to contain all the variables
          if (bvar[2] == ""){
              bvarname = bvar[1] # If there's only one frame then this is the longest exposure
          }
          else {          
              if (bvar[i] == "" || bvar[i+1] == "") {
              }
              else {
                  if (bvarexp[i] >= bvarexp[i+1]) {
                      bvarexpmax = bvarexp[i]
                      bvarexp[i+1] = bvarexp[i]
                      bvarname = bvar[i]
                      bvar[i+1] = bvar[i]
                  }
                  else {
                      bvarexpmax = bvarexp[i+1]
                      bvarname = bvar[i+1]
                  }
              }
          }
      }
      
      for (i=1; i<=max; i=i+1) { # 30 is a big enough number to contain all the variables
          if (vvar[2] == ""){
              vvarname = vvar[1] # If there's only one frame then this is the longest exposure
          }
          else {          
              if (vvar[i] == "" || vvar[i+1] == "") {
              }
              else {
                  if (vvarexp[i] >= vvarexp[i+1]) {
                      vvarexpmax = vvarexp[i]
                      vvarexp[i+1] = vvarexp[i]
                      vvarname = vvar[i]
                      vvar[i+1] = vvar[i]
                  }
                  else {
                      vvarexpmax = vvarexp[i+1]
                      vvarname = vvar[i+1]
                  }
              }
          }
      }
      
      for (i=1; i<=max; i=i+1) { # 30 is a big enough number to contain all the variables
          if (ivar[2] == ""){
              ivarname = ivar[1] # If there's only one frame then this is the longest exposure
          }
          else {          
              if (ivar[i] == "" || ivar[i+1] == "") {
              }
              else {
                  if (ivarexp[i] >= ivarexp[i+1]) {
                      ivarexpmax = ivarexp[i]
                      ivarexp[i+1] = ivarexp[i]
                      ivarname = ivar[i]
                      ivar[i+1] = ivar[i]
                  }
                  else {
                      ivarexpmax = ivarexp[i+1]
                      ivarname = ivar[i+1]
                  }
              }
          }
      }
      
      #
      # This structure identifies the frame with the longest exposure in each filter
      
      #
      # This structure puts the files with the longest exposure times first in the .mch files
      
      print (' \''//uvarname//'.txt                   \'      0.00      0.00   1.00000   0.00000   0.00000   1.00000', >> 'ufilter.mch')      
      print (' \''//bvarname//'.txt                   \'      0.00      0.00   1.00000   0.00000   0.00000   1.00000', >> 'bfilter.mch')      
      print (' \''//vvarname//'.txt                   \'      0.00      0.00   1.00000   0.00000   0.00000   1.00000', >> 'vfilter.mch')      
      print (' \''//ivarname//'.txt                   \'      0.00      0.00   1.00000   0.00000   0.00000   1.00000', >> 'ifilter.mch')
      
      i = 1
      list='list'
      while (fscan (list, infile) != EOF) {
      
          k = strlen(infile)
          if (substr (infile, k-5, k-1) == ".als.") {
              name = substr (infile, 1, k-6)
		          print ('')
		          print (' Image name: '//name)
          }
          else {
              if (substr (infile, k-5, k-1) == ".mag.") {
			            name = substr (infile, 1, k-6)
				          print ('')
				          print (' Image name: '//name)              
              }
              else {
		              print ('')
		              print (' Infile : '//infile)
		              print (' Image file extension found: '//substr (infile, k-4, k-1))
		              print (' Error in line 527. Check code.')
		              bye()
		          }    
				  }
				        
#          name=substr(infile, 1, stridx (".", infile) - 1)

          hselect.mode = "hl"
          hselect.images = name
          hselect.fields = "FILTERS"
          hselect.expr = yes
          hselect > "tempfile"
          flist = 'tempfile'
          while (fscan (flist,line) != EOF)
          del ("tempfile")
          filter = line

          if (filter=="U") {
              if (name == uvarname) { # Do nothing (ie: do not print the same name again)
              }
              else {
                  print (' \''//name//'.txt                   \'      0.00      0.00   1.00000   0.00000   0.00000   1.00000', >> 'ufilter.mch')
              }
          }
          else {
              if (filter=="B") {
                  if (name == bvarname) { # Do nothing (ie: do not print the same name again)
                  }
                  else {
                      print (' \''//name//'.txt                   \'      0.00      0.00   1.00000   0.00000   0.00000   1.00000', >> 'bfilter.mch')
                  }
              }    
              else {
                  if (filter=="V") {
                      if (name == vvarname) { # Do nothing (ie: do not print the same name again)
                      }
                      else {
                          print (' \''//name//'.txt                   \'      0.00      0.00   1.00000   0.00000   0.00000   1.00000', >> 'vfilter.mch')
                      }
                  } 
                  else {
                      if (filter=="I") {
                          if (name == ivarname) { # Do nothing (ie: do not print the same name again)
                          }
                          else {
                              print (' \''//name//'.txt                   \'      0.00      0.00   1.00000   0.00000   0.00000   1.00000', >> 'ifilter.mch')
                          }
                      } 
                  }
              }
          }

      i = i+1

      }
      #
      # End of 'This sctructure puts the files with the longest exposure times first in the .mch files'
          
      print ('')
      print ('')    
      print ('Longest U exposure: '//uvarname)
      print ('Longest B exposure: '//bvarname)
      print ('Longest V exposure: '//vvarname)
      print ('Longest I exposure: '//ivarname)
      print ('')

      if (modes != 1) {
		      print ('')
		      print (' Were the longest exposure time files correctly identified')
		      print (' (else manually edit .mch files)? (y/n)')
		      scan (check)
		      if (check) {
		      }
		      else {
		          print ('')
		          print (' Continue (have files been manually updated so that')
		          print (' the longest exposure is first)? Else, "no" means "exit" (y/n)')
		          scan (check)
		          if (check) {
		          }
		          else {
		              bye()
		          }    
		      }
      }
      else {
		      print ('Longest U exposure: '//uvarname, >> 'exptime.found')
		      print ('Longest B exposure: '//bvarname, >> 'exptime.found')
		      print ('Longest V exposure: '//vvarname, >> 'exptime.found')
		      print ('Longest I exposure: '//ivarname, >> 'exptime.found')
      }

# ---------------------------------------------------------------------------------------------
# End of 'Create .txt file to feed DAOMASTER's First run and identificate longest exposures and filters'
# ---------------------------------------------------------------------------------------------



# ---------------------------------------------------------------------------------------------
# Creating 'daom.mch' file'
# ---------------------------------------------------------------------------------------------

      o = 0
      # Search for "daom.mch"
      # ------------------------------------------------------------------------------------
      files ('*daom.mch', > 'tempfile')
      flist = 'tempfile'
      while (fscan (flist,line) != EOF) {
          o = o + 1
      }
      del ("tempfile")
      # ------------------------------------------------------------------------------------

      if (o!=1) {
      }
      else {
          if (o==1) {
              print ('')
              print (' "daom.mch" file found. Deleting')
              delete ('*daom.mch')
          }
          else {
          }
      }

      print ('')
      print (' Creating \'daom.mch\' file with the following order')
      print ('')
      print (' vfilter, bfilter, ufilter, ifilter')
      print ('')
      print (' \'vfilter.mag                   \'      0.00      0.00   1.00000   0.00000   0.00000   1.00000', >> 'daom.mch')
      print (' \'bfilter.mag                   \'      0.00      0.00   1.00000   0.00000   0.00000   1.00000', >> 'daom.mch')
      print (' \'ufilter.mag                   \'      0.00      0.00   1.00000   0.00000   0.00000   1.00000', >> 'daom.mch')
      print (' \'ifilter.mag                   \'      0.00      0.00   1.00000   0.00000   0.00000   1.00000', >> 'daom.mch')

# ---------------------------------------------------------------------------------------------
# End of 'Creating 'daom.mch' file'
# ---------------------------------------------------------------------------------------------



# ---------------------------------------------------------------------------------------------
# "DAOMASTER" run
# ---------------------------------------------------------------------------------------------

      auxi = no
      while (auxi==no) {
          print ('')
          print (' DAOMASTER RUN')
          
		      if (modes == 1) {
		          daomaster_auto
		      }
		      if (modes == 2) {
		          daomaster_semi
		      }  
		      if (modes == 3) {
		          daomaster_manual
		      }              

          if (modes != 1) {
              print ('')
              print (' Continue (else delete all files and do-over)? (y/n)')
              scan (check)
          }
          else {
              check = yes
          }    
          if (check) {
              auxi=yes
          }
          else {
              delete ("*filter.mag")
              delete ("*daom.raw")
          } 
      }    
       
          
# ---------------------------------------------------------------------------------------------
# End of '"DAOMASTER" run'
# ---------------------------------------------------------------------------------------------



# ---------------------------------------------------------------------------------------------
# Search for *filter.mag files, created after DAOMASTER run
# ---------------------------------------------------------------------------------------------

      auxi=no
      while (auxi==no) {

          o=0

          # Search for "ufilter.mag"
          # ------------------------------------------------------------------------------------
          files ('*ufilter.mag', > 'tempfile')
          flist = 'tempfile'
          while (fscan (flist,line) != EOF) {
              o = o + 1
          }
          del ("tempfile")
          # ------------------------------------------------------------------------------------

          # Search for "vfilter.mag"
          # ------------------------------------------------------------------------------------
          files ('*vfilter.mag', > 'tempfile')
          flist = 'tempfile'
          while (fscan (flist,line) != EOF) {
              o = o + 1
          }
          del ("tempfile")
          # ------------------------------------------------------------------------------------

          # Search for "bfilter.mag"
          # ------------------------------------------------------------------------------------
          files ('*bfilter.mag', > 'tempfile')
          flist = 'tempfile'
          while (fscan (flist,line) != EOF) {
              o = o + 1
          }
          del ("tempfile")
          # ------------------------------------------------------------------------------------

          # Search for "ifilter.mag"
          # ------------------------------------------------------------------------------------
          files ('*ifilter.mag', > 'tempfile')
          flist = 'tempfile'
          while (fscan (flist,line) != EOF) {
              o = o + 1
          }
          del ("tempfile")
          # ------------------------------------------------------------------------------------

          if (o!=4) {
              print('')
              print (' ****************************')
              print ('         Missing file(s).    ')
              print (' ****************************')
              print ('')
              print (' One or more \'*filter.mag\' files were not found.')
              print (' Search again (else exit)? (y/n)')
              scan (check)
              if (check) {
              }
              else {
                  bye()
              }   
              auxi=no
          }
          else {
              if (o==4) {
                  print ('')
                  print ('All "*filter.mag" files found.')
                  auxi=yes
              }
              else {
                  print (' Unknown error. Check code.')
                  bye()
              }
          }
      }

# ---------------------------------------------------------------------------------------------
# End of 'Search for *filter.mag files, created after DAOMASTER run'
# ---------------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------------
# Search for 'daom.raw' file, created after DAOMASTER run
# ---------------------------------------------------------------------------------------------

      auxi=no
      while (auxi==no) {

          o=0

          # Search for "daom.raw"
          # ------------------------------------------------------------------------------------
          files ('*daom.raw', > 'tempfile')
          flist = 'tempfile'
          while (fscan (flist,line) != EOF) {
              o = o + 1
          }
          del ("tempfile")
          # ------------------------------------------------------------------------------------

          if (o!=1) {
              print('')
              print (' ****************************')
              print ('   Missing "daom.raw" file.  ')
              print (' ****************************')
              print ('')
              print (' File \'daom.raw\' was not found.')
              print (' Search again (else exit)? (y/n)')
              scan (check)
              if (check) {
              }
              else {
                  bye()
              }              
              auxi=no
          }
          else {
              if (o==1) {
                  print ('')
                  print ('"daom.raw" file found.')
                  auxi=yes
              }
              else {
                  print (' Unknown error. Check code.')
                  bye()
              }
          }
      }

# ---------------------------------------------------------------------------------------------
# End of 'Search for 'daom.raw' file, created after DAOMASTER run'
# ---------------------------------------------------------------------------------------------



# ---------------------------------------------------------------------------------------------
# Extract AIRMASS from longest exposures for each filter
# ---------------------------------------------------------------------------------------------

      air = ""

          hselect.images = vvarname//'.fit'
          hselect.fields = "AIRMASS"
          hselect.expr = yes
          hselect.mode = "hl"
          hselect > 'vair'
          list2 = 'vair'
          while (fscan(list2, vairmass) != EOF)
          air = air // "  " // str(vairmass)

          hselect.images = bvarname//'.fit'
          hselect.fields = "AIRMASS"
          hselect.expr = yes
          hselect.mode = "hl"
          hselect > 'bair'
          list2 = 'bair'
          while (fscan(list2, bairmass) != EOF)
          air = air // "  " // str(bairmass)

          hselect.images = uvarname//'.fit'
          hselect.fields = "AIRMASS"
          hselect.expr = yes
          hselect.mode = "hl"
          hselect > 'uair'
          list2 = 'uair'
          while (fscan(list2, uairmass) != EOF)
          air = air // "  " // str(uairmass)

          hselect.images = ivarname//'.fit'
          hselect.fields = "AIRMASS"
          hselect.expr = yes
          hselect.mode = "hl"
          hselect > 'iair'
          list2 = 'iair'
          while (fscan(list2, iairmass) != EOF)
          air = air // "  " // str(iairmass)

      infile = 'daom.raw'
      #n = nfiltros  --> = 4 for CASLEO and CTIO
      flist = infile
      while (fscan(flist, line) != EOF) {
          line = line // air
          print (line, >> (substr(infile, 1, stridx (".", infile) - 1) // ".obs"))
      }
      print ('')
      print ('Output file = '// (substr(infile, 1, stridx (".", infile) - 1) // ".obs"))

      del ('vair')
      del ('uair')
      del ('bair')
      del ('iair')
      
# ---------------------------------------------------------------------------------------------
# End of 'Extract AIRMASS from longest exposures for each filter'
# ---------------------------------------------------------------------------------------------

      del ('*.mch')
      if (modes != 1) {
		      print ('')
		      print (' Delete *.txt files? (y/n)')
		      scan (check)
      }
      else {
          check = yes
      }    
      if (check) {
          del ('*.txt')
      }

# ------------------------------------------------------------------------
# 1º Invertfit task
# ------------------------------------------------------------------------   

       print ('                                                       ')
       print (' ----------------------------------------------------- ')
       print (' 1st Invertfit task (all equations)                    ')
       print (' ----------------------------------------------------- ')
       invertfit.mode = "hl"
       invertfit.observations = "daom.obs"
       invertfit.config = "noche.cfg"
       invertfit.parameters = "noche.ans"
       invertfit.calib = "night_obs.out.BV"  #  <---------------------------- FINAL OUTPUT FILE!!
       invertfit.print = "x,y"
       invertfit

       rename.field = 'all'
       rename ('noche.cfg', 'noche.cfg.BV')

       print('#N ID       XCENTER    YCENTER  VS      ERVS    BV      ERBV    UB      ERUB    VI      ERVI', >> 'auxiliar')
       print('#U ##       pixels     pixels   mag     mag     mag     mag     mag     mag     mag     mag', >> 'auxiliar')
       print('#F %-11d    %-8.3f     %-8.3f   %-8.3f  %-8.3f  %-8.3f  %-8.3f  %-8.3f  %-8.3f  %-8.3f  %-8.3f', >> 'auxiliar')
       
       concatenate.input_files = ('auxiliar,night_obs.out.BV') 
       concatenate.output_file = ('concatenated')
       concatenate.out_type = 'text'
       concatenate.append = yes
       concatenate.mode = 'hl'
       concatenate > 'concatenated'
       
       del ('auxiliar')
      
       txdump.mode = 'hl' 
       txdump.textfile = ('concatenated')
       txdump.headers = yes
       txdump.fields = 'ID,XCENTER,YCENTER,VS,ERVS,BV,ERBV,UB,ERUB,VI,ERVI'
       txdump.expr = 'VS!=INDEF'
       txdump > auxiliar.1
       
       del ('concatenated')
      
       txdump.mode = 'hl' 
       txdump.textfile = ('auxiliar.1')
       txdump.headers = yes
       txdump.fields = 'ID,XCENTER,YCENTER,VS,ERVS,BV,ERBV,UB,ERUB,VI,ERVI'
       txdump.expr = 'BV!=INDEF'
       txdump > auxiliar.2
      
       del ('auxiliar.1')
      
       txdump.mode = 'hl' 
       txdump.textfile = ('auxiliar.2')
       txdump.headers = yes
       txdump.fields = 'ID,XCENTER,YCENTER,VS,ERVS,BV,ERBV,UB,ERUB,VI,ERVI'
       txdump.expr = 'UB!=INDEF'
       txdump > auxiliar.3

       del ('auxiliar.2')
      
       print ('numero x y u-b eu-b b-v eb-v v ev v-i ev-i x', > 'shift_entrada.BV')
       txdump.mode = 'hl' 
       txdump.textfile = ('auxiliar.3')
       txdump.headers = no
       txdump.fields = 'ID,XCENTER,YCENTER,UB,ERUB,BV,ERBV,VS,ERVS,VI,ERVI'
       txdump.expr = 'VI!=INDEF'
       txdump >> shift_entrada.BV
      
       del ('auxiliar.3')      
       
# ------------------------------------------------------------------------
# 1º Invertfit task
# ------------------------------------------------------------------------


      print ('# Declare the catalog variables', >> 'noche.cfg')
      print ('', >> 'noche.cfg')
      print ('catalog', >> 'noche.cfg')
      print ('', >> 'noche.cfg')
      print ('V                2', >> 'noche.cfg')
      print ('error(V)         3', >> 'noche.cfg')
      print ('BV               4', >> 'noche.cfg')
      print ('error(BV)        5', >> 'noche.cfg')
      print ('UB               6', >> 'noche.cfg')
      print ('error(UB)        7', >> 'noche.cfg')
      print ('VI               8', >> 'noche.cfg')
      print ('error(VI)        9', >> 'noche.cfg')
      print ('', >> 'noche.cfg')
      print ('# Declare the observations file variables', >> 'noche.cfg')
      print ('', >> 'noche.cfg')
      print ('observations', >> 'noche.cfg')
      print ('', >> 'noche.cfg')
      print ('x             2              # x coordinate', >> 'noche.cfg')
      print ('y             3              # y coordinate', >> 'noche.cfg')
      print ('', >> 'noche.cfg')
      print ('mV            4              # instrumental magnitude in filter V', >> 'noche.cfg')
      print ('error(mV)     5              # magnitude error in filter V', >> 'noche.cfg')
      print ('', >> 'noche.cfg')
      print ('mB            6              # instrumental magnitude in filter B', >> 'noche.cfg')
      print ('error(mB)     7              # magnitude error in filter B', >> 'noche.cfg')
      print ('', >> 'noche.cfg')
      print ('mU            8              # instrumental magnitude in filter U', >> 'noche.cfg')
      print ('error(mU)     9              # magnitude error in filter U', >> 'noche.cfg')
      print ('', >> 'noche.cfg')
      print ('mI            10             # instrumental magnitude in filter I', >> 'noche.cfg')
      print ('error(mI)     11             # magnitude error en filter I', >> 'noche.cfg')
      print ('', >> 'noche.cfg')
      print ('XV            14             # airmass in filter V', >> 'noche.cfg')
      print ('XB            15             # airmass in filter B', >> 'noche.cfg')
      print ('XU            16             # airmass in filter U', >> 'noche.cfg')
      print ('XI            17             # airmass in filter I', >> 'noche.cfg')
      print ('', >> 'noche.cfg')
      print ('', >> 'noche.cfg')
      print ('# Sample transformation section for the Landolt UBVRI system', >> 'noche.cfg')
      print ('', >> 'noche.cfg')
      print ('transformation', >> 'noche.cfg')
      print ('', >> 'noche.cfg')
      print ('#fit   u1=0.0, u3=0.000', >> 'noche.cfg')
      print ('#const u4=0.0, u2=0.49', >> 'noche.cfg')
      print ('#UFIT : mU = (UB + BV + V) + u1 + u2 * XU + u3 * UB + u4 * UB * XU', >> 'noche.cfg')
      print ('', >> 'noche.cfg')
      print ('#fit   b1=0.0, b3=0.000', >> 'noche.cfg')
      print ('#const b4=0.0, b2=0.27', >> 'noche.cfg')
      print ('#BFIT : mB = (BV + V) + b1 + b2 * XB + b3 * BV + b4 * BV * XB', >> 'noche.cfg')
      print ('', >> 'noche.cfg')
      print ('#fit   v1=0.0, v3=0.000', >> 'noche.cfg')
      print ('#const v4=0.0, v2=0.12', >> 'noche.cfg')
      print ('#VFIT : mV = V + v1 + v2 * XV + v3 * BV + v4 * BV * XV', >> 'noche.cfg')
      print ('', >> 'noche.cfg')
      print ('fit   v21=0.0, v23=0.000', >> 'noche.cfg')
      print ('const v24=0.0, v22=0.12', >> 'noche.cfg')
      print ('ViFIT : mV = V + v21 + v22 * XV + v23 * VI + v24 * VI * XV', >> 'noche.cfg')
      print ('', >> 'noche.cfg')                  
      print ('#fit   r1=0.0, r2=0.08, r3=0.000', >> 'noche.cfg')
      print ('#const r4=0.0', >> 'noche.cfg')
      print ('#RFIT : mR = (V - VR)  + r1 + r2 * XR + r3 * VR + r4 * VR * XR', >> 'noche.cfg')
      print ('', >> 'noche.cfg')
      print ('fit   i1=0.0, i3=0.000', >> 'noche.cfg')
      print ('const i4=0.0, i2=0.02', >> 'noche.cfg')
      print ('IFIT : mI = (V - VI) + i1 + i2 * XI + i3 * VI + i4 * VI * XI', >> 'noche.cfg')

# ------------------------------------------------------------------------
# 2º Invertfit task
# ------------------------------------------------------------------------   

       print ('                                                       ')
       print (' ----------------------------------------------------- ')
       print (' 2nd Invertfit task (equation for V with (V-I)         ')
       print (' ----------------------------------------------------- ')
       invertfit.mode = "hl"
       invertfit.observations = "daom.obs"
       invertfit.config = "noche.cfg"
       invertfit.parameters = "noche.ans"
       invertfit.calib = "night_obs.out.VI"  #  <---------------------------- FINAL OUTPUT FILE!!
       invertfit.print = "x,y"
       invertfit

      rename.field = 'all'
      rename ('noche.cfg', 'noche.cfg.VI')
      
# ------------------------------------------------------------------------
# 2º Invertfit task
# ------------------------------------------------------------------------

       print ('')
       print (' FINAL output file\'s names are: \'night_obs.out.BV\' and \'night_obs.out.VI\' ')
       print ('')
       print (' FINAL output file for \'shift\' is: \'shift_entrada.BV\' ')
       print ('')


end

