################################################################################
#
#        ===========================================================
#               SCRIPT TO PERFORM ZEROCOMBINE, FLATCOMBINE, ETC..
#        ===========================================================
# 
#     To run properly this procedure must be executed inside the package:
#                           'noao.imred.ccdred'
#
#                            by Gabriel Perren 2009
#
################################################################################  



bool check, check2
int j,k,q,h,i,n,m,e,o
struct *flist, *flist2, *flist3
struct line,line2,line3,line4
string file_name,file_name2,file_name3,BIASSEC,TRIMSEC,biassec2,trimsec2
real gain, rdnoise


# ------------------------------------------------------------------------------------
# Control 0
# ------------------------------------------------------------------------------------
      print ('')
      print (' Before running this script, the task "setairmass" must be used')
      print (' inside the folder where ALL the frames (bias, flats and standars)')
      print (' are, in the following way:')
      print ('')
      print (' setarimass *.fit')
      print ('')
      print (' This task is located inside: noao/astutil')
      print ('')
      print ('Is this correct (else exit)? (y/n)')
      print ('')
      scan (check)
      if (check == yes) {
      }
      else {
          bye()
      }

# ------------------------------------------------------------------------------------
# Control 0
# ------------------------------------------------------------------------------------



# ------------------------------------------------------------------------------------
# Control 1
# ------------------------------------------------------------------------------------
      if (! defpac ("ccdred")) {
          print ('')
          print (' This script must be loaded inside the package noao/imred/ccdred')
          bye()
# ------------------------------------------------------------------------------------
# Control 1
# ------------------------------------------------------------------------------------



# ------------------------------------------------------------------------------------
# Activate following controls
# ------------------------------------------------------------------------------------

      print ('-----------------------------------------------------------------')
      print ('Activate remaining 3 controls? (y/n)')
      print ('If you dont know what this is, you should answer "yes"')
      print ('')
      scan (check)

      if (check == yes) {

          # ------------------------------------------------------------------------------------
          # 1º Control
          # ------------------------------------------------------------------------------------
          print ('-----------------------------------------------------------------')
          print ('1st Control')
          print ('TODOS LOS ARCHIVOS (FLATS, BIAS, FRAMES DE CUMULOS Y ESTANDARS)')
          print ('DEBEN ENCONTRARSE EN LA CARPETA DE TRABAJO')
          print ('')
          print ('Se encuentran todos los archivos en la carpeta de trabajo? (y/n)')
          print ('')
          scan (check2)
          if (check2 == yes) { # Do nothing
          }
          else {
              bye()
          }
          # ------------------------------------------------------------------------------------
          # 1º Control
          # ------------------------------------------------------------------------------------


          # ------------------------------------------------------------------------------------
          # 2º Control
          # ------------------------------------------------------------------------------------
          print ('-----------------------------------------------------------------')
          print ('2nd Control')
          print ('Esta tarea asume que los nombres de las imagenes de FLATS, BIAS,')
          print ('FRAMES DE CUMULOS y ESTANDARS se encuentran en el siguiente formato:')
          print (' ')
          print ('- para FLATS: "flatwx.fit" o "skywx.fit"; donde "w" representa la')
          print ('letra del filtro correspondiente (U, B, V, R, I) y "x" representa')
          print ('el numero de flat. De esta manera el tercer FLAT tomado en el filtro B')
          print ('tendria como nombre: "flatb3.fit" o "skyb3.fit".')
          print (' ')
          print ('- para BIAS: "biasx.fit"; donde "x" representa el numero de bias.')
          print ('De esta manera el bias numero 5 tendria como nombre: "bias5.fit".')
          print (' ')
          print ('- para FRAMES DE CUMULOS Y ESTANDARS: "abcdwx.fit"; donde "abcd"')
          print ('representa las cuatro letras que identifican al cumulo correspondiente')
          print ('(es indistinto las cuatro letras que sean), "w" representa la letra')
          print ('que identifica, al filtro y "x" representa el tiempo de exposicion.')
          print ('De esta manera un frame de 30 seg tomado en el filtro B para el cumulo')
          print ('Bochum 13 tendria como nombre: "bo13b30.fit".')
          print (' ')
          print ('Son correctos los nombres de los archivos? (y/n)')
          print ('')
          scan (check2)
          if (check2 == yes) { # Do nothnig
          }
          else {
              bye()
          }
          # ------------------------------------------------------------------------------------
          # 2º Control
          # ------------------------------------------------------------------------------------
    
    
          # ------------------------------------------------------------------------------------
          # 3º Control
          # ------------------------------------------------------------------------------------
          print ('-----------------------------------------------------------------')
          print ('3rd Control ')
          print ('Es necesario conocer los valores de TRIMSEC y BIASSEC')
          print ('(si es que son necesarios) en los pasos siguientes.')
          print (' ')
          print ('Conoce estos valores? (y/n)')
          print ('')
          scan (check2)
          if (check2 == yes) { # Do nothnig
          }
          else {
              bye()
          }
          # ------------------------------------------------------------------------------------
          # 3º Control
          # ------------------------------------------------------------------------------------
    
    

          # ------------------------------------------------------------------------------------
          # 4º Control
          # ------------------------------------------------------------------------------------
          #print ('-----------------------------------------------------------------')
          #print ('4th Control')
          #print ('IMPORTANTE')
          #print ('ANTES de correr este script debe haber ingresado en la')
          #print ('linea de comandos de IRAF (o debe estar configurado por')
          #print ('defecto) la sentencia:')
          #print ('')
          #print ('setinstrument=direct')
          #print (' ')
          #print ('Ha ingresado el comando o se encuentra configurado por defecto? (y/n)')
          #print ('')
          #scan (check2)
          #if (check2 == yes) { # Do nothnig
          #}
          #else {
          #    bye()
          #}

          setinstrument.instrument = "direct"

          # ------------------------------------------------------------------------------------
          # 4º Control
          # ------------------------------------------------------------------------------------
      }
      else { # Do nothnig
      }
      
# ------------------------------------------------------------------------------------
# End of 'Activate following controls'
# ------------------------------------------------------------------------------------



# ------------------------------------------------------------------------------------
# Add TRIMSEC, BIASSES, GAIN and RDNOISE to the image headers of all the files
# ------------------------------------------------------------------------------------

      o = 1
      h = 1
      i = 1

      while (o != 3) {
          o=1

          print ('Write BIASSEC to header?  1 = [25:1340,1301:1310], 2 = Other value, 3 = Do not write')
          print ('(Not needed if no "overscan" correction is to be performed)')
          print ('')
          scan (q)
          
          if (q == 2) {
              print ('')
              print ('Input BIASSEC with the following format: [num1:num2,num3:num4]')
              print ('')
              scan (BIASSEC)
          }
          else { # Do nothing
          }

          print ('')
          print ('Write TRIMSEC to header?')
          print ('(Not needed if no "trim" correction is to be performed)')
          print ('1 = [25:1340,1:1299], 2= [10:1339,1:1295], 3 = Other value, 4 = Do not write')
          print ('')
          scan (e)
          
          if (e == 3) {
              print ('Input TRIMSEC with the following format: [num1:num2,num3:num4]')
              print ('')
              scan (TRIMSEC)
          }
          else { # Do nothing
          }

          
          if (q == 1) {
              o = o+1
          }
          else {
              if (q == 2) {
                  o = o+1
              }
              else {
                  if (q == 3) {
                      o = o+1
                  }
                  else {
                      o = o+0
                      print ('')
                      print ('Wrong number. Try again..')
                      print ('')
                  }
              }
          }    
          
          if (e == 1) {
              o = o+1
          }
          else {
              if (e == 2) {
                  o = o+1
              }
              else {
                  if (e == 3) {
                      o = o+1
                  }
                  else {
                      if (e == 4) {
                          o = o+1
                      }
                      else {
                          o = o+0
                          print ('')
                          print ('Wrong number. Try again.')
                          print ('')
                      }
                  }
              }
          }          

      } # This bracket closes the 'o != 3' while


      print ('Adding TRIMSEC and BIASSEC...')

      if (q == 1) {
          hedit (images="*.fit", fields="BIASSEC", value="[25:1340,1301:1310]", add=yes, addonly=yes, delete=no, verify=no, update=yes)
          biassec2 = "[25:1340,1301:1310]"
      }
      else {
          if (q == 2) {
              hedit (images="*.fit", fields="BIASSEC", value=(BIASSEC), add=yes, addonly=no, delete=no, verify=no, update=yes)
              biassec2 = BIASSEC 
          }
          else {
              if (q == 3) {
                  h=2
                  hedit (images="*.fit", fields="BIASSEC", value="", add=yes, addonly=no, delete=no, verify=no, update=yes)
                  biassec2 = ""
              }
              else {
              }          
          }
      }


      if (e == 1) {
          hedit (images="*.fit", fields="TRIMSEC", value="[25:1340,1:1299]", add=yes, addonly=yes, delete=no, verify=no, update=yes)
          trimsec2 = "[25:1340,1:1299]"
      }
      else {
          if (e == 2) {
              hedit (images="*.fit", fields="TRIMSEC", value="[10:1339,1:1295]", add=yes, addonly=no, delete=no, verify=no, update=yes)
              trimsec2 = "[10:1339,1:1295]"
          }
          else {
              if (e == 3) {
                  hedit (images="*.fit", fields="TRIMSEC", value=(TRIMSEC), add=yes, addonly=no, delete=no, verify=no, update=yes)
                  trimsec2 = TRIMSEC
              }
              else {
                  if (e == 4) {
                    i=2
                    hedit (images="*.fit", fields="TRIMSEC", value="", add=yes, addonly=no, delete=no, verify=no, update=yes)
                    trimsec2 = ""
                  }
                  else {
                  }
              }
          }
      }

      print ('Adding GAIN y RDNOISE...')
      print ('')
      print ('Input GAIN value (CASLEO: 3.98e/ADU)')
      scan (gain)
      print ('Input RDNOISE value (CASLEO: 4.09e)')
      scan (rdnoise)
      hedit (images="*.fit", fields="GAIN", value=(gain), add=no, addonly=yes, delete=no, verify=no, update=yes)
      hedit (images="*.fit", fields="RDNOISE", value=(rdnoise), add=no, addonly=yes, delete=no, verify=no, update=yes)

      print (gain,rdnoise, > "gain_rdnoise") # This file will be used by the 'casleopsf' script later on
      
# ------------------------------------------------------------------------------------
# End of 'Add TRIMSEC, BIASSES, GAIN and RDNOISE to the image headers of all the files'
# ------------------------------------------------------------------------------------



# ------------------------------------------------------------------------------------
# Move BIAS frames and create folder to move FLATS later on
# ------------------------------------------------------------------------------------

      print ('Moving biad frames to folder bias/...')
      print ('')
      mkdir bias
      mkdir flats # In this folder we'll store the FLATS frames later on 
      mv bias*.fit bias/
      
# ------------------------------------------------------------------------------------
# End of 'Move BIAS frames and create folder to move FLATS later on'
# ------------------------------------------------------------------------------------



# ------------------------------------------------------------------------------------
# Rfits on BIAS frames to set datatype="ushort" (pixel value max: 64000)
# ------------------------------------------------------------------------------------

      print ('Rfits on BIAS frames to set datatype="ushort"...')
      files bias/*.fit > fitfiles
      rename (files="bias/*.fit", newname="fits", field="extn")
      files bias/*.fits > fitsfiles

      rfits (fits_file="@fitsfiles", file_list="0", iraf_file="@fitfiles", make_image=yes,
      long_header=no, short_header=yes, datatype="ushort", blank=0., scale=yes, oldirafname=no,
      offset=0)

      delete bias/*.fits
      delete fitfiles
      delete fitsfiles
      
# ------------------------------------------------------------------------------------
# End of 'Rfits on BIAS frames to set datatype="ushort" (pixel value max: 64000)'
# ------------------------------------------------------------------------------------



# ------------------------------------------------------------------------------------
# Add missing data to FLAT image headers
# ------------------------------------------------------------------------------------

      print ('Adding missing data to FLAT image headers')
      print ('')
      
      files sky*.* > temp.sky
      files flat*.* > temp.flats

      file_name2 = "temp.flats"
      file_name3 = "temp.sky"

      flist2 = file_name2
      flist3 = file_name3

      n=0
      while (fscan (flist2,line3) != EOF) {
          n = n + 1
      }

      m=0
      while (fscan (flist3,line4) != EOF) {
          m = m + 1
      }

      if (n == 0 && m != 0) {
          print ('No flat*.* files')
          print ('Using sky*.* files...')
          k=2
          sleep 2
      }
      else {
          if (m == 0 && n != 0) {
              print ('No sky*.* files')
              print ('Using flat*.* files...')
              k=1
              sleep 2
          }
          else {
              if (m == 0  && n == 0) {
              print ('No FLAT files of any kind found.')
              print ('Closing...')
              sleep 2
              bye()
              }
              else { # Do nothing
              }
          }
      } 
   
      delete temp.flats
      delete temp.sky

      if (k == 1) {
          hedit (images="flatu*.fit", fields="FILTERS", value="U", add=yes, addonly=no, delete=no, verify=no, update=yes)
          hedit (images="flatb*.fit", fields="FILTERS", value="B", add=yes, addonly=no, delete=no, verify=no, update=yes)
          hedit (images="flatv*.fit", fields="FILTERS", value="V", add=yes, addonly=no, delete=no, verify=no, update=yes)
          hedit (images="flati*.fit", fields="FILTERS", value="I", add=yes, addonly=no, delete=no, verify=no, update=yes)
          hedit (images="flat*.fit", fields="IMAGETYP", value="flat", add=yes, addonly=no, delete=no, verify=no, update=yes)
          mv flat*.fit flats/
      }
      else {
          if (k == 2){
              hedit (images="skyu*.fit", fields="FILTERS", value="U", add=yes, addonly=yes, delete=no, verify=no, update=yes)
              hedit (images="skyb*.fit", fields="FILTERS", value="B", add=yes, addonly=yes, delete=no, verify=no, update=yes)
              hedit (images="skyv*.fit", fields="FILTERS", value="V", add=yes, addonly=yes, delete=no, verify=no, update=yes)
              hedit (images="skyi*.fit", fields="FILTERS", value="I", add=yes, addonly=yes, delete=no, verify=no, update=yes)
              hedit (images="sky*.fit", fields="IMAGETYP", value="flat", add=yes, addonly=yes, delete=no, verify=no, update=yes)
              mv sky*.fit flats/
          }
          else {
              print ('No FLAT files of any kind found.') # Redundant control.
              print ('Closing...')
              sleep 2
              bye()
          }
      }

# ------------------------------------------------------------------------------------
# End of 'Add missing data to FLAT image headers'
# -----------------------------------------------------------------------------------



# ------------------------------------------------------------------------------------
# Rfits on FLAT frames to set datatype="ushort" (pixel value max: 64000)
# ------------------------------------------------------------------------------------

      print ('Rfits on FLAT frames...')
      print ('')
      files flats/*.fit > fitfiles
      rename (files="flats/*.fit", newname="fits", field="extn")
      files flats/*.fits > fitsfiles

      rfits (fits_file="@fitsfiles", file_list="0", iraf_file="@fitfiles", make_image=yes,
      long_header=no, short_header=yes, datatype="ushort", blank=0., scale=yes, oldirafname=no,
      offset=0)

      delete flats/*.fits
      delete fitfiles
      delete fitsfiles

# ------------------------------------------------------------------------------------
# End of 'Rfits on FLAT frames to set datatype="ushort" (pixel value max: 64000)'
# ------------------------------------------------------------------------------------



# ------------------------------------------------------------------------------------
# Add filter data to standards and clusters frame's image headers
# ------------------------------------------------------------------------------------

      print ('')
      print (' Adding filter data to standards and clusters frames image headers...')

      files *.fit > fitfiles

      file_name = "fitfiles"


      flist = file_name

      while (fscan (flist,line) != EOF) {
          j = strlen (line)
          if (substr (line, j-3, j) == ".fit") {
              line2 = substr (line, 5, 5)
              if (line2 == "u") {
                  hedit (images=(line), fields="FILTERS", value="U", add=yes, addonly=no, delete=no, verify=no, update=yes)
                  print (line, >> 'imagenes')
              }
              else {
                  if (line2 == "b") {
                      hedit (images=(line), fields="FILTERS", value="B", add=yes, addonly=no, delete=no, verify=no, update=yes)
                      print (line, >> 'imagenes')
                  }
                  else {
                      if (line2 == "v") {
                          hedit (images=(line), fields="FILTERS", value="V", add=yes, addonly=no, delete=no, verify=no, update=yes)
                          print (line, >> 'imagenes')
                      }
                      else {
                          if (line2 == "i") {
                              hedit (images=(line), fields="FILTERS", value="I", add=yes, addonly=no, delete=no, verify=no, update=yes)
                              print (line, >> 'imagenes')
                          }
                          else {
                              print ('ERROR IN FILENAME OF FILE: '//line)
                              print ('Filter value found: '//line2)
                              print ('Actual extension found: '//(substr (line, j-3, j))
                              print ('The file "imagenes" holds the names of the files correctly modified')
                              print ('Closing...')
                              bye()
                          }
                      }
                  }
              }
          }
          else {
              print ('ERROR IN FILENAME OF FILE: '//line)
              print ('Actual extension found: '//(substr (line, j-3, j))
              print ('The file "imagenes" holds the names of the files correctly modified')
              print ('Closing...')
              bye()
          }
      }

# ------------------------------------------------------------------------------------
# End of 'Add filter data to standards and clusters frame's image headers'
# ------------------------------------------------------------------------------------



# ------------------------------------------------------------------------------------
# Rfits on standards and clusters frames to set datatype="ushort" (pixel value max: 64000)
# ------------------------------------------------------------------------------------

      print ('')
      print ('Rfits on standards and clusters frames...')
      print ('')
      #files *.fit > fitfiles # Already done
      rename (files="*.fit", newname="fits", field="extn")
      files *.fits > fitsfiles

      rfits (fits_file="@fitsfiles", file_list="0", iraf_file="@fitfiles", make_image=yes,
      long_header=no, short_header=yes, datatype="ushort", blank=0., scale=yes, oldirafname=no,
      offset=0)

      delete *.fits
      delete fitfiles
      delete fitsfiles

# ------------------------------------------------------------------------------------
# End of 'Rfits on standards and clusters frames to set datatype="ushort" (pixel value max: 64000)'
# ------------------------------------------------------------------------------------



# ------------------------------------------------------------------------------------
# Zerocombine on BIAS frames
# ------------------------------------------------------------------------------------

      zerocombine ("bias/bias*.fit",
      output="Zero", combine="average", reject="minmax", ccdtype="zero", process=no,
      delete=no, clobber=no, scale="none", statsec="", nlow=0, nhigh=1, nkeep=1,
      mclip=yes, lsigma=3., hsigma=3., rdnoise="RDNOISE", gain="GAIN", snoise="0.",
      pclip=-0.5, blank=0.)

# ------------------------------------------------------------------------------------
# End of 'Zerocombine on BIAS frames'
# ------------------------------------------------------------------------------------



# ------------------------------------------------------------------------------------
# 1º 'ccdproc' run to apply 'biassec', 'trimsec' and 'zerocor'
# ------------------------------------------------------------------------------------

      #biassec (overscan)
      #1 = [25:1340,1301:1310], 2 = Otro valor, 3 = Do not write (h = 2)

      #trimsec (trim)
      #1 = [25:1340,1:1299], 2= [10:1339,1:1295], 3 = Otro valor, 4 = Do not write (i = 2)

      print ('')
      print ('Ccdproc run to apply "overscan", "trim" and "zerocor"')
      print ('')
      print ('Standards and cluster frames...')
      print ('')
      if (h == 2 && i == 2) {
          ccdproc ("*.fit",
          output=" ", ccdtype=" ", max_cache=0, noproc=no, fixpix=no, overscan=no,
          trim=no, zerocor=yes, darkcor=no, flatcor=no, illumcor=no, fringecor=no,
          readcor=no, scancor=no, readaxis="line", fixfile=" ",
          biassec="image", trimsec="image", zero="Zero", dark="",
          flat=" ", illum=" ", fringe="", minreplace=1., scantype="shortscan", nscan=1,
          interactive=no, function="chebyshev", order=4, sample="*", naverage=1,
          niterate=1, low_reject=3., high_reject=3., grow=1.)
      }
      else {
          if (h == 2 && i == 1) {
              ccdproc ("*.fit",
              output=" ", ccdtype=" ", max_cache=0, noproc=no, fixpix=no, overscan=no,
              trim=yes, zerocor=yes, darkcor=no, flatcor=no, illumcor=no, fringecor=no,
              readcor=no, scancor=no, readaxis="line", fixfile=" ",
              biassec="image", trimsec="image", zero="Zero", dark="",
              flat=" ", illum=" ", fringe="", minreplace=1., scantype="shortscan", nscan=1,
              interactive=no, function="chebyshev", order=4, sample="*", naverage=1,
              niterate=1, low_reject=3., high_reject=3., grow=1.)
          }
          else {
              if (h ==1 && i == 2) {
                  ccdproc ("*.fit",
                  output=" ", ccdtype=" ", max_cache=0, noproc=no, fixpix=no, overscan=yes,
                  trim=no, zerocor=yes, darkcor=no, flatcor=no, illumcor=no, fringecor=no,
                  readcor=no, scancor=no, readaxis="line", fixfile=" ",
                  biassec="image", trimsec="image", zero="Zero", dark="",
                  flat=" ", illum=" ", fringe="", minreplace=1., scantype="shortscan", nscan=1,
                  interactive=no, function="chebyshev", order=4, sample="*", naverage=1,
                  niterate=1, low_reject=3., high_reject=3., grow=1.)
              }
              else { #(h == 1 && i == 1)
                  ccdproc ("*.fit",
                  output=" ", ccdtype=" ", max_cache=0, noproc=no, fixpix=no, overscan=yes,
                  trim=yes, zerocor=yes, darkcor=no, flatcor=no, illumcor=no, fringecor=no,
                  readcor=no, scancor=no, readaxis="line", fixfile=" ",
                  biassec="image", trimsec="image", zero="Zero", dark="",
                  flat=" ", illum=" ", fringe="", minreplace=1., scantype="shortscan", nscan=1,
                  interactive=no, function="chebyshev", order=4, sample="*", naverage=1,
                  niterate=1, low_reject=3., high_reject=3., grow=1.)
              }
          }
      }

      print ('')
      print ('Flat frames...')
      print ('')
      if (h == 2 && i == 2) {
          ccdproc ("flats/*.fit",
          output=" ", ccdtype="flat", max_cache=0, noproc=no, fixpix=no, overscan=no,
          trim=no, zerocor=yes, darkcor=no, flatcor=no, illumcor=no, fringecor=no,
          readcor=no, scancor=no, readaxis="line", fixfile=" ",
          biassec="image", trimsec="image", zero="Zero", dark="",
          flat=" ", illum=" ", fringe="", minreplace=1., scantype="shortscan", nscan=1,
          interactive=no, function="chebyshev", order=4, sample="*", naverage=1,
          niterate=1, low_reject=3., high_reject=3., grow=1.)
      }
      else {
          if (h == 2 && i == 1) {
              ccdproc ("flats/*.fit",
              output=" ", ccdtype="flat", max_cache=0, noproc=no, fixpix=no, overscan=no,
              trim=yes, zerocor=yes, darkcor=no, flatcor=no, illumcor=no, fringecor=no,
              readcor=no, scancor=no, readaxis="line", fixfile=" ",
              biassec="image", trimsec="image", zero="Zero", dark="",
              flat=" ", illum=" ", fringe="", minreplace=1., scantype="shortscan", nscan=1,
              interactive=no, function="chebyshev", order=4, sample="*", naverage=1,
              niterate=1, low_reject=3., high_reject=3., grow=1.)
          }
          else {
              if (h ==1 && i == 2) {
                  ccdproc ("flats/*.fit",
                  output=" ", ccdtype="flat", max_cache=0, noproc=no, fixpix=no, overscan=yes,
                  trim=no, zerocor=yes, darkcor=no, flatcor=no, illumcor=no, fringecor=no,
                  readcor=no, scancor=no, readaxis="line", fixfile=" ",
                  biassec="image", trimsec="image", zero="Zero", dark="",
                  flat=" ", illum=" ", fringe="", minreplace=1., scantype="shortscan", nscan=1,
                  interactive=no, function="chebyshev", order=4, sample="*", naverage=1,
                  niterate=1, low_reject=3., high_reject=3., grow=1.)
              }
              else { #(h == 1 && i == 1)
                  ccdproc ("flats/*.fit",
                  output=" ", ccdtype="flat", max_cache=0, noproc=no, fixpix=no, overscan=yes,
                  trim=yes, zerocor=yes, darkcor=no, flatcor=no, illumcor=no, fringecor=no,
                  readcor=no, scancor=no, readaxis="line", fixfile=" ",
                  biassec="image", trimsec="image", zero="Zero", dark="",
                  flat=" ", illum=" ", fringe="", minreplace=1., scantype="shortscan", nscan=1,
                  interactive=no, function="chebyshev", order=4, sample="*", naverage=1,
                  niterate=1, low_reject=3., high_reject=3., grow=1.)
              }
          }
      }


# ------------------------------------------------------------------------------------
# End of '1º 'ccdproc' run to apply 'biassec', 'trimsec' and 'zerocor''
# ------------------------------------------------------------------------------------



# ------------------------------------------------------------------------------------
# Flatcombine on flat frames
# ------------------------------------------------------------------------------------

      unlearn ccdproc
      ccdproc.output=" "
      ccdproc.ccdtype="flat"
      ccdproc.max_cache=0
      ccdproc.noproc=no
      ccdproc.fixpix=no
      ccdproc.overscan=no
      ccdproc.trim=no
      ccdproc.zerocor=no
      ccdproc.darkcor=no
      ccdproc.flatcor=yes
      ccdproc.illumcor=no
      ccdproc.fringecor=no
      readcor=no
      ccdproc.scancor=no
      ccdproc.readaxis="line"
      ccdproc.fixfile=""
      ccdproc.biassec="image"
      ccdproc.trimsec="image"
      ccdproc.zero="Zero"
      ccdproc.dark=""
      ccdproc.flat=" "
      ccdproc.illum=" "
      ccdproc.fringe=""
      ccdproc.minreplace=1.
      ccdproc.scantype="shortscan"
      ccdproc.nscan=1
      ccdproc.interactive=no
      ccdproc.function="chebyshev"
      ccdproc.order=4
      ccdproc.sample="*"
      ccdproc.naverage=1
      ccdproc.niterate=1
      ccdproc.low_reject=3.
      ccdproc.high_reject=3.
      ccdproc.grow=1.

      print ('')
      print ('Flatcombine on flat frames')
      print ('')

      if (k == 1) {
          flatcombine ("flats/flat*.fit",
          output="Flat", combine="average", reject="crreject", ccdtype="flat", process=yes,
          subsets=yes, delete=no, clobber=no, scale="mode", statsec="", nlow=1, nhigh=1,
          nkeep=1, mclip=yes, lsigma=3., hsigma=3., rdnoise=(rdnoise), gain=(gain),
          snoise="0.", pclip=-0.5, blank=1.)
      }
      else {
          if (k == 2){
              flatcombine ("flats/sky*.fit",
              output="Flat", combine="average", reject="crreject", ccdtype="flat", process=yes,
              subsets=yes, delete=no, clobber=no, scale="mode", statsec="", nlow=1, nhigh=1,
              nkeep=1, mclip=yes, lsigma=3., hsigma=3., rdnoise=(rdnoise), gain=(gain),
              snoise="0.", pclip=-0.5, blank=1.)
          }
          else {
              print ('No flat files found.') # Redundant control
              print ('Closing...')
              bye()
          }
      }

# ------------------------------------------------------------------------------------
# End of 'Flatcombine on flat frames'
# ------------------------------------------------------------------------------------



# ------------------------------------------------------------------------------------
# 2º 'ccdproc' run to apply 'flatcor'
# ------------------------------------------------------------------------------------

      print ('')
      print ('Ccdproc run to apply "flatcor"')

      if (h == 2 && i == 2) {
          ccdproc ("*.fit",
          output=" ", ccdtype="object", max_cache=0, noproc=no, fixpix=no, overscan=no,
          trim=no, zerocor=no, darkcor=no, flatcor=yes, illumcor=no, fringecor=no,
          readcor=no, scancor=no, readaxis="column", fixfile=" ",
          biassec="", trimsec="", zero="Zero", dark="",
          flat="Flat*", illum=" ", fringe="", minreplace=1., scantype="shortscan",
          nscan=1, interactive=no, function="chebyshev", order=4, sample="*",
          naverage=1, niterate=1, low_reject=3., high_reject=3., grow=1.)
      }
      else {
          ccdproc ("*.fit",
          output=" ", ccdtype="object", max_cache=0, noproc=no, fixpix=no, overscan=no,
          trim=no, zerocor=no, darkcor=no, flatcor=yes, illumcor=no, fringecor=no,
          readcor=no, scancor=no, readaxis="column", fixfile=" ",
          biassec=(biassec2), trimsec=(trimsec2), zero="Zero", dark="",
          flat="Flat*", illum=" ", fringe="", minreplace=1., scantype="shortscan",
          nscan=1, interactive=no, function="chebyshev", order=4, sample="*",
          naverage=1, niterate=1, low_reject=3., high_reject=3., grow=1.)
      }

      mkdir calib
      mv Flat*.* calib/   # Moves the Flat* files created by Flatcor and the Zero.out file created by Zerocor
      mv Zero.* calib/    # to the folder "calib/"

# ------------------------------------------------------------------------------------
# End of '2º 'ccdproc' run to apply 'flatcor''
# ------------------------------------------------------------------------------------



# ------------------------------------------------------------------------------------
# Rotate the frames
# ------------------------------------------------------------------------------------

print ('')
print('Rotating the frames: [-*,*]')

imtranspose (input="*.fit[-*,*]", output="*.fit")

# ------------------------------------------------------------------------------------
# End of 'Rotate the frames'
# ------------------------------------------------------------------------------------



print ('')
print ('Script finished correctly! :)')
print ('')

