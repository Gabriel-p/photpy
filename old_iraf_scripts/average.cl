
procedure average ()

struct *file_var
struct *name


begin

     bool check
     struct imname
     real line, fwhm, stddev, sky, fwhm_a, stddev_a, sky_a
     int n, m



      files ('*_data', > 'temp.data')
      file_var = 'temp.data'
      fwhm_a=0
      stddev_a=0
      sky_a=0
      n=0
      print ('FWHM         STDDEV         SKY MEAN ', >> 'average')
      while (fscan (file_var,name) != EOF) {
          m=0
          while (fscan (name,line) != EOF ) {
              if (m==0) {
                  fwhm = line
#                  print ('FWHM: ', fwhm, >> ('average'))
              }
              else {
                  if (m==1) {
                      stddev=line
#                      print ('STDDEV: ', stddev, >> ('average'))
                  }
                  else {
                      sky=line
#                      print ('SKY MEAN: ', sky, >> ('average'))
                  }
              }
              m=m+1
          }
              print (fwhm,stddev,sky, >> 'average')
          fwhm_a = fwhm_a + fwhm
          stddev_a = stddev_a + stddev
          sky_a = sky_a + sky
          n = n+1
      }

      print ('\nAverage FWHM: ', fwhm_a/n, >> ('average'))      
      print ('Average STDDEV: ', stddev_a/n, >> ('average'))      
      print ('Average SKY MEAN: ', sky_a/n, >> ('average'))      
      print (' Average FWHM = '// fwhm_a/n)
      print (' Average STDDEV = '// stddev_a/n)
      print (' Average SKY MEAN = '// sky_a/n)

      del ('temp.data')


end
