
procedure star_average ()

struct *file_var
struct *name


begin

     real line 
     int n, m
     

      files ('*.als.1', > 'temp.data')
      file_var = 'temp.data'
      n=0

      while (fscan (file_var,name) != EOF) {
          m=0
          while (fscan (name,line) != EOF ) {
              m=m+1
          }
          print ((m-44)/2, >> 'star_average')
          n=n+(m-44)/2
      }
      print ('', >> 'star_average')
      print (n/8, >> 'star_average')
      
      del ('temp.data')

     
end
