################################################################################
#
#        ===========================================================
#                   OUTPUT ONLY STARS WITH NO INDEF VALUES
#        ===========================================================
# 
#     To run properly this procedure must be executed inside the package:
#                           'noao.digiphot.ptools'
#
#                            by Gabriel Perren 2009
################################################################################


procedure clean ()

real dummy 

begin

	    if (! defpac ("ptools")) {
	    print ('')
	        print (' This script must be loaded inside the package noao/digiphot/ptools')
	        bye()
	    }
	    else { # Do nothing
	    }

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
      
end
      
      
                  
