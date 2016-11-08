procedure count ()

struct *file_var
struct *file_var2
struct *name

begin

    real line
    int m
    struct name2
    

	    if (! defpac ("ptools")) {
	        print ('\n This script must be loaded inside the package noao/digiphot/ptools')
	        bye()
	    }

    files ('*.coo.1', > 'temp_count.coo')
    file_var = 'temp_count.coo'

    while (fscan (file_var,name2) != EOF) { # en 'name' queda guardado el nombre del achivo .coo

        while (fscan (file_var,name) != EOF) { # en 'name' queda guardado el nombre del achivo .coo

		        print (name2)
				    m=0
				    while (fscan (name,line) != EOF ) { # en 'line' queda guardada la linea
				        m=m+1
				        print (m)
				    }
				    m = m-41 # en 'm' queda guardada la cantidad de estrellas en este archivo
				    print (name2, m, >> ('coo_file'))
		    }
    }
    del ('temp_count.coo')


    files ('*.mag.1', > 'temp_count.mag')
    file_var = 'temp_count.mag'
    while (fscan (file_var,name) != EOF) { # en 'name' queda guardado el nombre del achivo .mag
		    txdump.mode = 'hl' 
		    txdump.textfile = name
		    txdump.headers = no
		    txdump.fields = 'MSKY'
		    txdump.expr = 'MAG[1]!=INDEF'
		    txdump > auxiliar

		    file_var2 = ('auxiliar')
		    m=0
		    while (fscan (file_var2,line) != EOF) {
		        m = m + 1
		    }
		    del ('auxiliar')
		    print (name, m, >> ('mag_file'))
    }
    del ('temp_count.mag')


    files ('*.als.1', > 'temp_count.als')
    file_var = 'temp_count.als'
    while (fscan (file_var,name) != EOF) { # en 'name' queda guardado el nombre del achivo .als
		    m=0
		    while (fscan (name,line) != EOF ) { # en 'line' queda guardada la linea
		        m=m+1
		    }
		    m = m-44
		    m= m/2 # en 'm' queda guardada la cantidad de estrellas en este archivo
		    print (name, m, >> ('als_file'))
    }
    del ('temp_count.als')


end
