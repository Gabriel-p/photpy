
procedure select (image,xmin1,xmax1,ymin1,ymax1)

file image {prompt = "Input file name"}
real xmin1 {prompt = "Input min x value"}
real xmax1 {prompt = "Input max x value"}
real ymin1 {prompt = "Input min y value"}
real ymax1 {prompt = "Input max y value"}


begin

      string imname, expression
      real xmin,xmax,ymin,ymax
      
      imname = image
      xmin = xmin1
      xmax = xmax1
      ymin = ymin1
      ymax = ymax1
      
      if (! defpac ("daophot")) {
          print ('')
          print (' This script must be loaded inside the package \'noao.digiphot.daophot\'')
          bye()
      }

      expression = "XCENTER < "//xmin
      expression = expression // " "// '|| XCENTER > '//xmax
      expression = expression // " "// '|| YCENTER < '//ymin
      expression = expression // " "// '|| YCENTER > '//ymax

      print ('\n Expression: '//expression)

      pselect.infiles = imname
      pselect.outfiles = imname//'_1'
      pselect.expr = expression
      pselect.mode = "hl"
      pselect      

      rename.field = 'all'
      rename (imname, imname//'_original')
      
end

