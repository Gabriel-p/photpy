
from imexam.imexamine import Imexamine
import numpy as np

plots = Imexamine()
data = np.random.rand(100,100) * np.ones((100,100))
plots.plot_line(35,45,data) #shows a matplotlib window with a plot
plots.save() #saves the current plot to file

# # You can also set the data attribute of the plots object and then just call many plots without specifying the data again:

# plots.set_data(data)
# plots.plot_line(35,45)
