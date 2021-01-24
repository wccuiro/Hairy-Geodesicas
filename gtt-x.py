import black_hole as bh
import matplotlib.pyplot as pyplot
import numpy as np

miBH = bh.black_hole(1, 12.52655373, 1.4)
print(miBH.horizonte_hairy_x(), miBH.horizonte_hairy(), miBH.horizonte())