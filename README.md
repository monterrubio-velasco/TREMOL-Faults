# TREMOL-Faults

This code is developed to simulate the aftershocks rupture. Our model considers the fault active systems as regions likely to fail. The model generates seismic statistical patterns than those observed in real series.
Further details are described in Monterrubio-Velasco, M., Zúñiga, F. R., Carrasco-Jiménez, J. C., Márquez-Ramírez, V., & de la Puente, J. (2019). Modeling active fault systems and seismic events by using a fiber bundle model–example case: the Northridge aftershock sequence. Solid Earth, 10(5), 1519-1540.

In this folder, we add the python codes and one input file to be the example "MatrizNorthridgeFaults180.dat".
This file contains the fault information of the Northridge active fault systems. The image is processed to be read and mapped as binary values to the model domain (see Fig. 1 of the paper)
The meaning of "180" in the input file name is because we are working in a square domain of 180 per 180 cells. This code is prepared to simulate square domains. We are working to generalized the dimensions to use in any domain shape.

To execute the code you simply run in the console the file:  

>> python Preprocessing.py

and the code starts to run ... 


You can play with different values of VectorP (P parameter in the paper), VectorPiFrac (\pi_{frac} parameter in the paper), VectorPiBG (\pi_{bkg} parameter in the paper)


We are working on optimizing the Python model version. At current status, the python version is also coded to work in an embarrassingly parallel workload.  Moreover, we have a Julia version code that runs faster than the Python version.

