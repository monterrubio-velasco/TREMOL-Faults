# TREMOL-Aftershocks and Faults systems

This code is developed to simulate the aftershocks rupture. Our model considers the fault active systems as regions likely to fail. The model generates seismic statistical patterns than those observed in real series.
Further details are described in Monterrubio-Velasco, M., Zúñiga, F. R., Carrasco-Jiménez, J. C., Márquez-Ramírez, V., & de la Puente, J. (2019). Modeling active fault systems and seismic events by using a fiber bundle model–example case: the Northridge aftershock sequence. Solid Earth, 10(5), 1519-1540.

In this folder, we add the python codes and one input file to be the example "MatrizNorthridgeFaults180.dat".
This file contains the fault information of the Northridge active fault systems. The image is processed to be read and mapped as binary values to the model domain (see Fig. 1 of the paper)
The meaning of "180" in the input file name is because we are working in a square domain of 180 per 180 cells. This code is prepared to simulate square domains. We are working to generalized the dimensions to use in any domain shape.

To execute the code you simply run in the console the file:  

>> python Preprocessing.py

and the code starts to run ... 

You can set with different values of VectorP (P parameter in the paper), VectorPiFrac (\pi_{frac} parameter in the paper), VectorPiBG (\pi_{bkg} parameter in the paper)

As output files we generated three differents datasets, in this example we have:

1) "VecPosi.dat": Is the initial load matrix considering the P parameter value. IN the example the size of this matrix is 180. This file is computed in the script Preprocessing.py (line 52) 

2) "vecWriteMagniHB1.dat": Computes the magnitude (row 1) and the cumulative number of events with magnitude larger or equal than that magnitude (row 2). This file is computed in the script calcuMagniInKM.py (line 78)

3) "VecDatosNR-180.0-0.0-0.9-0.65.dat": Is the raw data coming from the FBM algorithm. By columns the information provided is:
column 1: number of steps k
column 2: accumulated time, T_k
column 3: inter-event time, \delta_k
column 4: the flag that indicates if is normal or avalanche type rupture (0 or 1)
column 5: the sum of loads in all the system per each k-step
column 6: the sum of loads at the power \rho in all the system per each k-step
column 7: the load value of the selected cell that fails at each k-step
column 8: coordinate of the selected cell to fail in X_axis
column 9: coordinate of the selected cell to fail in Y_axis
column 10: the flag that indicates if the failure cell belongs to a fault or not in the domain. 
ignoring in this test

4) "CI_Northridge-180Test-FBMfaults.csv" This file contains the post-processing statistical analysis that emerges from the raw catalog. 

column 1: P
column 2: fhiFractura
column 3: Nbox
column 4: NumCeldasFrac
column 5: NumCeldasNonFrac

From column 5 to 16 is the information related with the Correlation Fractal Dimension (D0, D1, D2) computed in the script "CorreDimensionAvalanchesV11.py" 
  
From column 17 to 34 the parameters related with the magnitude are computed in the script "calcuMagniInKM.py"
From column 36 to 51 the parameters related with Hurst exponent are computed in the script "ParametrosResultantes.py"
From column 62 to 79 the parameters related with Omori-Utsu law are computed in the script "OmoriFit.py". The time series used to obtain these parameters are computed from  the script "separacionLACAS.py". Further details of this separation method can be read in "Monterrubio-Velasco et al., 2016. Two complementary stress release processes based on departures from Omori’s law. Geosciences Journal, 20(1), 41-55". 
From column 86 to 102 other values ara depicted, some of the most relevants are:

column 86: number of elements of the Magnitude series
column 89: area of each cell expresed in km²
   
  

We are working on optimizing the Python model version. At current status, the python version is also coded to work in an embarrassingly parallel workload.  Moreover, we have a Julia version code that runs faster than the Python version.

