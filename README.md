
### Updated R code and data from Gasparrini Biomet 2017

--------------------------------------------------------------------------------

R code reproducing the results published in the article:

Gasparrini A, Scheipl F, Armstrong B, and Kenward MG. A penalized framework for distributed lag non-linear models. Biometrics. 2017;73(3):938-948.[[freely available here](http://www.ag-myresearch.com/2017_gasparrini_biomet.html)]

This methodological article describes the extension of distributed lag linear and non-linear models (DLMs and DLNMs) thorugh generalized additive models via penalized splines. The methodology is implemented by embedding functions in the R packages [dlnm](https://github.com/gasparrini/dlnm) and mgcv. The code reproduces two examples of application in time series and survival analysis, respectively, and the results of the simulation study.

The software implementation and the use of the functions available in the [R package dlnm](https://github.com/gasparrini/dlnm) are described in detail in the vignette that accompanies the package (also available at the [CRAN page](https://cran.r-project.org/web/packages/dlnm/index.html))

--------------------------------------------------------------------------------

The code is provided in three different folders:

  * **example1** includes the code and data for replicating the first example illustrating an application in time series analysis:
    * *london.csv* stores the daily time series data from London in the period 1993--2006
    * *ex1_01.ext.R* performs the analysis using the external method
    * *ex1_01.int.R* performs the analysis using the internal method
    * *ex1_02.plots.R* reproduces the plots included in the article

  * **example2** includes the code and data for replicating the second example illustrating an application in survival analysis:
    * *uminers.csv* stores the data from the Colorado Plateau uranium miners cohort, including individual information for 3,347 male subjects
    * *ex2_01.prep.R* prepares the dataset
    * *ex2_02.ext.R* perform the analysis using the external method
    * *ex2_02.int.R* perform the analysis using the internal method
    * *ex2_03.plots.R* reproduces the plots included in the article

  * **simul** reproduces the simulation study
    * *sim_00.prep.R* simulates and prepares the data
    * *sim_01.ext.R* performs the simulations using the external method
    * *sim_01.int.R* performs the simulations using the internal method
    * *ex2_03.plots.R* summarizes and plots the results
  
Download as a ZIP file using the green button *Clone or download* above
