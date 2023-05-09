
#########################
### To generate the figures of the article we need obtain the results of the following tests:

# TEST_0 ---------------
This test generates the set of time series (of fixed size) with different number of changes (in different dimensions p) and checks the identically of results of FPOP and PELT for these simulations. 

For each dimension we write the results in the table, where there are two columns: the first - the number of segments, the second - TRUE, If the results of FPOP and PELT are identically, otherwise FALSE. These tables we also save in the files. 

The results of our testing are available in the archive “TEST_RESULTS.zip”.


# TEST_1 ---------------
This test generates the set of time series (of fixed size) in dimensions p =2,..,10 and evaluates the number of candidates of change stored over time by GeomFPOP with type R and S.

We write the results in a table with two columns: the first - time t, the second - the number of change point candidates at time t.

We also save these tables in files. 

These files are available in archive “TEST_RESULTS.zip”.


# TEST_2 ---------------
This test generates the set of time series of different size in dimensions p =2,3,4 and evaluates the run time of GeomFPOP with the types R and S.

For each dimension and for each type of GeomFPOP we write the results in the tables with two columns: the first - number of data points, the second - the work time of GeomFPOP. 

We also save these tables in files. 

These files are available in archive “TEST_RESULTS.zip”.


# TEST_3 and TEST_4 ---------------
These tests generate the set of time series of different size in dimensions p =2,..,10,100 and evaluate the run time of GeomFPOP(R-type random/random) and PELT.

For each dimension we write the results in the tables with two columns: the first - number of data points, the second - the work time of algorithm.

We also save these tables in files. 

These files are available in archive “TEST_RESULTS.zip”.


# TEST_5 ---------------
This test generates the set of time series of fixed size with changes in dimensions p =2,3,4 and evaluates the run time of GeomFPOP (R-type random/random) and PELT. 

For each dimension we write the results in the tables with two columns: the first - number of segments, the second - the work time of algorithm. 

We also save these tables in files. These files are available in archive “TEST_RESULTS.zip”.


# TEST_6 ---------------
This test generates the set of time series (of fixed size) in dimensions p =2,3,4 and evaluates the pruning efficiency of GeomFPOP with different update rules.

We write the results for each dimension and for each optimization strategy in a table with two columns: the first - time t, the second - the number of change point candidates at time t. 

We also save these tables in files. 

These files are available in archive “TEST_RESULTS.zip”.


# TEST_7 ---------------
This test generates the set of time series of different size in dimensions p =2,3,4 and evaluates the run time of GeomFPOP with different update rules.

We write the results for each dimension and for each optimization strategy in a table with two columns: the first - number of data points, the second - the work time of GeomFPOP. 

We also save these tables in files. 

These files are available in archive “TEST_RESULTS.zip”.



#########################
### Remark 1:  Attention, TEST_5 requires a long time to complete (more than a day).

### Remark 2: A detailed description of the tests can be found in the remark section of the corresponding files “*.R”

### Remark 3: To economize your time, you can use the archive “TEST_RESULTS.zip” in the folder “FIGURES and Test results”. In this archive there are the results of all tests, described above.

