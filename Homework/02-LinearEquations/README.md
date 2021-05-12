Using the command 'make' the program runs main liniarEquationSolver.c, writes out the answers and checks, generates timeDependence.png and opens this using the computer's default image viewer.

The answers and checks during the parts of the exercise can be found in 'out.txt'. 'data.txt' contains the dimensions and the runtimes, which are plottet on the graph 'timeDependence.png'.

**NB!** One may **never** just remove 'data.txt' from the directory, since the default make command expects 'data.txt' to be present in the directory and since there are no changes to linearEquationSolver.c, this is not run, thus not generating a new 'data.txt' file. Should one need to remove 'data.txt' one can use 'make clean' thus also cleaning 'out.txt' such that this shall be updated. Thus if one wants new times in 'data.txt' one should use 'make clean' first and then 'make'.
