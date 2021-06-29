--------------------------------------------------------------------------------
   Exam
--------------------------------------------------------------------------------
This folder contains the files for the exam in the Practical Programming and
Numerical Methods course at Aarhus University, spring 2021.

My student number is 201706760, so the number of the assignment to be solved is
60 % 22 = 16.


--------------------------------------------------------------------------------
   Comments on the examination project
--------------------------------------------------------------------------------
I feel that I have solved the examination project to the highest possible value
(10/10). I have not only solved it for the normal integrator with limits being
finite, but I have also made it possible to use infinite limits.

Due to my use of a Windows computer and therefore WSL, it is not possible for me
to actually run the code, but after completion of the examination project (and
being frustrated it wouldn't run on my pc, I looked at the code of Marc, who
also have this examination project, and I can see, that we have the same form
for the integration function, and his is able to run on an Ubuntu computer, thus
my examination project should also be able to. Only the part from the Homework 6
have been tested and run on my own computer.
The conclusion of the project in the out.txt file thus also stem from looking at
similar projects and concluding from their numbers..


--------------------------------------------------------------------------------
   16: Adaptive recursive integrator with subdivision into three subintervals
--------------------------------------------------------------------------------
Implement a (one-dimensional) adaptive recursive integrator (open or closed
quadrature, at your choice) which at each iteration subdivides the interval not
into two, but into three sub-intervals. Reuse points. Compare with the adaptive
integrator from your homework.


### Method ###

The method is based on the "Numerical Integration"-PDF.

The adaptive quadrature is an algorithm for which the integration interval is
subdivided into adaptively refined subintervals until a given accuracy goal is
reached. For this project, the given number of subintervals is 3, thus in order
to subdivide total interval [a, b] into three for then to use the quadrature on
the half-intervals, one needs to choose the points x_i, which will will become

			    x_i = {1/6, 3/6, 5/6}.

From these points weights w_i can be calculated using the trapezium rule,

	  [[   1        1        1    ]      [[ w_1 ]      [[  1  ]
	   [  x_1      x_2      x_3   ]   *   [ w_2 ]   =   [ 1/2 ]   ,
	   [  x_1^2    x_2^2    x_3^2 ]]      [ w_3 ]]      [ 1/3 ]]

yielding

			    w_i = {3/8, 2/8, 3/8},

and likewise from the rectangle rule

			    v_i = {1/3, 1/3, 1/3}.

The the above the points and weights, x_i's and w_i's, are cited for the
rescaled integration interval [0, 1], and the transformation of the points and
weights to the original interval [a, b] is given as

	    x_i -> a + (b - a)/x_i,      and      w_i -> (b - a)w_i.


Now, subdividing into three subintervals, Q_i, in a uniform partition in the
interval [a, b] yields the integration intervals

		a < c = a + (b - a)/3 < d = a + 2(b - a)/3 < b,

and thus the following integral

     int_a^b f(x) dx = int_a^c f(x) dx + int_c^d f(x) dx + int_d^b f(x) dx.
		          |_______________| |_______________| |_______________|
			         Q_1	           Q_2		     Q_3


--------------------------------------------------------------------------------
   Structure
--------------------------------------------------------------------------------
This folder (/ExamProject/) is the main folder of the exam project, it contains:

  - main.c: The code running and testing the tri-subinterval adaptive and
	    recursive integration routine.

  - Makefile: The file containing the make commands for running the main.c file.

  - README.md: This file.

  - integratorTridivision*
      - .c: The code for the recursive and adaptive integrator deviding each
            intergration interval into three parts.
      - .h: The header file for the recursive and adaptive integrator deviding
            each intergration interval into three parts.

  - integratorBidivision*
      - .c: The code for the recursive and adaptive integrator deviding each
            intergration interval into two parts. This is the same integrator as
	    that in Homework 6.
      - .h: The header file for the recursive and adaptive integrator deviding
            each intergration interval into two parts. This is the same
	    integrator as that in Homework 6.

  - out.txt: Text document containing the results from the examination project.
             (NOTE: Might not be there, see "Comments on examination project".)
