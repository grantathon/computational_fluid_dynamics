WIthin the currenty directory you have the project directory and CFD_documentation.pdf.  The project directory contains the final version of our project as it operates in its natural habitat--a CPU cluster.  The project can also be run on a single node, but it runs much faster on a CPU cluster.

Since the project's natural habitat is on a CPU cluster (e.g., lxlogin1.lrz.de), the following procedure MUST be implemented on that very same linux cluster with all of the contents in the project directory:

    1)  Copy project directory to lxlogin1.lrz.de and login.
        a)  For example => scp -r project/ [YOUR_USERNAME]@lxlogin1.lrz.de:
    2)  Enter the ns_sim_serial directory and make.
    3)  Copy over the generated "sim" executable to the monte_carlo directory.
    4)  Enter the monte_carlo directory and make.
    5)  Allocated space on the cluster.
        a)  For example => salloc --ntasks=8
    6)  After space has been allocated, run the program as explained in the submitted documentation (i.e., CFD_documentation.pdf).
        a)  For example => mpirun -np 8 ./sim_UQ 1 1 0 1 16 0.01225 0.00245 100 20