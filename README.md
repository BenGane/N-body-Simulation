## How to build the application
- To build the application, simply execute `bash make.sh`. This will create three executable files `nbody`, `nbody-gui` and `nbody-test`.
- To run the command line simulation, execute `./nbody <iterations> <dt> (-b <bodies> | -f <filename>) <number of threads>`. Note that `<number of threads>` is an optional argument. If it isn't provided, the number of threads defaults to 1.
- To run the GUI simulation, execute `./nbody-gui <resolution_width> <resolution_height> <iterations> <dt> (-b <bodies> | -f <filename>) <scale> <number of threads>`. Note that I have implemented a scaling feature whereby a scaling value of 1 makes the data fit the bounds of the display perfectly. Also note that `<number of threads>` is an optional argument. If it isn't provided, the number of threads defaults to 1.
- To run the suite of test cases outlined in `nbody_test.c`, execute `./nbody-test`.

## Benchmark
- To run the benchmark script, first build the application by executing `bash make.sh`, then execute `benchmark.sh`. This runs the `nbody` command line simulation using a different number of threads and compares their execution times.

## Errors Encountered
- During the development of the project, I encountered an issue where my energy calculations were fluctuating significantly between iterations. The root cause was in my kinetic energy function where I was using incorrectly using positional values to compute `velocity^2`. 
- I also encountered an issue where I my multi-threaded program was updating the position of bodies using `velocity(t+1)` rather than `velocity(t)`. To counteract this issue, I synchronized the threads using a barrier. The barrier makes sure that the correct velocity values are used when each of the threads are updating the positions of their assigned bodies.
- The last issue I encountered was that some of my calculations were returning `-nan`, i.e. not a number. This was due to the fact that my distance function was in some cases returning zero, and thus other functions were attempting to divide values by zero. To fix this, I changed my distance function so that it returns `1e-10` instead of zero, in the case that the actual distance between the two bodies is zero.

## Step Function
- My step function spawns `N` threads, each of which computes and updates the velocities/positions of `N/(number of bodies)` bodies in the system.
- Before any thread can update the positions of the bodies, they all have to finish calculating and storing `velocity(t+1)`. To enforce this restriction,  my step function initialises a barrier which is passed to and used by the `thread_step` worker function.