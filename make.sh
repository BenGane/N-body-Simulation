# Compiles the command line simulation into an executable file - nbody
gcc -O0 src/nbody_main.c src/nbody.c src/utility.c -o nbody -Wall -Werror -Wvla -lm -lpthread

# Compiles the GUI into an executable file - nbody-gui
gcc -O0 src/nbodygui_main.c src/nbody.c src/utility.c -o nbody-gui -Wall -Werror -Wvla -lm -lSDL2 -lSDL2main -lSDL2_gfx -pthread 

# Compiles the test cases into an executable file - nbody-test
gcc -fsanitize=address test/nbody_test.c src/nbody.c src/utility.c -o nbody-test -Wvla -lm -lpthread -lcmocka