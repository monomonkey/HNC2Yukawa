SHELL := /bin/bash
c_compiler = gcc     #Select the C compiler

### Compiling options ###
gsl = -lgsl -lgslcblas -lm
gflags = -ffast-math

### Variables ###
# Math #
math_aux = math_aux


### Directories ###
mainFolder = ./facdes2Y
mainFunc_dir = $(mainFolder)/mainFunc/
objects_dir = $(mainFolder)/objects/
structs_dir = $(mainFolder)/structures/
math_dir = $(mainFolder)/math/
main_dir = $(mainFolder)/mainFiles/

## Files ##
structs = structures
math_aux = math_aux
mainFunc = facdes2Y
mainFile = mainFile

### Libraries ###
math_libs = $(math_aux)


#All dependencies needed to create the test.exe file
all_dep = $(main_dir)$(mainFile).o $(main_dir)lib$(mainFunc).a\
$(objects_dir)$(mainFunc).o $(objects_dir)$(structs).o $(objects_dir)$(math_aux).o\

# Dependencies
structures_dep = $(structs_dir)$(structs).c $(structs_dir)$(structs).h
math_dep = $(math_dir)$(math_aux).c $(math_dir)$(math_aux).h
mainFunc_dep = $(mainFunc_dir)$(mainFunc).c $(mainFunc_dir)$(mainFunc).h


main.exe: $(all_dep)
	gcc $(main_dir)$(mainFile).o -L$(main_dir). -l$(mainFunc) -o $@ $(gsl) $(gflags)

#Building main.o
$(main_dir)$(mainFile).o:$(mainFile).c
	gcc -c $(mainFile).c -o $@ -L$(main_dir). -l$(mainFunc) $(gsl) $(gflags)

#Building allIn library
$(main_dir)lib$(mainFunc).a:$(objects_dir)*.o
	ar rcs $@ $^

#Building mainFunction library
$(objects_dir)$(mainFunc).o: $(mainFunc_dep)
	gcc -c $(mainFunc_dir)$(mainFunc).c -o $@ $(gsl) $(gflags)

#Building structural libraries
$(objects_dir)$(structs).o: $(structures_dep)
	gcc -c $(structs_dir)$(structs).c -o $@ $(gsl) $(gflags)


#Building math libraries
$(objects_dir)$(math_aux).o: $(math_dep)
	gcc -c $(math_dir)$(math_aux).c -o $@ $(gsl) $(gflags)


clean:
	rm $(objects_dir)*.o
	rm ./*.dat
	rm ./*.exe

