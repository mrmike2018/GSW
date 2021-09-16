IFlags = -I$(HOME)/sw/include
LFlags = -L$(HOME)/sw/lib

all:
	g++ ${IFlags} GSW.cpp -O3 -o GSW_output ${LFlags} -lntl -lgmp -lstdc++
	
run:
	./GSW_output