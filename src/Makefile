
FC = g++
FCLAGS = -g -O3 -Wall -Werror -std=c++11


SOURCE = volume.cpp cell_volume.cpp sampler.cpp main.cpp

all: main

main: $(SOURCE)
	$(FC) $(FCLAGS) -o $@ $(SOURCE) 

clean:
	-rm -f *.txt main

