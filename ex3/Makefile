
FC = g++
FCLAGS = -g -O3 -Wall -Werror -std=c++11

SOURCE = /Users/albertobassi/Documents/Molecular-Dynamics/functions.cpp main.cpp

all: main

main: $(SOURCE)
	$(FC) $(FCLAGS) -o $@ $(SOURCE) 

clean:
	-rm -f *.txt main

