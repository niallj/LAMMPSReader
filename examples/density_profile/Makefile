CC = g++ -Wall --std=c++0x
AR = ar

EXEC = density_profile

LIBRARY_PATH = /home/niall/lib/
INCLUDE_PATH =  /home/niall/include/

.PHONY = clean

all: density_profile.cpp
	$(CC) -I$(INCLUDE_PATH) -L$(LIBRARY_PATH) -o $(EXEC) density_profile.cpp -llammpsreader

clean:
	rm -fv *.o
	rm -fv density_profile
