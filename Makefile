CC = g++
AR = ar

SOURCE = lammpsreader.cpp
HEADER = lammpsreader.h
OBJ = LAMMPSReader.o
TARGET = liblammpsreader.a


.PHONY = clean

default: lib

clean:
	rm $(OBJ) $(TARGET)

lib: $(SOURCE) $(HEADER)
	$(CC) -c $(SOURCE) -o $(OBJ)
	$(AR) rcs $(TARGET) $(OBJ)
