CC = g++
AR = ar

SOURCE = lammpsreader.cpp
HEADER = lammpsreader.h
OBJ = LAMMPSReader.o
TARGET = liblammpsreader.a

INSTALL_PATH = ~/lib/
INCLUDE_PATH =  ~/include/


.PHONY = clean install uninstall

default: lib

install: lib
	cp -v $(TARGET) $(INSTALL_PATH)
	cp -v $(HEADER) $(INCLUDE_PATH)

uninstall:
	rm -fv $(INSTALL_PATH)/$(TARGET)
	rm -fv $(INCLUDE_PATH)/$(HEADER)

clean:
	rm $(OBJ) $(TARGET)

lib: $(SOURCE) $(HEADER)
	$(CC) -c $(SOURCE) -o $(OBJ)
	$(AR) rcs $(TARGET) $(OBJ)
