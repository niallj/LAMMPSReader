/*density_profile.cpp
  Intended as a simple example of the LAMMPSReader class
  LAMMPSReader is a wrapper to extract data from LAMMPS dump files, text and binary.
  Copyright (C) 2013 Niall Jackson <niall.jackson@gmail.com>
  This program contains no LAMMPS source code.
  More LAMMPS information may be found at http://lammps.sandia.gov

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along
  with this program; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#include <fstream> //provides file I/O
#include <iostream> //provides terminal I/O
#include <sstream> //provides string->int conversion
#include <string> //provides string comparison functions

#include <lammpsreader.h> //declares the LAMMPSReader class

//the LAMMPSReader class is kept in the LAMMPSReaderNS "namespace"
//we bring this namespace into our current scope to save typing
//it means we can just say LAMMPSReader, instead of LAMMPSReaderNS::LAMMPSReader
using namespace LAMMPSReaderNS;

//declaration of our callback object 
//LAMMPSReader uses Callbacks
//When ReadFrame() is called to read a timestep, there are four key events.
//Each of these events triggers a call to one of the functions in the Callback
//class. You must define an object which inherits from the Callback class, and
//define the functions for those events that you want to implement.
//the four events are as follows:
// Header of timestep read: call StartOfTimestep(). Box dimensions, timestep number,
//                          and number of atoms are available at this stage
//                          call BoxBounds() to give box information directly to the user
// Each time an atom is read: call AtomLine with the populated AtomData struct
// End of current frame reached: call EndOfTimestep()
class Cback : public Callback {
private:
  int axis;
  int nbins;
  int ntimesteps;
  
  int *bins;
  double bin_width;
  double bin_volume;
  
  bool memory_allocated;
  double lo_offset;
public:
  Cback(int, int);
  ~Cback();
  void AtomLine(const AtomData& ad, LAMMPSReader *lr);
  void BoxBounds(char boundaries[3][2], double lo[3], double hi[3]);
  void StartOfTimestep(LAMMPSReader* lr);
  void EndOfTimestep(LAMMPSReader* lr);
  void print_histogram();
};

//Constructor: special function called whenever a new object of this type
//is created
Cback::Cback(int a, int b) : axis(a), nbins(b), ntimesteps(0), memory_allocated(false), lo_offset(0) {};

//Destructor: called whenever an object of this type is destroyed
//mainly for giving back any memory we asked format
Cback::~Cback() {
  if(memory_allocated) {
    delete[] bins;
  }
}

//definition of the AtomLine hook
//ad is an AtomData struct. It is an object containing atom properties.
//ad.x, ad.vx, ad.id, etc.
//Obviously, only those properties requested from the file are provided.
//lr is a pointer back to the LAMMPSReader object that called this function
//this is useful is we want to access the number of atoms, or some such
void Cback::AtomLine(const AtomData& ad, LAMMPSReader *lr) {
  //select the relevant component of the position vector
  double pos= 0;
  if(axis == 0) {
    pos= ad.x;
  } else if(axis == 1) {
    pos= ad.y;
  } else if (axis == 2) {
    pos= ad.z;
  }
  
  //shift this so that it starts at zero
  pos-= lo_offset;
  //work out which bin we're in
  int i= pos/bin_width;
  //increment the total in that bin
  bins[i]++;
}

//definition of the BoxBounds hook
//boundaries contains the boundary states
//boundaries[0][0] is the lower x boundary, boundaries[2][1] is the upper z boundary
//same values as in the lammps manual: p for periodic, f for fixed, etc.
//if a boundary is periodic, LAMMPSReader will automatically wrap the co-ordinates
//this can be switched off by setting LAMMPSReader::wrap to false.
//lo contains the lower boundaries (LAMMPS boxes do not need to start at zero)
//hi contains the upper boundaries
void Cback::BoxBounds(char boundaries[3][2], double lo[3], double hi[3]) {
  //if this is the first time we've run, use the box dimensions to work
  //out how many bins we need for the density histogram
  //if the box dimensions change during the run, then the number of bins
  //needed may change. THIS MAY CAUSE A SEGFAULT. A proper program 
  //would handle this, but this is only a simple example.
  if(!memory_allocated) {
    bins= new int[nbins];
    bin_width= (hi[axis]-lo[axis])/nbins;
    bin_volume= (hi[0]-lo[0])*(hi[1]-lo[1])*(hi[2]-lo[2])/nbins;
    lo_offset= lo[axis];
    memory_allocated= true;
  }
}

//definition of the StartOfTimestep hook
//lr is a pointer back to the LAMMPSReader object that called this function
//this is useful is we want to access the number of atoms, or some such
void Cback::StartOfTimestep(LAMMPSReader* lr) {};

//definition of the EndOfTimestep hook
//lr is a pointer back to the LAMMPSReader object that called this function
//this is useful is we want to access the number of atoms, or some such
void Cback::EndOfTimestep(LAMMPSReader* lr) {
  //increment the number of timesteps so that we can
  //normalise by it
  ntimesteps++;
};

//prints out the density histogram at the end
void Cback::print_histogram() {
  //creates an ofstream (output file stream)
  //will be created if non existent, overwritten otherwiser
  std::ofstream o("density.dat");
  
  for(int i= 0; i < nbins; i++) {
    //write the centre of the bin, and the density
    o << (i+0.5) * bin_width << " " << bins[i]/(bin_volume*ntimesteps) << std::endl;
  }
  
  o.close();
}


// entry point
// argc is the number of command line arguments
// argv is an array of C-style strings containing the arguments
// e.g. if the command given is /home/niall/density_profile dump.lammpstrj x 50
// then argc = 4
// argv[0] = /home/niall/density_profile
// argv[1] = dump.lammpstrj
// argv[2] = x
// argv[3] = 50
int main(int argc, char **argv) {
  //we access the standard output through std::cout.
  //the << operator means "redirect the string to std::cout"
  //std::endl is the newline character
  std::cout << "Density profile calculator" << std::endl;
  
  //we need at least four arguments
  if(argc < 4) {
    //we access the standard error through std::cerr
    std::cerr << "Usage: " << argv[0] << " <path to trajectory> <axis along which to draw profile> <number of bins>" << std::endl;
    std::cerr << "<axis> should be x, y or z." << std::endl;
    
    //terminate the program with non-zero return value (error)
    return 1;
  }
  
  //work out which axis we want
  //create a string object out of argument 2
  std::string axis(argv[2]);
  int iaxis= 0;
  if(axis.compare("x") == 0) {
    iaxis= 0;
  } else if(axis.compare("y") == 0) {
    iaxis= 1;
  } else if (axis.compare("z") == 0) {
    iaxis= 2;
  } else {
    std::cerr << "Unrecognised axis. Should be x, y or z." << std::endl;
    return 1;
  }
  
  //get the number of bins required
  //the "stringstream" is used to convert the string containing the number of bins into an integer
  std::stringstream ss(argv[3]);
  int bins;
  //if for some reason this fails, throw an error.
  if(!(ss >> bins)) {
    std::cerr << "For some reason, " << argv[3] << " is not an acceptable number of bins. It must be convertible to an integer!" << std::endl;
    return 1;
  }
  
  //create a new LAMMPSReader object
  LAMMPSReader lr;
  
  //try to open the file named by argv[1] (the first command line argument)
  //the second argument specifies whether or not this is a binary file
  //for the sake of example, we'll say that it is
  //the second argument has a default value of false
  //so if working with a text file, argument 2 is not required
  //the return value is a boolean: true on success, false on failure
  bool result= lr.open(argv[1], true);
  
  //if we failed to open the file for some reason
  if(!result) {
    //print an error message
    std::cerr << "Failed to open dump file (" << argv[1] << "). Stopping." << std::endl;
    //terminate with an error
    return 1;
  }
  
  //create an instance of the Cback object that we defined earlier
  //we have to create a pointer to it for the interface with LAMMPSReader to work
  Cback *c= new Cback(iaxis, bins);
  
  //ReadFrame reads an entire timestep of data from the file. It returns true if the read
  //was successful, and false otherwise. This line can be read as "while there are still more 
  //timesteps to read, read the next timestep".
  //
  //the behaviour of argument 1 is different for binary and text files
  //
  //for BINARY files
  //argument one is a string which specifies EVERY atom field in the dump file, in the correct order
  //e.g., if your dump command is dump d1 all custom 1 out.lammpstrj.bin id type mol x y z, then 
  //argument one MUST be "id type mol x y z". The LAMMPS binary format is not "self-labelling", like
  //the text format is.
  //
  //for TEXT files
  //
  //argument one is a string which specifies only those fields that you wish to extract from the dump
  //file, in any order. e.g. if your dump command is dump d1 all custom 1 out.lammpstrj.bin id type mol x y z,
  // and you want to read the fields id, type, and z, argument one should be "id type z".
  //
  //argument 2 specifies the Callback object to be used when reading the file (see the explanation of the Callback class
  //above)
  
  int i= 0; //loop counter
  while(lr.ReadFrame("id type x y z", c)) {
    i++; //increment the count of the total timesteps read
  }
  
  //once all timesteps are read, call the function 
  //that prints out the histogram
  c->print_histogram();
  
  std::cout << "Read " << i << " timesteps." << std::endl;
  
  //get rid of the Cback object that we created
  delete c;
}