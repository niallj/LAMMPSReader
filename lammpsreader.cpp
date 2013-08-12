/*
    lammpsreader.cpp
    LAMMPSReader is a C++ class for reading LAMMPS dump files
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

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "lammpsreader.h"

namespace LAMMPSReaderNS {
  std::vector<std::string> explode(std::string s) {
    //tokenize the string
    std::string sub;
    std::vector<std::string> v;
    size_t tok_pos= s.find(" ");
    while(tok_pos != std::string::npos) {
      sub= s.substr(0, tok_pos);
      s= s.substr(tok_pos+1, std::string::npos);
      if(!sub.empty()) {
	v.push_back(sub);
      }
      tok_pos= s.find(" ");
    }
    if(!s.empty()) {
      v.push_back(s);
    }
    return v;	
  }
	
  LAMMPSReader::LAMMPSReader() {
    //initialise the variables
    last_tstep= -1;
    n_atoms= 0;
    for(int i= 0; i < 3; i++) {
      box_lo[i]= 0.0;
      box_hi[i]= 0.0;
      //mark the boundaries as u, for unset, for now
      boundaries[i][0]= 'u';
      boundaries[i][1]= 'u';
    }
  }
  
  LAMMPSReader::~LAMMPSReader() {
    close();
  }

  bool LAMMPSReader::open(const std::string& filename, bool bin) {
    if(file.is_open()) {
      file.close();
    }
    if(bin) {
      file.open(filename.c_str(), std::ios::binary);
    } else {
      file.open(filename.c_str());
    }
    if(!file.is_open()) {
      std::cerr << "Error! Failed to open file " << filename << std::endl;
      return false;
    }
    curfile= filename;
    binary= bin;
    return true;
  }

  void LAMMPSReader::close() {
    if(file.is_open()) {
      file.close();
    }
    curfile= "";
  }

  bool LAMMPSReader::ReadFrame(const std::string& s, Callback *c) {
    //if this is a text file, s tells us which properties the user
    //wants us to extract from the file
    //if this is a binary file, s tells us ALL of the properties in the dump file
    //s tells us what properties the user wants us to extract from the file
    std::vector<std::string> args= explode(s);
    if(!file.is_open()) {
      std::cerr << "LAMMPSReader::ReadFrame() called while no file is open." << std::endl;
      return false;
    }
    if(binary) {
      //this is a binary file, which is handled a little differently
      return ReadBinaryFrame(args, c);
    }
    if(file.eof()) {
      //we're already at the end, so no more to read
      return false;
    }
    std::string line;
    bool insideTstep= false;
    //the map is sorted by key, which makes output look awkward
    //because it's not in the same order as in the file.
    //store a vector as well, which are sorted in the order by which values were pushed on
    std::map<std::string, int> columns;
    std::vector<std::string> avail_columns;
    std::ifstream::streampos line_start= file.tellg();
    while(std::getline(file, line)) {
      //process any information about the frame
      //tokenize the string
      std::vector<std::string> v= explode(line);
      if(v[0].compare("ITEM:") == 0) {
	if(v[1].compare("TIMESTEP") == 0) {
	  //the next line contains the current timestep
	  if(insideTstep) {
	    //but we're already in a timestep, so seeing this line means we've hit the end of the timestep
	    //go back one line in the file, then return from this function
	    c->EndOfTimestep(this);
	    file.seekg(line_start);
	    return true;
	  } else {
	    c->StartOfTimestep(this);
	    insideTstep= true;
	  }
	  if(std::getline(file, line)) {
	    last_tstep= atoi(line.c_str());
	  } else {
	    std::cerr << "ERROR: Failed to read a timestep after an ITEM: TIMESTEP line. (" << curfile << ")" << std::endl;
	    return false;
	  }
  } else if(v[1].compare("NUMBER") == 0) {
	  //the next line contains the number of atoms
	  if(std::getline(file, line)) {
	    n_atoms= atoi(line.c_str());
	  } else {
	    std::cerr << "ERROR: Failed to read a timestep after an ITEM: TIMESTEP line. (" << curfile << ")" << std::endl;
	    return false;
	  }				
	} else if(v[1].compare("BOX") == 0) {
	  //the remaining 3 tokens on this line specify the nature of the boundaries
	  if(v.size() < 6) {
	    std::cerr << "ERROR: Malformed ITEM: BOX BOUNDS line. Expected 6 tokens, only found " << v.size() << ". The offending line is: " << std::endl;
	    std::cerr << line <<  " (" << curfile << ")" <<std::endl;
	    return false;
	  } else {
	    for(int i= 0; i < 3; i++) {
	      boundaries[i][0]= v[i+3][0];
	      boundaries[i][1]= v[i+3][1];
	    }
	  }
	  //the next 3 lines contain box dimensions
	  for(int i= 0; (i < 3) && std::getline(file, line); i++) {
	    //we should have two tokens in each line, a lower bound and an upper bound
	    std::vector<std::string> tokens= explode(line);
	    if(tokens.size() < 2) {
	      std::cerr << "ERROR: Malformed box bounds line. Expected 2 tokens, only found " << v.size() << ". The offending line is: " << std::endl;
	      std::cerr << line  << "(" << curfile << ")" << std::endl;
	      return false;
	    } else {
	      box_lo[i]= atof(tokens[0].c_str());
	      box_hi[i]= atof(tokens[1].c_str());
	    }
	  }
	  c->BoxBounds(boundaries, box_lo, box_hi);
	} else if(v[1].compare("ATOMS") == 0) {
	  //the remaining tokens on this line tell us what data we're going to get
	  for(int i= 2; i < (int)v.size(); i++) {
	    columns[v[i]]= i - 2;
	    avail_columns.push_back(v[i]);
	  }
	  //now check that the user didn't request a field that we don't have
	  for(std::vector<std::string>::iterator it= args.begin(); it < args.end(); it++) {
	    if(columns.count(*it) == 0) {
	      //one of the requested columns isn't in the file
	      std::cerr << "ERROR: '" << *it << "' was requested from the dump file, but it doesn't appear to exist. The available data in this frame (tstep = " << last_tstep << ") are: ";
	      for(std::vector<std::string>::iterator jt= avail_columns.begin(); jt < avail_columns.end(); jt++) {
		std::cerr << *jt << " ";
	      }
	      std::cerr << " (" << curfile << ")" << std::endl;
	      return false;
	    }
	  }
	}
      } else {
	//atom data line
	AtomData ad;
	//make sure that everything is zeroed
	//if the user does something silly (like accessing a field they haven't requested), they'll just see a zero
	memset(&ad, 0, sizeof(AtomData));

	if(v.size() != avail_columns.size()) {
	  std::cerr << "ERROR: Mismatch between the number of columns reported and the number of columns read. The LAMMPS header lines indicate " << avail_columns.size() << " columns, but only " << v.size() << " were read. (" << curfile << ")" << std::endl;
	  return false;
	}

	//process the columns that the user wants
	for(std::vector<std::string>::iterator it= args.begin(); it < args.end(); it++) {
	  std::string property= *it;
	  int col= columns[property];
	  if(!updateAtomData(ad, property, v[col])) {
	    std::cerr << "ERROR: LAMMPSReader doesn't know what to do with the property '" << property << "'. This is a shortcoming in LAMMPSReader. ReadFrame will now return false as a precautionary measure, to prevent the return of uninitialised variables." << std::endl;
	    std::cerr << "Aside for the technically minded: To correct this error, modify the private function LAMMPSReader::updateAtomData(AtomData&, std::string, std::string). Add a new Assign line, of the form 'Assign(" << property << ", f)', where f is the function that converts a C-style string to the same type as " << property << " i.e. atoi, atof, etc. (" << curfile << ")" << std::endl;
	    return false;
	  }
	  
	  //check the PBCs
	  //LAMMPS only updates them on reneighbouring steps
	  if((property.compare("x") == 0) && (boundaries[0][0] == 'p') && (ad.x < box_lo[0])) {
	    ad.x+= (box_hi[0] - box_lo[0]);
	  } else if((property.compare("x") == 0) && (boundaries[0][1] == 'p') && (ad.x >= box_hi[0])) {
	    ad.x-= (box_hi[0] - box_lo[0]);
	  } else if((property.compare("y") == 0) && (boundaries[1][0] == 'p') && (ad.y < box_lo[1])) {
	    ad.y+= (box_hi[1] - box_lo[1]);
	  } else if((property.compare("y") == 0) && (boundaries[1][1] == 'p') && (ad.y >= box_hi[1])) {
	    ad.y-= (box_hi[1] - box_lo[1]);
	  } else if((property.compare("z") == 0) && (boundaries[2][0] == 'p') && (ad.z < box_lo[2])) {
	    ad.z+= (box_hi[2] - box_lo[2]);
	  } else if((property.compare("z") == 0) && (boundaries[2][1] == 'p') && (ad.z >= box_hi[2])) {
	    ad.z-= (box_hi[2] - box_lo[2]);
	  } else if((property.compare("xs") == 0) && (boundaries[0][0] == 'p') && (ad.xs <0.0)) {
	    ad.xs+= 1.0;
	  } else if((property.compare("xs") == 0) && (boundaries[0][1] == 'p') && (ad.xs >= 1.0)) {
	    ad.xs-= 1.0;
	  } else if((property.compare("ys") == 0) && (boundaries[1][0] == 'p') && (ad.ys < 0.0)) {
	    ad.ys+= 1.0;
	  } else if((property.compare("ys") == 0) && (boundaries[1][1] == 'p') && (ad.ys >= 1.0)) {
	    ad.ys-= 1.0;
	  } else if((property.compare("zs") == 0) && (boundaries[2][0] == 'p') && (ad.zs < 0.0)) {
	    ad.zs+= 1.0;
	  } else if((property.compare("zs") == 0) && (boundaries[2][1] == 'p') && (ad.zs >= 1.0)) {
	    ad.z-= 1.0;
	  }
	}

	//pass this atom data onto the callback function that the user provided
	c->AtomLine(ad, this);
      }
      line_start= file.tellg();
    }
    //when we hit the end of the file, we've also read a new timestep
    c->EndOfTimestep(this);
    return true;
  }
  
  bool LAMMPSReader::ReadBinaryFrame(const std::vector<std::string>& args, Callback* c) {
    file.read(ubi.buf, sizeof(int64_t));
    //we do a quick check here to make sure that we haven't hit the end of the file
    if(file.fail()) {
      return false;
    }
    last_tstep= static_cast<int>(ubi.i);
    
    file.read(ubi.buf, sizeof(int64_t));
    n_atoms= static_cast<int>(ubi.i);
    
    file.read(ui.buf, sizeof(int));
    if(ui.i) {
      std::cerr << "ERROR: LAMMPSReader does not currently support triclinic boxes." << std::endl;
      return false;
    }
    
    for(int i= 0; i < 3; i++) {
      file.read(ui.buf, sizeof(int));
      if(ui.i == 0) {
	boundaries[i][0]= 'p';
      } else if(ui.i == 1) {
	boundaries[i][0]= 'f';
      } else if(ui.i == 2) {
	boundaries[i][0]= 's';
      } else if(ui.i == 3) {
	boundaries[i][0]= 'm';
      }
      file.read(ui.buf, sizeof(int));
      if(ui.i == 0) {
	boundaries[i][1]= 'p';
      } else if(ui.i == 1) {
	boundaries[i][1]= 'f';
      } else if(ui.i == 2) {
	boundaries[i][1]= 's';
      } else if(ui.i == 3) {
	boundaries[i][1]= 'm';
      }
    }
    
    double box[6];
    for(int i= 0; i < 6; i++) {
      file.read(ud.buf, sizeof(double));
      box[i]= ud.d;
    }
    box_lo[0]= box[0];
    box_lo[1]= box[2];
    box_lo[2]= box[4];
    box_hi[0]= box[1];
    box_hi[1]= box[3];
    box_hi[2]= box[5];
    
    file.read(ui.buf, sizeof(int));
    unsigned int fields_per_atom= ui.i;
    if(fields_per_atom != args.size()) {
      std::cerr << "ERROR: LAMMPSReader was told to expect " << args.size() << " fields per atom from the binary file, but the file reports that there are " << fields_per_atom << ". Remember that when reading binary files, the argument passed to ReadFrame() must specify EVERY field in the dump file." << std::endl;
      return false;
    }
    
    property *fields= new property[fields_per_atom];
    for(unsigned int i= 0; i < args.size(); i++) {
      fields[i]= string_to_property(args[i]);
      if(fields[i] == NULL_PROPERTY) {
	std::cerr << "ERROR: LAMMPSReader doesn't know what to do with the property " << args[i] << std::endl;
	delete[] fields;
	return false;
      }
    }
    
    //Atom data comes in processor blocks!
    //we get the number of processors first
    //then the number of doubles from proc 1
    //then the data from proc 1
    //then the number of doubles from proc 2
    //then the data from proc 2
    //etc.
    file.read(ui.buf, sizeof(int)); //the number of processors used.
    int nprocs= ui.i;
    
    //check that we haven't flagged any errors in the file
    if(file.fail()) {
      std::cerr << "ERROR: LAMMPSReader encountered an error when reading the binary file. This suggests that either your binary file is corrupted, or is of a different format. The LAMMPSReader README file explains the format that it expects to encounter." << std::endl;
      return false;
    }
    
    //if that's fine, handle the start of timestep and box hooks
    c->StartOfTimestep(this);
    c->BoxBounds(boundaries, box_lo, box_hi);
    int atoms_total= 0;
    for(int i= 0; i < nprocs; i++) {
      file.read(ui.buf, sizeof(int)); //buffer size per atom.
      int bufsize= ui.i;
      unsigned int field= 0;
      for(int j= 0; j < bufsize; j++) {
	AtomData ad;
	memset(&ad, 0, sizeof(AtomData));
	file.read(ud.buf, sizeof(double));
	double val= ud.d;
	int ival= static_cast<int>(ud.d);
	
	switch(fields[field]) {
	  case ID:
	    ad.id= ival;
	    break;
	  case TYPE:
	    ad.type= ival;
	    break;
	  case MOL:
	    ad.mol= ival;
	    break;
	  case MASS:
	    ad.mass= val;
	    break;
	  case X:
	    ad.x= val;
	    if(boundaries[0][0] == 'p' && ad.x < box_lo[0]) {
	      ad.x+= (box_hi[0] - box_lo[0]);
	    } else if(boundaries[0][1] == 'p' && ad.x >= box_hi[0]) {
	      ad.x-= (box_hi[0] - box_lo[0]);
	    }
	    break;
	  case Y:
	    ad.y= val;
	    if(boundaries[1][0] == 'p' && ad.y < box_lo[1]) {
	      ad.y+= (box_hi[1] - box_lo[1]);
	    } else if(boundaries[1][1] == 'p' && ad.y >= box_hi[1]) {
	      ad.y-= (box_hi[1] - box_lo[1]);
	    }
	    break;
	  case Z:
	    ad.z= val;
	    if(boundaries[2][0] == 'p' && ad.z < box_lo[2]) {
	      ad.z+= (box_hi[2] - box_lo[2]);
	    } else if(boundaries[2][1] == 'p' && ad.z >= box_hi[2]) {
	      ad.z-= (box_hi[2] - box_lo[2]);
	    }
	    break;
	  case XS:
	    ad.xs= val;
	    if(boundaries[0][0] == 'p' && ad.xs < 0.0) {
	      ad.xs+= 1.0;
	    } else if(boundaries[0][1] == 'p' && ad.xs >= 1.0) {
	      ad.xs-= 1.0;
	    }
	    break;
	  case YS:
	    ad.ys= val;
	    if(boundaries[1][0] == 'p' && ad.ys < 0.0) {
	      ad.ys+= 1.0;
	    } else if(boundaries[1][1] == 'p' && ad.ys >= 1.0) {
	      ad.ys-= 1.0;
	    }
	    break;
	  case ZS:
	    ad.zs= val;
	    if(boundaries[2][0] == 'p' && ad.zs < 0.0) {
	      ad.zs+= 1.0;
	    } else if(boundaries[2][1] == 'p' && ad.zs >= 1.0) {
	      ad.zs-= 1.0;
	    }
	    break;
	  case XU:
	    ad.xu= val;
	    break;
	  case YU:
	    ad.yu= val;
	    break;
	  case ZU:
	    ad.zu= val;
	    break;
	  case XSU:
	    ad.xsu= val;
	    break;
	  case YSU:
	    ad.ysu= val;
	    break;
	  case ZSU:
	    ad.zsu= val;
	    break;
	  case IX:
	    ad.ix= ival;
	    break;
	  case IY:
	    ad.iy= ival;
	    break;
	  case IZ:
	    ad.iz= ival;
	    break;
	  case VX:
	    ad.vx= val;
	    break;
	  case VY:
	    ad.vy= val;
	    break;
	  case VZ:
	    ad.vz= val;
	    break;
	  case FX:
	    ad.fx= val;
	    break;
	  case FY:
	    ad.fy= val;
	    break;
	  case FZ:
	    ad.fz= val;
	    break;
	  case MUX:
	    ad.mux= val;
	    break;
	  case MUY:
	    ad.muy= val;
	    break;
	  case MUZ:
	    ad.muz= val;
	    break;
	  case MU:
	    ad.mu= val;
	    break;
	  case Q:
	    ad.q= val;
	    break;
	  case NULL_PROPERTY:
	    //we can't, but we include it for completeness
	    delete[] fields;
	    return false;
	    break;
	}
	
	if(++field == fields_per_atom) {
	  c->AtomLine(ad, this);
	  memset(&ad, 0, sizeof(AtomData));
	  atoms_total++;
	  field= 0;
	}
      }
    }
    
    if(atoms_total != n_atoms) {
      std::cerr << "Error: total number of atoms provided by the file (" << atoms_total << ") doesn't match the number in the header (" << n_atoms << ")!" << std::endl;
      delete[] fields;
      return false;
    }
    
    //if we made it this far, do the end of timestep hook
    c->EndOfTimestep(this);
    
    delete[] fields;
    return true;
  }
  
  LAMMPSReader::property LAMMPSReader::string_to_property(const std::string& s) {
    #define AddProperty(prop, str) if(s.compare(str) == 0) { return prop; }
    AddProperty(ID, "id");
    AddProperty(TYPE, "type");
    AddProperty(MOL, "mol");
    AddProperty(MASS, "mass");
    AddProperty(X, "x");
    AddProperty(Y, "y");
    AddProperty(Z, "z");
    AddProperty(XS, "xs");
    AddProperty(YS, "ys");
    AddProperty(ZS, "zs");
    AddProperty(XU, "xu");
    AddProperty(YU, "yu");
    AddProperty(ZU, "zu");
    AddProperty(XSU, "xsu");
    AddProperty(YSU, "ysu");
    AddProperty(ZSU, "zsu");
    AddProperty(IX, "ix");
    AddProperty(IY, "iy");
    AddProperty(IZ, "iz");
    AddProperty(VX, "vx");
    AddProperty(VY, "vy");
    AddProperty(VZ, "vz");
    AddProperty(FX, "fx");
    AddProperty(FY, "fy");
    AddProperty(FZ, "fz");
    AddProperty(Q, "q");
    AddProperty(MUX, "mux");
    AddProperty(MUY, "muy");
    AddProperty(MUZ, "muz");
    AddProperty(MU, "mu");
    return NULL_PROPERTY;
  }



  bool LAMMPSReader::updateAtomData(AtomData& ad, const std::string& prop, const std::string& val) {
    #define Assign(tag,conv) if(prop.compare(#tag) == 0) { ad.tag= conv(val.c_str()); return true; }
    Assign(id, atoi)
    Assign(type, atoi)
    Assign(mol, atoi)
    Assign(mass, atof)
    Assign(x, atof)
    Assign(y, atof)
    Assign(z, atof)
    Assign(xs, atof)
    Assign(ys, atof)
    Assign(zs, atof)
    Assign(xu, atof)
    Assign(yu, atof)
    Assign(zu, atof)
    Assign(xsu, atof)
    Assign(ysu, atof)
    Assign(zsu, atof)
    Assign(ix, atoi)
    Assign(iy, atoi)
    Assign(iz, atoi)
    Assign(vx, atof)
    Assign(vy, atof)
    Assign(vz, atof)
    Assign(fx, atof)
    Assign(fy, atof)
    Assign(fz, atof)
    Assign(q, atof)
    Assign(mux, atof)
    Assign(muy, atof)
    Assign(muz, atof)
    return false;
  }
};
