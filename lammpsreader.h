/*
    lammpsreader.h
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

#ifndef LAMMPSREADER_H
#define LAMMPSREADER_H

#include <fstream>
#include <string>
#include <vector>

namespace LAMMPSReader {
  
  std::vector<std::string> explode(std::string);
  
  struct AtomData {
    int id;
    int type;
    int mol;
    double mass;
    double x, y, z;
    double xs, ys, zs;
    double xu, yu, zu;
    double xsu, ysu, zsu;
    int ix, iy, iz;
    double vx, vy, vz;
    double fx, fy, fz;
    double mux, muy, muz;
    double mu;
    double q;
  };

  class LAMMPSReader;

  class Callback {
  public:
    virtual void AtomLine(const AtomData&, LAMMPSReader*) {};
    virtual void BoxBounds(char[3][2], double[3], double[3]) {};
    virtual void StartOfTimestep(LAMMPSReader*) {};
    virtual void EndOfTimestep(LAMMPSReader*) {};
  };
  
  class LAMMPSReader {
  public:
    char boundaries[3][2];

    double box_lo[3];
    double box_hi[3];

    int last_tstep;
    int n_atoms;
    
    LAMMPSReader();
    ~LAMMPSReader();

    bool open(const std::string&, bool bin= false);
    void close();

    bool ReadFrame(const std::string&, Callback *c);
  private:
    bool binary;
    std::ifstream file;
    std::string curfile;
    bool updateAtomData(AtomData&, const std::string&, const std::string&);
    bool ReadBinaryFrame(const std::vector<std::string>&, Callback*);
    
    //because enums are int-based and count from 0, PROPERTIES gives the number of other properties
    enum property {ID, TYPE, MOL, MASS, X, Y, Z, XS, YS, ZS, XU, YU, ZU,
    XSU, YSU, ZSU, IX, IY, IZ, VX, VY, VZ, FX, FY, FZ, Q, MUX, MUY, MUZ, 
    MU, NULL_PROPERTY};
    property string_to_property(const std::string& s);
    
    //the following unions are used for reading binary files
    
    union bigint_ {
      char buf[sizeof(int64_t)];
      int64_t i;
    } ubi;

    union int_ {
      char buf[sizeof(int)];
      int i;
    } ui;

    union double_ {
      char buf[sizeof(double)];
      double d;
    } ud;
  };
}

#endif
