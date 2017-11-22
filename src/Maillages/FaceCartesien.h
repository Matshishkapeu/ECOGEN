//  
//       ,---.     ,--,    .---.     ,--,    ,---.    .-. .-. 
//       | .-'   .' .')   / .-. )  .' .'     | .-'    |  \| | 
//       | `-.   |  |(_)  | | |(_) |  |  __  | `-.    |   | | 
//       | .-'   \  \     | | | |  \  \ ( _) | .-'    | |\  | 
//       |  `--.  \  `-.  \ `-' /   \  `-) ) |  `--.  | | |)| 
//       /( __.'   \____\  )---'    )\____/  /( __.'  /(  (_) 
//      (__)              (_)      (__)     (__)     (__)     
//
//  This file is part of ECOGEN.
//
//  ECOGEN is the legal property of its developers, whose names 
//  are listed in the copyright file included with this source 
//  distribution.
//
//  ECOGEN is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published 
//  by the Free Software Foundation, either version 3 of the License, 
//  or (at your option) any later version.
//  
//  ECOGEN is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//  GNU General Public License for more details.
//  
//  You should have received a copy of the GNU General Public License
//  along with ECOGEN (file LICENSE).  
//  If not, see <http://www.gnu.org/licenses/>.

#ifndef FACECARTESIEN_H
#define FACECARTESIEN_H

#include "Face.h"
class FaceCartesien : public Face
{
public:
  FaceCartesien();
  virtual ~FaceCartesien();

  virtual void setSurface(const double &surface);
  virtual void initialiseAutres(const double &surface, const Coord &normale, const Coord &tangente, const Coord &binormale);
  virtual void setPos(const double &X, const double &Y, const double &Z);
  virtual void setNormale(const double &X, const double &Y, const double &Z);
  virtual void setTangente(const double &X, const double &Y, const double &Z);
  virtual void setBinormale(const double &X, const double &Y, const double &Z);

  //Pour methode AMR
  virtual Face* creerNouvelleFace();

protected:
};

#endif // FACECARTESIEN_H