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

#include "FaceCartesien.h"

using namespace std;

//***********************************************************************

FaceCartesien::FaceCartesien(){}

//***********************************************************************

FaceCartesien::~FaceCartesien(){}

//***********************************************************************

void FaceCartesien::setSurface(const double &surface)
{
  m_surface = surface;
}

//***********************************************************************

void FaceCartesien::initialiseAutres(const double &surface, const Coord &normale, const Coord &tangente, const Coord &binormale)
{
  m_surface = surface;
  m_normale = normale;
  m_tangente = tangente;
  m_binormale = binormale;
}

//***********************************************************************

void FaceCartesien::setPos(const double &X, const double &Y, const double &Z)
{
  m_position.setXYZ(X, Y, Z);
}

//***********************************************************************

void FaceCartesien::setNormale(const double &X, const double &Y, const double &Z)
{
  m_normale.setXYZ(X, Y, Z);
}

//***********************************************************************

void FaceCartesien::setTangente(const double &X, const double &Y, const double &Z)
{
  m_tangente.setXYZ(X, Y, Z);
}

//***********************************************************************

void FaceCartesien::setBinormale(const double &X, const double &Y, const double &Z)
{
  m_binormale.setXYZ(X, Y, Z);
}

//****************************************************************************
//***************************** Methode AMR **********************************
//****************************************************************************

Face* FaceCartesien::creerNouvelleFace()
{
  return new FaceCartesien;
}

