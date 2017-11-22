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

#include "ElementCartesien.h"

using namespace std;

//***********************************************************************

ElementCartesien::ElementCartesien() : m_elementsEnfants(0) {}

//***********************************************************************

ElementCartesien::~ElementCartesien()
{
  for (unsigned int i = 0; i < m_elementsEnfants.size(); i++) {
    delete m_elementsEnfants[i];
  }
  m_elementsEnfants.clear();
}

//***********************************************************************

void ElementCartesien::setVolume(const double &volume)
{
  m_volume = volume;
}

//***********************************************************************

void ElementCartesien::setLCFL(const double &lCFL)
{
  m_lCFL = lCFL;
}

//***********************************************************************

void ElementCartesien::setPos(const double &X, const double &Y, const double &Z)
{
  m_position.setXYZ(X, Y, Z);
}

//***********************************************************************


void ElementCartesien::setPos(const Coord &pos)
{
  m_position = pos;
}

//***********************************************************************

void ElementCartesien::setPosX(const double &X)
{
  m_position.setX(X);
}

//***********************************************************************

void ElementCartesien::setPosY(const double &Y)
{
  m_position.setY(Y);
}

//***********************************************************************

void ElementCartesien::setPosZ(const double &Z)
{
  m_position.setZ(Z);
}

//****************************************************************************
//***************************** Methode AMR **********************************
//****************************************************************************

void ElementCartesien::creerElementEnfant()
{
  m_elementsEnfants.push_back(new ElementCartesien);
}

//****************************************************************************

Element* ElementCartesien::getElementEnfant(const int &numeroEnfant)
{
  return m_elementsEnfants[numeroEnfant];
}

//****************************************************************************

void ElementCartesien::finaliseElementsEnfants()
{
  for (unsigned int i = 0; i < m_elementsEnfants.size(); i++) {
    delete m_elementsEnfants[i];
  }
  m_elementsEnfants.clear();
}
