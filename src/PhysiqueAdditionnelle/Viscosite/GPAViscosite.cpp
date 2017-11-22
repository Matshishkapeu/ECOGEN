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

#include "GPAViscosite.h"
#include <iostream>

using namespace std;

//***********************************************************************

GPAViscosite::GPAViscosite(){}

//***********************************************************************

GPAViscosite::GPAViscosite(PhysAdd* physAdd) : GrandeursPhysAdd(physAdd), m_gradU(0.), m_gradV(0.), m_gradW(0.)
{}

//***********************************************************************

GPAViscosite::~GPAViscosite(){}

//***********************************************************************

void GPAViscosite::calculGrandeurs(Cellule* cell)
{
  m_gradU = cell->calculGradient("u", 0); //0 : On prend la premiere phase
  m_gradV = cell->calculGradient("v", 0);
  m_gradW = cell->calculGradient("w", 0);
}

//***********************************************************************

void GPAViscosite::setGrad(const Coord &grad, int num)
{
  switch (num) {
  case 1: m_gradU = grad; break;
  case 2: m_gradV = grad; break;
  case 3: m_gradW = grad; break;
  default: Erreurs::messageErreur("Erreur dans GPAViscosite::setGrad valeur de num non definie"); break;
  }
}

//***********************************************************************

Coord GPAViscosite::getGrad(int num) const
{
  switch (num) {
  case 1: return m_gradU; break;
  case 2: return m_gradV; break;
  case 3: return m_gradW; break;
  default: Erreurs::messageErreur("Erreur dans GPAViscosite::getGrad valeur de num non definie"); break;
  }
  return 0;
}

//***********************************************************************