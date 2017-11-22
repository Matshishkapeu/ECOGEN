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

#include "CondLimAbs.h"

using namespace std;

CondLimAbs CondLimDefaut;

//****************************************************************************

CondLimAbs::CondLimAbs(){}

//****************************************************************************

CondLimAbs::CondLimAbs(const CondLimAbs& Source, const int lvl) : CondLim(Source)
{
  m_lvl = lvl;
}

//****************************************************************************

CondLimAbs::CondLimAbs(int numPhysique) : CondLim(numPhysique)
{}

//****************************************************************************

CondLimAbs::~CondLimAbs() {}

//****************************************************************************

void CondLimAbs::creeLimite(BordDeMaille **face)
{
  *face = new CondLimAbs(*(this));
}

//****************************************************************************

void CondLimAbs::resolRiemannLimite(Cellule &cellGauche, const int & nombrePhases, const double & dxGauche, double & dtMax)
{
  m_mod->resolRiemannInterne(cellGauche, cellGauche, nombrePhases, dxGauche, dxGauche, dtMax);
}

//****************************************************************************

void CondLimAbs::resolRiemannTransportLimite(Cellule &cellGauche, const int & nombreTransports) const
{
	m_mod->resolRiemannTransportInterne(cellGauche, cellGauche, nombreTransports);
}

//****************************************************************************
//******************************Methode AMR***********************************
//****************************************************************************

void CondLimAbs::creerBordEnfant()
{
  m_bordsEnfants.push_back(new CondLimAbs(*this, m_lvl + 1));
}

//***********************************************************************