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

#include "CondLimMur.h"

using namespace std;

//****************************************************************************

CondLimMur::CondLimMur(){}

//****************************************************************************

CondLimMur::CondLimMur(const CondLimMur& Source, const int lvl) : CondLim(Source)
{
  m_lvl = lvl;
}

//****************************************************************************

CondLimMur::CondLimMur(int numPhysique) : CondLim(numPhysique)
{}

//****************************************************************************

CondLimMur::~CondLimMur(){}

//****************************************************************************

void CondLimMur::creeLimite(BordDeMaille **face)
{
  *face = new CondLimMur(*(this));
}

//****************************************************************************

void CondLimMur::resolRiemannLimite(Cellule &cellGauche, const int & nombrePhases, const double & dxGauche, double & dtMax)
{
  m_mod->resolRiemannMur(cellGauche, nombrePhases, dxGauche, dtMax);
}

//****************************************************************************

void CondLimMur::resolRiemannTransportLimite(Cellule &cellGauche, const int & nombreTransports) const
{
	m_mod->resolRiemannTransportMur(nombreTransports);
}

//****************************************************************************
//******************************Methode AMR***********************************
//****************************************************************************

void CondLimMur::creerBordEnfant()
{
  m_bordsEnfants.push_back(new CondLimMur(*this, m_lvl + 1));
}

//****************************************************************************