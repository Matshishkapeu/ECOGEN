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

#include "Transport.h"
#include <fstream>

using namespace std;
using namespace tinyxml2;

Transport* fluxTempTransport;

//***********************************************************************

Transport::Transport() : m_valeur(0.)
{}

//***********************************************************************

Transport::~Transport(){}

//***********************************************************************

void Transport::setValeur(double valeur)
{
  m_valeur = valeur;
}

//***********************************************************************

double Transport::getValeur() const
{
  return m_valeur;
}

//***********************************************************************

void Transport::resolRiemann(double transportGauche, double transportDroite, double sM)
{
	if (sM > 0.) { m_valeur = transportGauche*sM; }
	else { m_valeur = transportDroite*sM; }
}

//***********************************************************************

void Transport::resolRiemannMur()
{
	m_valeur = 0.;
}

//***********************************************************************

void Transport::resolRiemannInj(double transportGauche, double sM, double valeurTransport)
{
  if (sM > 0.) { m_valeur = transportGauche*sM; }
  else { m_valeur = valeurTransport*sM; }
}

//***********************************************************************

void Transport::resolRiemannRes(double transportGauche, double sM, double valeurTransport)
{
  if (sM > 0.) { m_valeur = transportGauche*sM; }
  else { m_valeur = valeurTransport*sM; }
}

//***********************************************************************

void Transport::resolRiemannSortie(double transportGauche, double sM, double valeurTransport)
{
	if (sM > 0.) { m_valeur = transportGauche*sM; }
  else { m_valeur = valeurTransport*sM; }
}

//***********************************************************************

void Transport::ajoutFlux(double coefA, const int num)
{
  m_valeur += coefA*fluxTempTransport[num].m_valeur;
}

//***********************************************************************

void Transport::retireFlux(double coefA, const int num)
{
  m_valeur -= coefA*fluxTempTransport[num].m_valeur;
}

//***********************************************************************

void Transport::ajoutNonCons(double coefA, double transport, const double sM)
{
  m_valeur += -coefA*transport*sM;
}

//***********************************************************************

void Transport::retireNonCons(double coefA, double transport, const double sM)
{
  m_valeur -= -coefA*transport*sM;
}

//***********************************************************************

void Transport::multiplie(double scalaire) 
{
  m_valeur *= scalaire;
}

//***********************************************************************

void Transport::ajoute(double scalaire)
{
  m_valeur += scalaire;
}

//***************************************************************************

void Transport::changeSigne()
{
  m_valeur = -m_valeur;
}

//***************************************************************************

void Transport::calculPentesTransport(const double valeurGauche, const double valeurDroite, const double &distance)
{
  m_valeur = (valeurDroite - valeurGauche) / distance;
}

//***********************************************************************

void Transport::extrapole(const double &pente, const double &distance)
{
  m_valeur += pente * distance;
}