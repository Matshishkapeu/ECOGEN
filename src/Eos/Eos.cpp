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

#include <iostream>
#include <string>
#include "Eos.h"

using namespace std;
using namespace tinyxml2;

//***********************************************************************

Eos::Eos(){}

//***********************************************************************

Eos::Eos(int &numero) :
m_numero(numero), m_mu(-1.), m_lambda(-1.)
{
  numero++;
}

//***********************************************************************

Eos::Eos(string nom, int &numero) :
m_nom(nom), m_numero(numero), m_mu(-1.), m_lambda(-1.)
{
  numero++;
}

//***********************************************************************

Eos::~Eos(){}

//***********************************************************************

void Eos::lectureParamPhysiques(XMLNode *element, string nomFichier)
{
  XMLError erreur;

  XMLElement *sousElement(element->FirstChildElement("paramPhysiques"));
  if (sousElement != NULL) {
    //Recuperation des donnees
    erreur = sousElement->QueryDoubleAttribute("mu", &m_mu);
    if (erreur != XML_NO_ERROR) m_mu = -1.;
    erreur = sousElement->QueryDoubleAttribute("lambda", &m_lambda);
    if (erreur != XML_NO_ERROR) m_lambda = -1.;
  }
}

//***********************************************************************

void Eos::visualise() const
{
    cout << "Materiau : " << m_nom << endl;
}

//***********************************************************************

string Eos::retourneNom() const
{
  return m_nom;
}

//***********************************************************************

int Eos::getNumero() const
{
  return m_numero;
}

//***********************************************************************

double Eos::calculEnthalpieTotale(const double &densite, const double &pression, const double &vitesse) const
{
  return this->calculEnergie(densite, pression) + pression / densite + 0.5*vitesse*vitesse;
}

//***********************************************************************

double Eos::getMu() const { return m_mu; }

//***********************************************************************

double Eos::getLambda() const { return m_lambda; }

//***********************************************************************