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

#include <vector>
#include "DGSphere.h"
#include <iostream>

using namespace std;
using namespace tinyxml2;

//***************************************************************

DGSphere::DGSphere(){}

//***************************************************************
//Constructeur geometrie a partir d une lecture au format XML
//type = sphere
//ex :  <donneesSphere rayon = "0.5">
//        <centre x = "1." y = "0.5" z = "0.5"/>
//      </donneesSphere>
DGSphere::DGSphere(string nom, vector<Phase*> vecPhases, Melange *melange, vector<Transport> vecTransports, XMLElement *element, string nomFichier) :
DomaineGeometrique(nom, vecPhases, melange, vecTransports)
{
  XMLElement *sousElement(element->FirstChildElement("donneesSphere"));
  if (sousElement == NULL) throw ErreurXMLElement("donneesSphere", nomFichier, __FILE__, __LINE__);
  //Recuperation des attributs
  //--------------------------
  XMLError erreur;
  //Rayon de la sphere
  erreur = sousElement->QueryDoubleAttribute("rayon", &m_rayon);
  if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut("rayon", nomFichier, __FILE__, __LINE__);
  //Position centre
  double x(0.), y(0.), z(0.);
  XMLElement *coin(sousElement->FirstChildElement("centre"));
  if (coin == NULL) throw ErreurXMLElement("centre", nomFichier, __FILE__, __LINE__);
  erreur = coin->QueryDoubleAttribute("x", &x);
  erreur = coin->QueryDoubleAttribute("y", &y);
  erreur = coin->QueryDoubleAttribute("z", &z);
  m_posCentre.setXYZ(x, y, z);
}

//***************************************************************

DGSphere::DGSphere(std::string nom, vector<Phase*> vecPhases, Melange *melange, vector<Transport> vecTransports, const Coord &posCentre, const double &rayon) :
DomaineGeometrique(nom, vecPhases, melange, vecTransports), m_posCentre(posCentre), m_rayon(rayon)
{}

//***************************************************************

DGSphere::~DGSphere(){}

//***************************************************************

bool DGSphere::appartient(Coord &posMaille) const
{
  double somme;
  somme = pow(posMaille.getX() - m_posCentre.getX(), 2.) 
        + pow(posMaille.getY() - m_posCentre.getY(), 2.) 
        + pow(posMaille.getZ() - m_posCentre.getZ(), 2.);
  if (somme <= m_rayon*m_rayon) return true;
  return false;
}