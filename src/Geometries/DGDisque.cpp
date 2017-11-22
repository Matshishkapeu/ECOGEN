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
#include "DGDisque.h"
#include <iostream>

using namespace std;
using namespace tinyxml2;

//***************************************************************

DGDisque::DGDisque(){}

//***************************************************************
//Constructeur geometrie a partir d une lecture au format XML
//ex : <donneesDisque axe1="x" axe2="y" rayon="0.5">
//        <centre x = "0." y = "0." z = "0." />
//     </donneesDisque>
DGDisque::DGDisque(string nom, vector<Phase*> vecPhases, Melange *melange, vector<Transport> vecTransports, XMLElement *element, string nomFichier) :
DomaineGeometrique(nom, vecPhases, melange, vecTransports)
{
  XMLElement *sousElement(element->FirstChildElement("donneesDisque"));
  if (sousElement == NULL) throw ErreurXMLElement("donneesDisque", nomFichier, __FILE__, __LINE__);
  //Recuperation des attributs
  //--------------------------
  XMLError erreur;
  //Rayon
  erreur = sousElement->QueryDoubleAttribute("rayon", &m_rayon);
  if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut("rayon", nomFichier, __FILE__, __LINE__);
  //Axe1
  string axe(sousElement->Attribute("axe1"));
  Outils::majuscule(axe);
  if (axe == "X"){ m_axe1 = X; }
  else if (axe == "Y"){ m_axe1 = Y; }
  else if (axe == "Z"){ m_axe1 = Z; }
  else { throw ErreurXMLAttribut("axe1", nomFichier, __FILE__, __LINE__); }
  //Axe2
  axe = sousElement->Attribute("axe2");
  Outils::majuscule(axe);
  if (axe == "X"){ m_axe2 = X; }
  else if (axe == "Y"){ m_axe2 = Y; }
  else if (axe == "Z"){ m_axe2 = Z; }
  else { throw ErreurXMLAttribut("axe2", nomFichier, __FILE__, __LINE__); }
  //Centre disque
  double x(0.), y(0.), z(0.);
  XMLElement *centre(sousElement->FirstChildElement("centre"));
  if (centre == NULL) throw ErreurXMLElement("centre", nomFichier, __FILE__, __LINE__);
  erreur = centre->QueryDoubleAttribute("x", &x);
  erreur = centre->QueryDoubleAttribute("y", &y);
  erreur = centre->QueryDoubleAttribute("z", &z);
  m_posCentre.setXYZ(x, y, z);
}

//***************************************************************

DGDisque::DGDisque(string nom, vector<Phase*> vecPhases, Melange *melange, vector<Transport> vecTransports, const Coord &posCentre, const double &rayon, const Axe &axe1, const Axe &axe2) :
DomaineGeometrique(nom, vecPhases, melange, vecTransports), m_posCentre(posCentre), m_rayon(rayon), m_axe1(axe1), m_axe2(axe2)
{}

//***************************************************************

DGDisque::~DGDisque(){}

//***************************************************************

bool DGDisque::appartient(Coord &posMaille) const
{
  double somme(0.);
  vector<Axe> axes;
  axes.push_back(m_axe1);
  axes.push_back(m_axe2);

  for (unsigned int i = 0; i < axes.size(); i++)
  {
    switch (axes[i])
    {
    case X:
      somme += pow(posMaille.getX() - m_posCentre.getX(), 2.); break;
    case Y:
      somme += pow(posMaille.getY() - m_posCentre.getY(), 2.); break;
    case Z:
      somme += pow(posMaille.getZ() - m_posCentre.getZ(), 2.); break;
    }
  }
  
  if (somme <= m_rayon*m_rayon) return true;
  return false;
}