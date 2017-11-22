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

#include "DGDemiEspace.h"

using namespace std;
using namespace tinyxml2;

//***************************************************************

DGDemiEspace::DGDemiEspace(){}

//***************************************************************
/*!
 *  Constructeur geometrie a partir d une lecture au format XML
 *  ex : <donneesDemiEspace axe="x" origine="0.5" direction="positif"/>
 */
DGDemiEspace::DGDemiEspace(string nom, vector<Phase*> vecPhases, Melange *melange, vector<Transport> vecTransports, XMLElement *element, string nomFichier) :
  DomaineGeometrique(nom, vecPhases, melange, vecTransports)
{
  XMLElement *sousElement(element->FirstChildElement("donneesDemiEspace"));
  if (sousElement == NULL) throw ErreurXMLElement("donneesDemiEspace", nomFichier, __FILE__, __LINE__);
  //Recuperation des attributs
  //--------------------------
  XMLError erreur;
  //Origine
  erreur = sousElement->QueryDoubleAttribute("origine", &m_position);
  if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut("origine", nomFichier, __FILE__, __LINE__);
  //Axe
  string axe(sousElement->Attribute("axe"));
  Outils::majuscule(axe);
  if      (axe == "X"){ m_axe = X; }
  else if (axe == "Y"){ m_axe = Y; }
  else if (axe == "Z"){ m_axe = Z; }
  else { throw ErreurXMLAttribut("axe", nomFichier, __FILE__, __LINE__); }
  //Direction
  string direction(sousElement->Attribute("direction"));
  Outils::majuscule(direction);
  if      (direction == "POSITIF"){ m_sens = 1; }
  else if (direction == "NEGATIF"){ m_sens = -1; }
  else { throw ErreurXMLAttribut("direction", nomFichier, __FILE__, __LINE__); }
}

//***************************************************************

DGDemiEspace::DGDemiEspace(string nom, vector<Phase*> vecPhases, Melange *melange, vector<Transport> vecTransports, const double &position, const Axe &axe, const int &sens) :
DomaineGeometrique(nom, vecPhases, melange, vecTransports),m_position(position), m_axe(axe), m_sens(sens)
{}

//***************************************************************

DGDemiEspace::~DGDemiEspace(){}

//***************************************************************

bool DGDemiEspace::appartient(Coord &posMaille) const
{
  bool result(false);
  switch (m_axe)
  {
  case X:
    if (m_sens >= 0){if (posMaille.getX() >= m_position) result = true;}
    else{ if (posMaille.getX() <= m_position) result = true; }
    break;
  case Y:
    if (m_sens >= 0){ if (posMaille.getY() >= m_position) result = true; }
    else{ if (posMaille.getY() <= m_position) result = true; }
    break;
  case Z:
    if (m_sens >= 0){ if (posMaille.getZ() >= m_position) result = true; }
    else{ if (posMaille.getZ() <= m_position) result = true; }
    break;
  default:
    result = false; //Si appartient pas renvoie false
    break;
  }
  return result;
}

//***************************************************************