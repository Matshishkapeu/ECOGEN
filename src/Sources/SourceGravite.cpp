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

#include "SourceGravite.h"

using namespace std;
using namespace tinyxml2;

//***********************************************************************

SourceGravite::SourceGravite(){}

//***********************************************************************
/*!
 *  Constructeur source a partir d une lecture au format XML
 *  ex : <donneesGravite g="9.81" axe="x" direction="positif"/>
 */
SourceGravite::SourceGravite(XMLElement *element, string nomFichier)
{
  XMLElement *sousElement(element->FirstChildElement("donneesGravite"));
  if (sousElement == NULL) throw ErreurXMLElement("donneesGravite", nomFichier, __FILE__, __LINE__);
  //Recuperation des attributs
  //--------------------------
  XMLError erreur;
  //Valeur de gravite
  erreur = sousElement->QueryDoubleAttribute("g", &m_g);
  if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut("g", nomFichier, __FILE__, __LINE__);
  //Axe d application
  string axe(sousElement->Attribute("axe"));
  Outils::majuscule(axe);
  if (axe == "X"){ m_axe = X; }
  else if (axe == "Y"){ m_axe = Y; }
  else if (axe == "Z"){ m_axe = Z; }
  else { throw ErreurXMLAttribut("axe", nomFichier, __FILE__, __LINE__); }
  //Direction
  string direction(sousElement->Attribute("direction"));
  Outils::majuscule(direction);
  if (direction == "POSITIF"){ m_sens = 1; }
  else if (direction == "NEGATIF"){ m_sens = -1; }
  else { throw ErreurXMLAttribut("direction", nomFichier, __FILE__, __LINE__); }
}

//***********************************************************************

SourceGravite::~SourceGravite(){}

//***********************************************************************

void SourceGravite::prepareSource(Cellule *cell, const int nombrePhases)
{
  //Force de gravite
  double rho = cell->getMelange()->getDensite();
  Coord Fg(0., 0., 0.); //KS//FP// Penser a changer ca pour eviter de creer et detruire a chaque fois des objets de type Coord (gain en temps de calcul)
  switch (m_axe) {
  case X: Fg.setX(m_sens*rho*m_g); break;
  case Y: Fg.setY(m_sens*rho*m_g); break;
  case Z: Fg.setZ(m_sens*rho*m_g); break;
  default: Erreurs::messageErreur("nom axe inconnu dans SourceGravite::prepareSource");
  }
  cell->getCons()->prepareSourceGravite(cell, nombrePhases, Fg);
}

//***********************************************************************