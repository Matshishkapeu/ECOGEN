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

#include "SourceAxi.h"

using namespace std;
using namespace tinyxml2;

//***********************************************************************

SourceAxi::SourceAxi(){}

//***********************************************************************
/*!
*  Constructeur source a partir d une lecture au format XML
*  ex : <donneesGravite g="9.81" axe="x" direction="positif"/>
*/
SourceAxi::SourceAxi(XMLElement *element, string nomFichier)
{
  XMLElement *sousElement(element->FirstChildElement("donneesAxi"));
  if (sousElement == NULL) throw ErreurXMLElement("donneesAxi", nomFichier, __FILE__, __LINE__);
  //Recuperation des attributs
  //--------------------------
  //Axe d application
  string axe(sousElement->Attribute("axe"));
  Outils::majuscule(axe);
  if (axe == "YM") { m_axe = Y; }
  else { throw ErreurXMLAttribut("axe", nomFichier, __FILE__, __LINE__); }
}


//***********************************************************************

SourceAxi::~SourceAxi(){}

//***********************************************************************

void SourceAxi::prepareSource(Cellule *cell, const int nombrePhases)
{
  double r, v;
  switch (m_axe) {
  //case X: r = cell->getPosition().getX(); v = cell->getMelange()->getVitesse().getX(); break;
  case Y: r = cell->getPosition().getY(); v = cell->getMelange()->getVitesse().getY(); break;
  //case Z: r = cell->getPosition().getZ(); v = cell->getMelange()->getVitesse().getZ(); break;
  default: Erreurs::messageErreur("nom axe inconnu dans SourceAxi::prepareSource");
  }
  cell->getCons()->prepareSourceAxi(cell, nombrePhases, r, v);
}

//***********************************************************************
