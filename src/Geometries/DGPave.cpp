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
#include "DGPave.h"
#include <iostream>

using namespace std;
using namespace tinyxml2;

//***************************************************************

DGPave::DGPave(){}

//***************************************************************
//Constructeur geometrie a partir d une lecture au format XML
//type=pave
//ex :	<donneesPave lAxeX="0.3" lAxeY="0.2" lAxeZ="0.4">
//        <posCoinInferieur x = "0.4" y = "0.5" z = "0." />
//      </donneesPave>
DGPave::DGPave(string nom, vector<Phase*> vecPhases, Melange *melange, vector<Transport> vecTransports, XMLElement *element, string nomFichier) :
DomaineGeometrique(nom, vecPhases, melange, vecTransports)
{
  XMLElement *sousElement(element->FirstChildElement("donneesPave"));
  if (sousElement == NULL) throw ErreurXMLElement("donneesPave", nomFichier, __FILE__, __LINE__);
  //Recuperation des attributs
  //--------------------------
  XMLError erreur;
  //Longueur le long de l axe X
  erreur = sousElement->QueryDoubleAttribute("lAxeX", &m_lX);
  if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut("lAxeX", nomFichier, __FILE__, __LINE__);
  //Longueur le long de l axe Y
  erreur = sousElement->QueryDoubleAttribute("lAxeY", &m_lY);
  if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut("lAxeX", nomFichier, __FILE__, __LINE__);
  //Longueur le long de l axe Z
  erreur = sousElement->QueryDoubleAttribute("lAxeZ", &m_lZ);
  if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut("lAxeZ", nomFichier, __FILE__, __LINE__);
  //Position coin inferieur
  double x(0.), y(0.), z(0.);
  XMLElement *coin(sousElement->FirstChildElement("posCoinInferieur"));
  if (coin == NULL) throw ErreurXMLElement("posCoinInferieur", nomFichier, __FILE__, __LINE__);
  erreur = coin->QueryDoubleAttribute("x", &x);
  erreur = coin->QueryDoubleAttribute("y", &y);
  erreur = coin->QueryDoubleAttribute("z", &z);
  m_posXmYmZm.setXYZ(x, y, z);
}

//***************************************************************

DGPave::DGPave(string nom, vector<Phase*> vecPhases, Melange *melange, vector<Transport> vecTransports, const Coord &posXmYmZm, const double &lX, const double &lY, const double &lZ) :
DomaineGeometrique(nom, vecPhases, melange, vecTransports), m_posXmYmZm(posXmYmZm), m_lX(lX), m_lY(lY), m_lZ(lZ)
{}

//***************************************************************

DGPave::~DGPave(){}

//***************************************************************

bool DGPave::appartient(Coord &posMaille) const
{
  if ((posMaille.getX() - m_posXmYmZm.getX())<0) return false;
  if ((posMaille.getX() - m_posXmYmZm.getX())>m_lX) return false;
  if ((posMaille.getY() - m_posXmYmZm.getY())<0) return false;
  if ((posMaille.getY() - m_posXmYmZm.getY())>m_lY) return false;
  if ((posMaille.getZ() - m_posXmYmZm.getZ())<0) return false;
  if ((posMaille.getZ() - m_posXmYmZm.getZ())>m_lZ) return false;
  return true;
}