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
#include "DGRectangle.h"
#include <iostream>

using namespace std;
using namespace tinyxml2;

//***************************************************************

DGRectangle::DGRectangle(){}

//***************************************************************
//Constructeur geometrie a partir d une lecture au format XML
//type = rectangle
// ex : <donneesRectangle axe1 = "x" axe2 = "y" lAxe1 = "0.3" lAxe2 = "0.2">
//        <posCoinInferieur x = "0.4" y = "0.5" z = "0."/>
//      </donneesRectangle>
DGRectangle::DGRectangle(string nom, vector<Phase*> vecPhases, Melange *melange, vector<Transport> vecTransports, XMLElement *element, string nomFichier) :
DomaineGeometrique(nom, vecPhases, melange, vecTransports)
{
  XMLElement *sousElement(element->FirstChildElement("donneesRectangle"));
  if (sousElement == NULL) throw ErreurXMLElement("donneesRectangle", nomFichier, __FILE__, __LINE__);
  //Recuperation des attributs
  //--------------------------
  XMLError erreur;
  //Longueur le long de l axe 1
  erreur = sousElement->QueryDoubleAttribute("lAxe1", &m_lAxe1);
  if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut("lAxe1", nomFichier, __FILE__, __LINE__);
  //Longueur le long de l axe 2
  erreur = sousElement->QueryDoubleAttribute("lAxe2", &m_lAxe2);
  if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut("lAxe2", nomFichier, __FILE__, __LINE__);
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
  //Position coin inferieur
  double x(0.), y(0.), z(0.);
  XMLElement *coin(sousElement->FirstChildElement("posCoinInferieur"));
  if (coin == NULL) throw ErreurXMLElement("posCoinInferieur", nomFichier, __FILE__, __LINE__);
  erreur = coin->QueryDoubleAttribute("x", &x);
  erreur = coin->QueryDoubleAttribute("y", &y);
  erreur = coin->QueryDoubleAttribute("z", &z);
  m_posBasGauche.setXYZ(x, y, z);
}

//***************************************************************

DGRectangle::DGRectangle(string nom, vector<Phase*> vecPhases, Melange *melange, vector<Transport> vecTransports, const Coord &posBasGauche, const double &lAxe1, const double &lAxe2, const Axe &axe1, const Axe &axe2) :
DomaineGeometrique(nom, vecPhases, melange, vecTransports), m_posBasGauche(posBasGauche), m_lAxe1(lAxe1), m_lAxe2(lAxe2), m_axe1(axe1), m_axe2(axe2)
{}

//***************************************************************

DGRectangle::~DGRectangle(){}

//***************************************************************

bool DGRectangle::appartient(Coord &posMaille) const
{
  double somme(0.);
  vector<Axe> axes;
  axes.push_back(m_axe1);
  axes.push_back(m_axe2);
  vector<double> longueurs;
  longueurs.push_back(m_lAxe1);
  longueurs.push_back(m_lAxe2);

  for (unsigned int i = 0; i < axes.size(); i++)
  {
    switch (axes[i])
    {
    case X:
      if (posMaille.getX() - m_posBasGauche.getX()<0) return false;
      if (posMaille.getX() - m_posBasGauche.getX()>longueurs[i]) return false;
    case Y:
      if (posMaille.getY() - m_posBasGauche.getY()<0) return false;
      if (posMaille.getY() - m_posBasGauche.getY()>longueurs[i]) return false;
    case Z:
      if (posMaille.getZ() - m_posBasGauche.getZ()<0) return false;
      if (posMaille.getZ() - m_posBasGauche.getZ()>longueurs[i]) return false;
    }
  }
  return true;
}