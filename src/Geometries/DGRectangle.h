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

#pragma once
#include "DomaineGeometrique.h"

class DGRectangle :
  public DomaineGeometrique
{
public:
  DGRectangle();
  DGRectangle(std::string nom, std::vector<Phase*> vecPhases, Melange *melange, std::vector<Transport> vecTransports, tinyxml2::XMLElement *element, std::string nomFichier = "Fichier Inconnu");
  DGRectangle(std::string nom, std::vector<Phase*> vecPhases, Melange *melange, std::vector<Transport> vecTransports, const Coord &posBasGauche, const double &lAxe1, const double &lAxe2, const Axe &axe1 = X, const Axe &axe2 = Y);
  virtual ~DGRectangle();

  virtual bool appartient(Coord &posMaille) const;
private:
  Coord m_posBasGauche;
  Axe m_axe1, m_axe2;
  double m_lAxe1, m_lAxe2;
};

