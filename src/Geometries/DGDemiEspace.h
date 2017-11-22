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

class DGDemiEspace : public DomaineGeometrique
{
  public:
    DGDemiEspace();
    DGDemiEspace(std::string nom, std::vector<Phase*> vecPhases, Melange *melange, std::vector<Transport> vecTransports, tinyxml2::XMLElement *element, std::string nomFichier="Fichier Inconnu");
    DGDemiEspace(std::string nom, std::vector<Phase*> vecPhases, Melange *melange, std::vector<Transport> vecTransports, const double &position, const Axe &axe, const int &sens);
    virtual ~DGDemiEspace();

    virtual bool appartient(Coord &posMaille) const;

  private:
    double m_position;
    Axe m_axe; 
    int m_sens; //positif pour le sens positif, negatif pour le sens negatif
};

