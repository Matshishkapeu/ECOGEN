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

#ifndef DOMAINEGEOMETRIQUE_H
#define DOMAINEGEOMETRIQUE_H

#include <string>
#include <vector>
#include <cmath>
#include "../Maths/Coord.h"
#include "../Modeles/Phase.h"
#include "../libTierces/tinyxml2.h"
#include "../Erreurs.h"
#include "../Outils.h"

class DomaineGeometrique; //Predeclaration de la classe DomaineGeometrique pour pouvoir inclure Cellule.h
#include "../Cellule.h"

class DomaineGeometrique
{
public:
  DomaineGeometrique();
  DomaineGeometrique(std::string nom, std::vector<Phase*> vecPhases, Melange *melange, std::vector<Transport> vecTransports);
  virtual ~DomaineGeometrique();

  virtual bool appartient(Coord &posMaille) const = 0;
  void rempli(Cellule *cellule, const int &nombrePhases, const int &nombreTransports) const;

private:
  std::string m_nom;
  int m_nombrePhases;
  int m_nombreTransports;
  Phase **m_vecPhases;
  Melange *m_melange;
  Transport *m_vecTransports;
};

#endif // DOMAINEGEOMETRIQUE_H

