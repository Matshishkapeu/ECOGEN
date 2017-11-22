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

#ifndef SORTIESCOUPEGNU_H
#define SORTIESCOUPEGNU_H

#include "SortiesGNU.h"
#include "../Maths/OGDroite.h"
#include "../Maths/OGPlan.h"

class SortiesCoupeGNU : public SortiesGNU
{
public:
  SortiesCoupeGNU();
  SortiesCoupeGNU(std::string casTest, std::string run, tinyxml2::XMLElement *element, std::string nomFichier, TypeOG type, Entrees *entree);
  virtual ~SortiesCoupeGNU();

  virtual void ecritSolutionSpecifique(Maillage *maillage, std::vector<Cellule *> *cellulesLvl);

  virtual void prepareSortiesInfos() {}; //Aucune infos a ecrire
  virtual void ecritInfos() {};

private:
  ObjetGeometrique *m_objet; //droite ou plan de coupe
};

#endif //SORTIESCOUPEGNU_H