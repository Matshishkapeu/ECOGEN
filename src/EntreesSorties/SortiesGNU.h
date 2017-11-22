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

#ifndef SORTIESGNU_H
#define SORTIESGNU_H

#include "Sorties.h"
class SortiesGNU : public Sorties
{
public:
  SortiesGNU();
  SortiesGNU(std::string casTest, std::string run, tinyxml2::XMLElement *element, std::string nomFichier, Entrees *entree);
  virtual ~SortiesGNU();

  virtual void prepareSortieSpecifique() {};  //Rien a faire pour cette sortie
  virtual void ecritSolutionSpecifique(Maillage *maillage, std::vector<Cellule *> *cellulesLvl);

protected:

  void ecritScriptGnuplot(const int &dim);

  std::string creationNomFichierGNU(const char* nom, int lvl = -1, int proc = -1, int numFichier = -1, std::string nomVariable = "defaut") const;
  void ecritureBlocGnuplot(std::ofstream &fluxFichier, int &indice, const int &dim);

  std::string m_nomFichierVisu;
};

#endif //SORTIESGNU_H