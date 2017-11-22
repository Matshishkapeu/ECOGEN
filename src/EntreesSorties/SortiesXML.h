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

#ifndef SORTIESXML_H
#define SORTIESXML_H

#include "Sorties.h"

class SortiesXML :  public Sorties
{
public:
  SortiesXML();
  SortiesXML(std::string casTest, std::string run, tinyxml2::XMLElement *element, std::string nomFichier, Entrees *entree);
  virtual ~SortiesXML();

  virtual void prepareSortieSpecifique();
  virtual void ecritSolutionSpecifique(Maillage *maillage, std::vector<Cellule *> *cellulesLvl);

protected:

  std::string creationNomFichierXML(const char* nom, Maillage *maillage=0, int lvl=-1, int proc=-1, int numFichier=-1, std::string nomVariable ="defaut");

  void ecritSolutionXML(Maillage *maillage, std::vector<Cellule *> *cellulesLvl);
  void ecritCollectionXML(Maillage *maillage);
  void ecritDonneesPhysiquesXML(Maillage *maillage, std::vector<Cellule *> *cellulesLvl, std::ofstream &fluxFichier, const int &lvl, bool parallel = false);

  //Dependant du type de maillage
  void ecritMaillageRectilinearXML(Maillage *maillage, std::vector<Cellule *> *cellulesLvl, std::ofstream &fluxFichier, bool parallel = false);
  void ecritMaillageUnstructuredXML(Maillage *maillage, std::vector<Cellule *> *cellulesLvl, std::ofstream &fluxFichier, const int &lvl, bool parallel = false);
  void ecritMaillagePolyDataXML(Maillage *maillage, std::vector<Cellule *> *cellulesLvl, std::ofstream &fluxFichier, const int &lvl, bool parallel = false);
  void ecritFinFichierRectilinearXML(std::ofstream &fluxFichier, bool parallel = false);
  void ecritFinFichierUnstructuredXML(std::ofstream &fluxFichier, bool parallel = false);
  void ecritFinFichierPolyDataXML(std::ofstream &fluxFichier, bool parallel = false);

  //Non utilise
  void ecritFichierParallelXML(Maillage *maillage, std::vector<Cellule *> *cellulesLvl);
};

#endif //SORTIESXML_H