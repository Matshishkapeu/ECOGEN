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

#ifndef MAILLAGENONSTRUCT_H
#define MAILLAGENONSTRUCT_H

#include <fstream>
#include <vector>
#include <stdint.h>

#include "Maillage.h"
#include "MaillageNonStruct/EnteteElements.h"
#include "../EntreesSorties/IO.h"

class MaillageNonStruct : public Maillage
{
public:
  MaillageNonStruct(const std::string &fichierMaillage);
  ~MaillageNonStruct();

  virtual void attributLimites(std::vector<CondLim*> &condLim);
  virtual int initialiseGeometrie(Cellule ***cellules, BordDeMaille ***bord, bool pretraitementParallele = true, std::string ordreCalcul = "ORDRE1");
  virtual void effetsMaillage(BordDeMaille **face, const int &nombrePhases) const{};
  virtual void repriseCalcul(Cellule **cellules, const int &nombrePhases, const int &nombreTransports, const int &repriseFichier, bool ecritVTK, bool ecritXML = false, bool ecritBinaire = false, bool cree = false);
  virtual std::string quiSuisJe() const { return 0; };

  //Ecriture
  virtual void ecritEntetePiece(std::ofstream &fluxFichier, std::vector<Cellule *> *cellulesLvl, int lvl = 0) const;
  virtual void recupereNoeuds(std::vector<double> &jeuDonnees, int lvl = 0) const;
  virtual void recupereConnectivite(std::vector<double> &jeuDonnees, int lvl = 0) const;
  virtual void recupereOffsets(std::vector<double> &jeuDonnees, int lvl = 0) const;
  virtual void recupereTypeCell(std::vector<double> &jeuDonnees, int lvl = 0) const;
  virtual void recupereDonnees(std::vector<Cellule *> *cellulesLvl, std::vector<double> &jeuDonnees, const int var, int phase, int lvl = 0) const;

  //Pour methode AMR
  virtual void procedureRaffinementInitialisation(std::vector<Cellule *> *cellulesLvl, std::vector<BordDeMaille *> *bordsLvl,
    const std::vector<PhysAdd*> &physAdd, Modele *modele, int &nbMaillesTotalAMR, std::vector<DomaineGeometrique*> &domaines, Cellule **cellules, Eos **eos) { nbMaillesTotalAMR = m_nombreCellulesCalcul; };
  virtual void procedureRaffinement(std::vector<Cellule *> *cellulesLvl, std::vector<BordDeMaille *> *bordsLvl, const int &lvl,
    const std::vector<PhysAdd*> &physAdd, Modele *modele, int &nbMaillesTotalAMR, Cellule **cellules, Eos **eos) {};

private:
  virtual void initialiseGeometrieMonoCPU(Cellule ***cellules, BordDeMaille ***bord, std::string ordreCalcul);
  virtual void initialiseGeometrieParallele(Cellule ***cellules, BordDeMaille ***bord, std::string ordreCalcul);
  void pretraitementFichierMaillageGmsh();
  void lectureGeometrieGmsh(std::vector<ElementNS*>** voisinsNoeuds);
  void lectureGeometrieGmshParallele();
  void lectureElementGmsh(const Coord *TableauNoeuds, std::ifstream &fichierMaillage, ElementNS **element);

  void repriseXML(Cellule **cellules, const int &nombrePhases) const;

  void rechercheElementsArrieres(ElementNS *element, FaceNS *face, BordDeMaille *bord, std::vector<ElementNS *> voisins, Cellule **cellules) const;
  void rechercheElementsAvants(ElementNS *element, FaceNS *face, BordDeMaille *bord, std::vector<ElementNS *> voisins, Cellule **cellules) const;

  std::string m_fichierMaillage;  /*nom du fichier de maillage lu*/
  std::string m_nomMaillage;

  int m_nombreNoeuds;               /*nombre de noeuds definissant le domaine geometrique*/
  int m_nombreNoeudsInternes;       /*nombre de noeuds interne (hors fantomes)*/
  Coord *m_noeuds;                  /*Tableau des coordonnees des noeuds du domaine geometrique*/
  int m_nombreElementsInternes;     /*Nombre d'elements de dimension n de calcul internes */
  int m_nombreElementsFantomes;     /*Nombre d'elements fantomes de dimension n pour calcul parallele */
  int m_nombreElementsCommunicants; /*Nombre reel d'elements communicants*/
  ElementNS **m_elements;           /*Tableau des elements geometriques internes*/
  FaceNS **m_faces;                 /*Tableau des face geometriques*/
  std::vector<CondLim*> m_lim;      /*Tableau des conditions aux limites*/

  int m_nombreFacesInternes;        /*nombre de faces entre deux cellules de calcul*/
  int m_nombreFacesLimites;         /*nombre de faces entre une cellule de calcul et une limite*/
  int m_nombreFacesParallele;       /*nombre de faces entre une cellule de calcul et une ghost cell*/

  int m_nombreCellulesFantomes;     /*nombre de cellules fantomes*/

  int m_nombreElements0D;
  int m_nombreElements1D;
  int m_nombreElements2D;
  int m_nombreElements3D;
  int m_nombreSegments;
  int m_nombreTriangles;
  int m_nombreQuadrangles;
  int m_nombreTetraedres;
  int m_nombrePyramides;
  int m_nombrePoints;
  int m_nombreHexaedres;
};

#endif // MAILLAGENONSTRUCT_H