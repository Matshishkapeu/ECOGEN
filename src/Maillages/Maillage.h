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

#ifndef MAILLAGE_H
#define MAILLAGE_H

#include <ctime>
#include <fstream>
#include <vector>
#include "../libTierces/tinyxml2.h"
#include "../BordDeMaille.h"
#include "../Ordre2/BordDeMailleO2.h"
#include "../Ordre2/CelluleO2.h"
#include "../Ordre2/CelluleO2Ghost.h"
#include "../CondLims/EnteteCondLim.h"
#include "../Maths/Coord.h"
#include "../Parallel.h"
#include "../PhysiqueAdditionnelle/EnteteGPA.h"
#include "../Maths/ObjetGeometrique.h"

class Maillage
{
public:
  Maillage();
  virtual ~Maillage();

  virtual void attributLimites(std::vector<CondLim*> &condLim) = 0;
  virtual int initialiseGeometrie(Cellule ***cellules, BordDeMaille ***bord, bool pretraitementParallele = true, std::string ordreCalcul = "ORDRE1") = 0; //!< renvoi le nombre de dimensions (1,2 ou 3)
  virtual void effetsMaillage(BordDeMaille **face, const int &nombrePhases) const = 0;
  //virtual void ecritSolution(Cellule **cellules, std::vector<Cellule *> *cellulesLvl, const int &nombrePhases, const int &nombreTransports, const int &lvlMax, std::string const &fichier,
  //  bool ecritVTK, std::string variableConstanteCoupe1, std::string variableConstanteCoupe2, const double &valeurCoupe1, const double &valeurCoupe2,
  //  bool coupe1Dde2D = false, bool coupe1Dde3D = false, bool coupe2Dde3D = false, bool ecritXML = false, bool ecritBinaire = false, bool cree = false) const = 0;
  virtual void repriseCalcul(Cellule **cellules, const int &nombrePhases, const int &nombreTransports, const int &repriseFichier, bool ecritVTK, bool ecritXML = false, bool ecritBinaire = false, bool cree = false) { Erreurs::messageErreur("reprise non prevue pour maillage considere"); };
  virtual std::string quiSuisJe() const { Erreurs::messageErreur("quiSuisJe pas prevu pour le maillage demande"); return 0; };

  //Accesseur
  int getGeometrie() const { return m_geometrie; };
  int getNombreCellules() const;
  int getNombreCellulesTotal() const;
  int getNombreFaces() const;
  int getNumFichier() const;
  virtual double getdX() const { return 0; };
  virtual double getdY() const { return 0; };
  virtual double getdZ() const { return 0; };
  TypeM getType() const { return m_type; };
  virtual int getLvlMax() const { return 0; };

  //Ecriture
  void ecritSolutionGnuplot(std::vector<Cellule *> *cellulesLvl, std::ofstream &fluxFichier, ObjetGeometrique *objet = 0) const;

  virtual void ecritEntetePiece(std::ofstream &fluxFichier, std::vector<Cellule *> *cellulesLvl, int lvl = 0) const { Erreurs::messageErreur("ecritEntetePiece non prevu pour maillage considere"); };
  virtual std::string recupereChaineExtent(int rangLocal, bool global = false) const { Erreurs::messageErreur("recupereChaineExtent non prevu pour maillage considere"); return 0; };
  virtual void recupereCoord(std::vector<Cellule *> *cellulesLvl, std::vector<double> &jeuDonnees, Axe axe) const { Erreurs::messageErreur("recupereCoord non prevu pour maillage considere"); };
  virtual void recupereNoeuds(std::vector<double> &jeuDonnees, int lvl = 0) const { Erreurs::messageErreur("recupereNoeuds non prevu pour maillage considere"); };
  virtual void recupereConnectivite(std::vector<double> &jeuDonnees, int lvl = 0) const { Erreurs::messageErreur("recupereConnectivite non prevu pour maillage considere"); };
  virtual void recupereOffsets(std::vector<double> &jeuDonnees, int lvl = 0) const { Erreurs::messageErreur("recupereOffsets non prevu pour maillage considere"); };
  virtual void recupereTypeCell(std::vector<double> &jeuDonnees, int lvl = 0) const { Erreurs::messageErreur("recupereTypeCell non prevu pour maillage considere"); };
  virtual void recupereDonnees(std::vector<Cellule *> *cellulesLvl, std::vector<double> &jeuDonnees, const int var, int phase, int lvl = 0) const { Erreurs::messageErreur("recupereDonnees non prevu pour maillage considere"); };
  
  //Pour methode AMR
  virtual void genereTableauxCellulesBordsLvl(Cellule **cellules, BordDeMaille **bord, std::vector<Cellule *> **cellulesLvl,
    std::vector<BordDeMaille *> **bordsLvl);
  virtual void procedureRaffinementInitialisation(std::vector<Cellule *> *cellulesLvl, std::vector<BordDeMaille *> *bordsLvl,
    const std::vector<PhysAdd*> &physAdd, Modele *modele, int &nbMaillesTotalAMR, std::vector<DomaineGeometrique*> &domaines, Cellule **cellules, Eos **eos) { Erreurs::messageErreur("procedureRaffinementInitialisation pas prevu pour le maillage demande"); };
  virtual void procedureRaffinement(std::vector<Cellule *> *cellulesLvl, std::vector<BordDeMaille *> *bordsLvl, const int &lvl,
    const std::vector<PhysAdd*> &physAdd, Modele *modele, int &nbMaillesTotalAMR, Cellule **cellules, Eos **eos) { Erreurs::messageErreur("procedureRaffinement pas prevu pour le maillage demande"); };

	//Pour parallele
	virtual void initialiseCommunicationsPersistantes(const int nombrePhases, const int nombreTransports, Cellule **cellules, std::string ordreCalcul);
	virtual void communicationsPrimitives(Cellule **cellules, Eos **eos, const int &lvl, Prim type = vecPhases);
	virtual void communicationsPentes(Cellule **cellules, const int &lvl);
	virtual void communicationsScalaire(Cellule **cellules, std::string nomScalaire, const int &lvl);
	virtual void communicationsVecteur(Cellule **cellules, std::string nomVecteur, const int &dim, const int &lvl, int num, int indice);
	virtual void communicationsPhysAdd(const std::vector<PhysAdd*> &physAdd, Cellule **cellules, const int &lvl);
  virtual void communicationsTransports(Cellule **cellules, const int &lvl);
	virtual void finaliseParallele(const int &lvlMax);
  
protected:
  mutable char m_decoupe;
  mutable int m_numFichier;

  int m_geometrie;                  /*indicateur 2D/3D*/
  int m_nombreElements;             /*Nombre d'elements au total (cellules de calculs internes de dimension n + elements limites de dimension n-1 + ghost cells de dimensions n)*/
  int m_nombreFacesTotal;           /*Nombre de faces entre deux cellules ou entre une cellule et une limite*/
  int m_nombreCellulesCalcul;       /*Nombre de cellules de calcul internes au domaine*/
	int m_nombreCellulesTotal;        /*Cellules de calcul internes + cellules fantomes dediees aux communications parallele*/

  std::vector<Cellule *> **m_cellules;
  TypeM m_type;

	int m_nombrePhases;
	int m_nombreTransports;
};
#endif // MAILLAGE_H