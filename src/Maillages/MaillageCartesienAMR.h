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

#ifndef MAILLAGECARTESIENAMR_H
#define MAILLAGECARTESIENAMR_H

#include "MaillageCartesien.h"

class MaillageCartesienAMR : public MaillageCartesien
{
public:
  MaillageCartesienAMR(double lX = 1., int nombreMaillesX = 100, double lY = 1., int nombreMaillesY = 1, double lZ = 1., int nombreMaillesZ = 1,
		int lvlMax = 0, double critereVar = 1.e10, bool varRho = false, bool varP = false, bool varU = false, bool varAlpha = false, double xiSplit = 1., double xiJoin = 1.);
  virtual ~MaillageCartesienAMR();

	virtual void genereTableauxCellulesBordsLvl(Cellule **cellules, BordDeMaille **bord, std::vector<Cellule *> **cellulesLvl,
		std::vector<BordDeMaille *> **bordsLvl);
  virtual void procedureRaffinementInitialisation(std::vector<Cellule *> *cellulesLvl, std::vector<BordDeMaille *> *bordsLvl,
		const std::vector<PhysAdd*> &physAdd, Modele *modele, int &nbMaillesTotalAMR, std::vector<DomaineGeometrique*> &domaines,
		Cellule **cellules, Eos **eos);
  virtual void procedureRaffinement(std::vector<Cellule *> *cellulesLvl, std::vector<BordDeMaille *> *bordsLvl, const int &lvl,
    const std::vector<PhysAdd*> &physAdd, Modele *modele, int &nbMaillesTotalAMR, Cellule **cellules, Eos **eos);
  virtual std::string quiSuisJe() const;

  virtual void ecritEntetePiece(std::ofstream &fluxFichier, std::vector<Cellule *> *cellulesLvl, int lvl = 0) const;
  virtual void recupereNoeuds(std::vector<double> &jeuDonnees, int lvl = 0) const;
  virtual void recupereConnectivite(std::vector<double> &jeuDonnees, int lvl = 0) const;
  virtual void recupereOffsets(std::vector<double> &jeuDonnees, int lvl = 0) const;
  virtual void recupereTypeCell(std::vector<double> &jeuDonnees, int lvl = 0) const;
  virtual void recupereDonnees(std::vector<Cellule *> *cellulesLvl, std::vector<double> &jeuDonnees, const int var, int phase, int lvl = 0) const;

  //Accesseurs
  virtual int getLvlMax() const { return m_lvlMax; };

	//Pour parallele
	virtual void initialiseCommunicationsPersistantes(const int nombrePhases, const int nombreTransports, Cellule **cellules, std::string ordreCalcul);
	virtual void communicationsPrimitives(Cellule **cellules, Eos **eos, const int &lvl, Prim type = vecPhases);
	virtual void communicationsPentes(Cellule **cellules, const int &lvl);
	virtual void communicationsScalaire(Cellule **cellules, std::string nomScalaire, const int &lvl);
	virtual void communicationsVecteur(Cellule **cellules, std::string nomVecteur, const int &dim, const int &lvl, int num, int indice);
	virtual void communicationsPhysAdd(const std::vector<PhysAdd*> &physAdd, Cellule **cellules, const int &lvl);
  virtual void communicationsTransports(Cellule **cellules, const int &lvl);
	virtual void finaliseParallele(const int &lvlMax);

private:
  double m_L;                                //!<Longueur utile pour le calcul du pas de temps diffusif
  int m_lvlMax;                              //!<Niveau maximal sur l arbre AMR (si m_lvlMax = 0, pas d AMR)
	double m_critereVar;                       //!<Valeur du critere a depasser sur la variation d'une variable pour le (de)raffinement (met xi=1.)
	bool m_varRho, m_varP, m_varU, m_varAlpha; //!<Choix sur quelle variation on (de)raffine
	double m_xiSplit, m_xiJoin;                //!<Valeur de xi pour split ou join les mailles
	std::vector<Cellule *> *m_cellulesLvlGhost;//!<Tableau de vecteurs contenant les cellules fantomes, un vecteur par niveau.

};

#endif // MAILLAGECARTESIENAMR_H