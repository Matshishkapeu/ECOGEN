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

#ifndef MAILLAGECARTESIEN_H
#define MAILLAGECARTESIEN_H

#include "Maillage.h"
#include "ElementCartesien.h"
#include "FaceCartesien.h"

class MaillageCartesien : public Maillage
{
public:
  MaillageCartesien(double lX = 1., int nombreMaillesX = 100, double lY = 1., int nombreMaillesY = 1, double lZ = 1., int nombreMaillesZ = 1);
  virtual ~MaillageCartesien();

  virtual void attributLimites(std::vector<CondLim*> &condLim);
  void recupereIJK(const int &indice, int &i, int &j, int &k) const;
  void construitIGlobal(const int &i, const int &j, const int &k, int &indice) const;
  virtual int initialiseGeometrie(Cellule ***cellules, BordDeMaille ***bord, bool pretraitementParallele, std::string ordreCalcul);
  virtual void initialiseGeometrieMonoCPU(Cellule ***cellules, BordDeMaille ***bord, std::string ordreCalcul);
  virtual void initialiseGeometrieParallele(Cellule ***cellules, BordDeMaille ***bord, std::string ordreCalcul);
  virtual void effetsMaillage(BordDeMaille **face, const int &nombrePhases) const {};
  void decoupageParallele();
  virtual std::string quiSuisJe() const;

  //Ecriture
  virtual std::string recupereChaineExtent(int rangLocal, bool global = false) const;
  virtual void recupereCoord(std::vector<Cellule *> *cellulesLvl, std::vector<double> &jeuDonnees, Axe axe) const;
  virtual void recupereDonnees(std::vector<Cellule *> *cellulesLvl, std::vector<double> &jeuDonnees, const int var, int phase, int lvl = 0) const;

  //Accesseur
  virtual double getdX() const;
  virtual double getdY() const;
  virtual double getdZ() const;
  double getVolume() const; //Temporaire car proprietes des mailles normalement...
  double getLCFL() const;

  //Pour methode AMR
  virtual void procedureRaffinementInitialisation(std::vector<Cellule *> *cellulesLvl, std::vector<BordDeMaille *> *bordsLvl,
    const std::vector<PhysAdd*> &physAdd, Modele *modele, int &nbMaillesTotalAMR, std::vector<DomaineGeometrique*> &domaines, Cellule **cellules, Eos **eos) { nbMaillesTotalAMR = m_nombreCellulesCalcul; };
  virtual void procedureRaffinement(std::vector<Cellule *> *cellulesLvl, std::vector<BordDeMaille *> *bordsLvl, const int &lvl,
    const std::vector<PhysAdd*> &physAdd, Modele *modele, int &nbMaillesTotalAMR, Cellule **cellules, Eos **eos) {};

protected:
  ElementCartesien *m_elements;
  FaceCartesien *m_faces;

  double m_lX;
  double m_lY;
  double m_lZ;
  Coord m_origine;
  int m_nombreMaillesX;
  int m_nombreMaillesY;
  int m_nombreMaillesZ;
  int m_nombreMaillesXGlobal;
  int m_nombreMaillesYGlobal;
  int m_nombreMaillesZGlobal;
  int m_numNoeudDeb;  //<! Numero du premier noeud du CPU dans l ensemble du domaine (necessaire pour ecriture XML //)
  int m_numNoeudFin;  //<! Numero du dernier noeud du CPU dans l ensemble du domaine (necessaire pour ecriture XML //)
  double m_dX;
  double m_dY;
  double m_dZ;
  double m_L;         // Pour methode AMR : longueur utile pour le calcul du pas de temps diffusif
  double m_volume;
  double m_lCFL;

  CondLim *m_limXm;
  CondLim *m_limXp;
  CondLim *m_limYm;
  CondLim *m_limYp;
  CondLim *m_limZm;
  CondLim *m_limZp;
};

#endif // MAILLAGECARTESIEN_H