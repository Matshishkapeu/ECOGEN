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

#ifndef FLUX_H
#define FLUX_H

class Flux; //Predeclaration de la classe Flux pour pouvoir inclure Cellule.h

#include "Phase.h"
#include "../Cellule.h"
#include "../Outils.h"

class Flux
{
  public:
    /** Default constructor */
    Flux();
    /** Default destructor */
    virtual ~Flux();

    virtual void afficheFlux() const { Erreurs::messageErreur("afficheFlux non implemente pour modele utilise"); };
    virtual void ajoutFlux(double coefA, const int &nombrePhases){ Erreurs::messageErreur("ajoutFlux non implemente pour modele utilise"); };
    virtual void retireFlux(double coefA, const int &nombrePhases){ Erreurs::messageErreur("retireFlux non implemente pour modele utilise"); };
    virtual void multiplie(double scalaire, const int &nombrePhases){ Erreurs::messageErreur("multiplie non implemente pour modele utilise"); };
    virtual void miseEnTampon(Cellule &cell, const int &nombrePhases){ Erreurs::messageErreur("miseEnTampon non implemente pour modele utilise"); };
    virtual void construitCons(Phase **phases, const int &nombrePhases, Melange *melange) { Erreurs::messageErreur("construitCons non implemente pour modele utilise"); };
    virtual void construitPrim(Phase **phases, Melange *melange, const int &nombrePhases) { Erreurs::messageErreur("construitPrim non implemente pour modele utilise"); };
    virtual void miseAZero(const int &nombrePhases){ Erreurs::messageErreur("miseAZero non implemente pour modele utilise"); };
    virtual void miseAZeroFluxTemp(const int &nombrePhases) { Erreurs::messageErreur("miseAZero non implemente pour modele utilise"); };
    virtual void ajoutNonCons(double coefA, const Cellule *cell, const int &nombrePhases){ Erreurs::messageErreur("ajoutNonCons non implemente pour modele utilise"); };
    virtual void retireNonCons(double coefA, const Cellule *cell, const int &nombrePhases){ Erreurs::messageErreur("retireNonCons non implemente pour modele utilise"); };
    virtual void relaxPressions(Cellule *cell, const int &nombrePhases, Prim type = vecPhases) const { Erreurs::messageErreur("relaxPressions non implemente pour modele utilise"); };
    virtual void relaxPTMu(Cellule *cell, const int &nombrePhases, Prim type = vecPhases) const { Erreurs::messageErreur("relaxPTMu non implemente pour modele utilise"); };
    virtual double calculTsat(const double &pression, const int &nombrePhases, double *dTsat = 0) const { Erreurs::messageErreur("calculTsat non implemente pour modele utilise"); return 0.; }; //FP//TODO// Methode au mauvais endroit -> temporaire pour aller vite!
    virtual void correctionEnergie(Cellule *cell, const int &nombrePhases, Prim type = vecPhases) const { Erreurs::messageErreur("correctionEnergie non implemente pour modele utilise"); };

    virtual void prepareSourceAxi(Cellule *cell, const int &nombrePhases, double &r, double &v) { Erreurs::messageErreur("axi non implemente pour modele utilise"); };
    virtual void prepareSourceGravite(Cellule *cell, const int &nombrePhases, Coord &Fg){ Erreurs::messageErreur("gravite non implemente pour modele utilise"); };
    virtual void integreTermeSource(Cellule *cell, const double &dt, const int &nombrePhases){ Erreurs::messageErreur("integreTermeSource non implemente pour modele utilise"); };

    virtual void ajoutTuyere1D(Coord &normale, double const &surface, Cellule *cell, const int &nombrePhases){ Erreurs::messageErreur("ajoutTuyere1D non implemente pour modele utilise"); };
    virtual void retireTuyere1D(Coord &normale, double const &surface, Cellule *cell, const int &nombrePhases){ Erreurs::messageErreur("retireTuyere1D non implemente pour modele utilise"); };

    // Accesseurs
    virtual double getAlpha(const int &numPhase) const { return 0.; };
    virtual double getMasse(const int &numPhase) const { return 0.; };
    virtual double getEnergie(const int &numPhase) const { return 0.; };
    virtual Coord getQdm() const { return 0.; };
    virtual double getMasseMel() const { return 0.; };
    virtual double getEnergieMel() const { return 0.; };
    virtual void setCons(const Flux *cons, const int &nombrePhases) { Erreurs::messageErreur("setCons non implemente pour modele utilise"); };

  protected:
    double  m_sM;
  private:
};

#endif // FLUX_H
