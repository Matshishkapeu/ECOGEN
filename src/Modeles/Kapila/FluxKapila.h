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

#ifndef FLUXKAPILA_H
#define FLUXKAPILA_H

#include "../Flux.h"
#include <iostream>

class FluxKapila : public Flux
{
  public:
    /** Default constructor */
    FluxKapila();
    FluxKapila(const int &nombrePhases);
    /** Default destructor */
    virtual ~FluxKapila();

    virtual void afficheFlux() const;
    virtual void ajoutFlux(double coefA, const int &nombrePhases);
    virtual void retireFlux(double coefA, const int &nombrePhases);
    virtual void multiplie(double scalaire, const int &nombrePhases);
    virtual void miseEnTampon(Cellule &cell, const int &nombrePhases);
    virtual void construitCons(Phase **phases, const int &nombrePhases, Melange *melange);
    virtual void construitPrim(Phase **phases, Melange *melange, const int &nombrePhases);
    virtual void miseAZero(const int &nombrePhases);
    virtual void miseAZeroFluxTemp(const int &nombrePhases);
    virtual void ajoutNonCons(double coefA, const Cellule *cell, const int &nombrePhases);
    virtual void retireNonCons(double coefA, const Cellule *cell, const int &nombrePhases);
    virtual void relaxPressions(Cellule *cell, const int &nombrePhases, Prim type = vecPhases) const;
    virtual void relaxPTMu(Cellule *cell, const int &nombrePhases, Prim type = vecPhases) const;
    virtual double calculTsat(const double &pression, const int &nombrePhases, double *dTsat = 0) const;
    virtual void correctionEnergie(Cellule *cell, const int &nombrePhases, Prim type = vecPhases) const;

    virtual void prepareSourceAxi(Cellule *cell, const int &nombrePhases, double &r, double &v);
    virtual void prepareSourceGravite(Cellule *cell, const int &nombrePhases, Coord &Fg);
    virtual void integreTermeSource(Cellule *cell, const double &dt, const int &nombrePhases);

    virtual void ajoutTuyere1D(Coord &normale, double const &surface, Cellule *cell, const int &nombrePhases){};
    virtual void retireTuyere1D(Coord &normale, double const &surface, Cellule *cell, const int &nombrePhases){};

    // Accesseurs
    virtual double getAlpha(const int &numPhase) const;
    virtual double getMasse(const int &numPhase) const;
    virtual double getEnergie(const int &numPhase) const;
    virtual Coord getQdm() const;
    virtual double getEnergieMel() const;
    virtual void setCons(const Flux *cons, const int &nombrePhases);

protected:
    double *m_alpha; //!<taille : nombre de phases
    double *m_masse; //!<taille : nombre de phases
    double *m_energ; //!<taille : nombre de phases
    Coord m_qdm;
    double m_energMelange;

  private:

    friend class ModKapila;
    // A modifier en fonction du besoin, exemple: si on veut faire une classe PAKViscosite, on doit ajouter friend class PAKViscosite.
    friend class PAKCapillarite;
    friend class PAKViscosite;
    friend class PAKConductivite;

};

extern FluxKapila *fluxTempKapila;
extern FluxKapila *sourceConsKap;

#endif // FLUXKAPILA_H
