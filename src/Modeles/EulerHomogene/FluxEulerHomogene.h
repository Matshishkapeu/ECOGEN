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

#ifndef FLUX_EULER_HOMOGENE_H
#define FLUX_EULER_HOMOGENE_H

#include <iostream>
#include "../Flux.h"

class FluxEulerHomogene;

#include "ModEulerHomogene.h"

class FluxEulerHomogene : public Flux
{
  public:
    /** Default constructor */
    FluxEulerHomogene();
    FluxEulerHomogene(ModEulerHomogene *modele);
    /** Default destructor */
    virtual ~FluxEulerHomogene();

    virtual void afficheFlux() const;
    virtual void ajoutFlux(double coefA, const int &nombrePhases);
    virtual void retireFlux(double coefA, const int &nombrePhases);
    virtual void multiplie(double scalaire, const int &nombrePhases);
    virtual void miseEnTampon(Cellule &cell, const int &nombrePhases);
    virtual void construitCons(Phase **phase, const int &nombrePhases, Melange *melange);
    virtual void construitPrim(Phase **phase, Melange *melange, const int &nombrePhases);
    virtual void miseAZero(const int &nombrePhases);
    virtual void ajoutNonCons(double coefA, const Cellule *cell, const int &nombrePhases){};
    virtual void retireNonCons(double coefA, const Cellule *cell, const int &nombrePhases){};
    virtual void relaxPressions(Cellule *cell, const int &nombrePhases, Prim type = vecPhases) const {};
    virtual void correctionEnergie(Cellule *cell, const int &nombrePhases, Prim type = vecPhases) const{};
    
    virtual void ajoutTuyere1D(Coord &normale, double const &surface, Cellule *cell, const int &nombrePhases);
    virtual void retireTuyere1D(Coord &normale, double const &surface, Cellule *cell, const int &nombrePhases);
    
    // Accesseurs
    virtual Coord getQdm() const;
    virtual double getMasseMel() const;
    virtual double getEnergieMel() const;
    virtual void setCons(const Flux *cons, const int &nombrePhases);

  protected:
    double m_masse;
    Coord m_qdm;
    double m_energ;
    ModEulerHomogene *m_modele;

  private:

    friend class ModEulerHomogene;
    //friend class PAEuler; // A modifier en fonction du besoin, exemple: si on veut faire une classe PAEViscosite, on doit ajouter friend class PAEViscosite.

};

extern FluxEulerHomogene fluxTempEulerHomogene;

#endif // FLUX_EULER_HOMOGENE_H


