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

#ifndef MODEULERHOMOGENE_H
#define MODEULERHOMOGENE_H

#include "../Modele.h"
#include "../../Cellule.h"
#include "MelEulerHomogene.h"

class ModEulerHomogene;

#include "FluxEulerHomogene.h"

class ModEulerHomogene : public Modele
{
  public:
    /** Default constructor */
    ModEulerHomogene(const int &nombreTransports, const int liquide = 1, const int vapeur =2);
    /** Default destructor */
    virtual ~ModEulerHomogene();

    virtual void alloueCons(Flux **cons, const int &nombrePhases);
    virtual void allouePhase(Phase **phase);
    virtual void alloueMelange(Melange **melange);
	  virtual void alloueEos(Cellule &cell, const int &nombrePhases) {};

  	// Methodes pour hydrodynamique
    virtual void resolRiemannInterne(Cellule &cellGauche, Cellule &cellDroite, const int &nombrePhases, const double &dxGauche, const double &dxDroite, double &dtMax) const; // Riemann entre deux maille de calcul
    //virtual void resolRiemannMur(Cellule &cellGauche, const int &nombrePhases, const double &dxGauche, double &dtMax) const; // Riemann entre maille gauche et mur
    //virtual void resolRiemannInj(Cellule &cellGauche, const int &nombrePhases, const double &dxGauche, double &dtMax, const double m0, const double *ak0, const double *rhok0, const double *pk0) const; // Riemann pour injection
    //virtual void resolRiemannRes(Cellule &cellGauche, const int &nombrePhases, const double &dxGauche, double &dtMax, const double *ak0, const double *rhok0, const double *pk0) const; // Riemann pour reservoir
    //virtual void resolRiemannSortie(Cellule &cellGauche, const int &nombrePhases, const double &dxGauche, double &dtMax, const double p0, double *debitSurf) const; // Riemann pour sortie a pression imposee
    virtual double getSM();

    virtual void projectionRepereAbsolu(const Coord normale, const Coord tangente, const Coord binormale) const;

	  // Methodes pour transport
	  //KS//DEV// Pas encore fait pour EulerHomogene

    //Accesseur
    int getLiq();
    int getVap();
    virtual std::string quiSuisJe() const;

  protected:

  private:
    static const std::string NOM;
    int m_liq;
    int m_vap;

    friend class FluxEulerHomogene;
};

#endif // MODEULERHOMOGENE_H
