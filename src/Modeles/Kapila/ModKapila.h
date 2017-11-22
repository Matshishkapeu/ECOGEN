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

#ifndef MODKAPILA_H
#define MODKAPILA_H

#include "../Modele.h"
#include "../../Cellule.h"
#include "FluxKapila.h"
#include "MelKapila.h"

class ModKapila : public Modele
{
  public:
    /** Default constructor */
    ModKapila(int &nombreTransports, const int &nombrePhases);
    /** Default destructor */
    virtual ~ModKapila();

    virtual void alloueCons(Flux **cons, const int &nombrePhases);
    virtual void allouePhase(Phase **phase);
    virtual void alloueMelange(Melange **melange);
    virtual void alloueEos(Cellule &cell, const int &nombrePhases);

    //Methodes pour hydrodynamique
    virtual void resolRiemannInterne(Cellule &cellGauche, Cellule &cellDroite, const int &nombrePhases, const double &dxGauche, const double &dxDroite, double &dtMax) const; // Riemann entre deux maille de calcul
    virtual void resolRiemannMur(Cellule &cellGauche, const int &nombrePhases, const double &dxGauche, double &dtMax) const; // Riemann entre maille gauche et mur
    virtual void resolRiemannInj(Cellule &cellGauche, const int &nombrePhases, const double &dxGauche, double &dtMax, const double m0, const double *ak0, const double *rhok0, const double *pk0) const; // Riemann pour injection
    virtual void resolRiemannRes(Cellule &cellGauche, const int &nombrePhases, const double &dxGauche, double &dtMax, const double *ak0, const double *rhok0, const double *pk0) const; // Riemann pour reservoir
    virtual void resolRiemannSortie(Cellule &cellGauche, const int &nombrePhases, const double &dxGauche, double &dtMax, const double p0, double *debitSurf) const; // Riemann pour sortie a pression imposee

    // Methodes pour transport
    virtual void resolRiemannTransportInterne(Cellule &cellGauche, Cellule &cellDroite, const int &nombreTransports);
    virtual void resolRiemannTransportMur(const int &nombreTransports);
    virtual void resolRiemannTransportInj(Cellule &cellGauche, const int &nombreTransports, double *valeurTransports);
    virtual void resolRiemannTransportRes(Cellule &cellGauche, const int &nombreTransports, double *valeurTransports);
    virtual void resolRiemannTransportSortie(Cellule &cellGauche, const int &nombreTransports, double *m_valeurTransport);
    virtual double getSM();

    virtual void projectionRepereAbsolu(const Coord normale, const Coord tangente, const Coord binormale) const;

    //Accesseur
    virtual std::string quiSuisJe() const;

  protected:
  
  private:
    static const std::string NOM;
    double* m_rhokStar;  //!< tableaux temporaires utiles pour la resolution des pb de Riemann
    double* m_pkStar;    //!< tableaux temporaires utiles pour la resolution des pb de Riemann
    double* m_ekStar;    //!< tableaux temporaires utiles pour la resolution des pb de Riemann
    double* m_EkStar;    //!< tableaux temporaires utiles pour la resolution des pb de Riemann
    double* m_vkStar;    //!< tableaux temporaires utiles pour la resolution des pb de Riemann
    double* m_YkStar;    //!< tableaux temporaires utiles pour la resolution des pb de Riemann
    double* m_Hk0;       //!< tableaux temporaires utiles pour la resolution des pb de Riemann
    double* m_Yk0;       //!< tableaux temporaires utiles pour la resolution des pb de Riemann
};

#endif // MODKAPILA_H
