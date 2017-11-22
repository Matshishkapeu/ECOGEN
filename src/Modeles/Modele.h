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

#ifndef MODELE_H
#define MODELE_H

class Modele; //Predeclaration de la classe Modele pour pouvoir inclure Flux.h

#include "Flux.h"
#include "../Maths/Coord.h"
#include "../Erreurs.h"

class Modele
{
  public:
    /** Default constructor */
    Modele();
    Modele(const std::string &nom, const int &nombreTransports);
    /** Default destructor */
    virtual ~Modele();

    virtual void alloueCons(Flux **cons, const int &nombrePhases) { Erreurs::messageErreur("alloueCons non implemente pour modele utilise"); };
    virtual void allouePhase(Phase **phase) { Erreurs::messageErreur("allouePhase non implemente pour modele utilise"); };
    virtual void alloueMelange(Melange **melange) { Erreurs::messageErreur("alloueMelange non implemente pour modele utilise"); };
    virtual void alloueEos(Cellule &cell, const int &nombrePhases) { Erreurs::messageErreur("alloueEos non implemente pour modele utilise"); };

    //Methodes pour hydrodynamique
    virtual void resolRiemannInterne(Cellule &cellGauche, Cellule &cellDroite, const int &nombrePhases, const double &dxGauche, const double &dxDroite, double &dtMax) const { Erreurs::messageErreur("resolRiemannInterne non implemente pour modele utilise"); }; // Riemann entre deux maille de calcul
    virtual void resolRiemannMur(Cellule &cellGauche, const int &nombrePhases, const double &dxGauche, double &dtMax) const { Erreurs::messageErreur("resolRiemannMur non implemente pour modele utilise"); }; // Riemann entre maille gauche et mur
    virtual void resolRiemannInj(Cellule &cellGauche, const int &nombrePhases, const double &dxGauche, double &dtMax, const double m0, const double *ak0, const double *rhok0, const double *pk0) const { Erreurs::messageErreur("resolRiemannInj non implemente pour modele utilise"); }; // Riemann pour injection
    virtual void resolRiemannRes(Cellule &cellGauche, const int &nombrePhases, const double &dxGauche, double &dtMax, const double *ak0, const double *rhok0, const double *pk0) const { Erreurs::messageErreur("resolRiemannRes non implemente pour modele utilise"); }; // Riemann pour reservoir
    virtual void resolRiemannSortie(Cellule &cellGauche, const int &nombrePhases, const double &dxGauche, double &dtMax, const double p0, double *debitSurf) const { Erreurs::messageErreur("resolRiemannSortie non implemente pour modele utilise"); }; // Riemann pour sortie

    //Methodes pour transport
    virtual void resolRiemannTransportInterne(Cellule &cellGauche, Cellule &cellDroite, const int &nombreTransports) { Erreurs::messageErreur("resolRiemannTransportInterne non implemente pour modele utilise"); };
    virtual void resolRiemannTransportMur(const int &nombreTransports) { Erreurs::messageErreur("resolRiemannTransportMur non implemente pour modele utilise"); };
    virtual void resolRiemannTransportInj(Cellule &cellGauche, const int &nombreTransports, double *valeurTransports) { Erreurs::messageErreur("resolRiemannTransportInj non implemente pour modele utilise"); };
    virtual void resolRiemannTransportRes(Cellule &cellGauche, const int &nombreTransports, double *valeurTransports) { Erreurs::messageErreur("resolRiemannTransportRes non implemente pour modele utilise"); };
    virtual void resolRiemannTransportSortie(Cellule &cellGauche, const int &nombreTransports, double *m_valeurTransport) { Erreurs::messageErreur("resolRiemannTransportSortie non implemente pour modele utilise"); };
    virtual double getSM() { Erreurs::messageErreur("getSM non implemente pour modele utilise"); return 0;};

    virtual void projectionRepereAbsolu(const Coord normale, const Coord tangente, const Coord binormale) const { Erreurs::messageErreur("projectionRepereAbsolu non implemente pour modele utilise"); };

    //Accesseur
    void afficheInfos() const;
    virtual std::string quiSuisJe() const { return 0; };

  protected:
    std::string m_nom;
    
  private:
};

#endif // MODELE_H
