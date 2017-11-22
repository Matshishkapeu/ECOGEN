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

#ifndef TRANSPORT_H
#define TRANSPORT_H

#include "../libTierces/tinyxml2.h"
#include <fstream>

class Transport
{
  public:
    Transport();
    virtual ~Transport();

    void setValeur(double valeur);
    double getValeur() const;

   	void resolRiemann(double transportGauche, double transportDroite, double sM);
    void resolRiemannMur();
    void resolRiemannInj(double transportGauche, double sM, double valeurTransport);
    void resolRiemannRes(double transportGauche, double sM, double valeurTransport);
    void resolRiemannSortie(double transportGauche, double sM, double valeurTransport);
    void ajoutFlux(double coefA, const int num);
    void retireFlux(double coefA, const int num);
    void ajoutNonCons(double coefA, double transport, const double sM);
    void retireNonCons(double coefA, double transport, const double sM);
    void multiplie(double scalaire);
    void ajoute(double scalaire);
    void changeSigne();
    
    //Specifique ordre2
    void calculPentesTransport(const double valeurGauche, const double valeurDroite, const double &distance);
    void extrapole(const double &pente, const double &distance);

  private:
    double m_valeur;
};

extern Transport* fluxTempTransport;

#endif // TRANSPORT_H