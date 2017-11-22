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

#ifndef CONDLIMSORTIE_H
#define CONDLIMSORTIE_H

#include "CondLim.h"


class CondLimSortie : public CondLim
{
  public:
    CondLimSortie();
    CondLimSortie(int numPhysique, tinyxml2::XMLElement *element, int &nombrePhases, int &nombreTransports, std::vector<std::string> nomTransports, std::string nomFichier);
    CondLimSortie(double p0);
    CondLimSortie(const CondLimSortie& Source, const int lvl = 0); //Constructeur de copie (utile pour AMR)
    virtual ~CondLimSortie();

    virtual void creeLimite(BordDeMaille **face);
    virtual void resolRiemannLimite(Cellule &cellGauche, const int &nombrePhases, const double &dxGauche, double &dtMax);
    virtual void resolRiemannTransportLimite(Cellule &cellGauche, const int &nombreTransports) const;

    virtual int quiSuisJe() const { return 3; };
    virtual void afficheInfos();

    //Accesseur
    //virtual double getDebit(int numPhase) const;

    //Pour methode AMR
    virtual void creerBordEnfant();  /*!< Creer un bord enfant (non initialise) */

  protected:
  private:
    double m_p0;
    int m_nombreTransports;
    int m_nombrePhases;
    double *m_valeurTransport;

    double *m_debits;
};

#endif // CONDLIMSORTIE_H
