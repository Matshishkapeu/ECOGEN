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

#ifndef CONDLIMINJ_H
#define CONDLIMINJ_H

#include "CondLim.h"


class CondLimInj : public CondLim
{
  public:
    CondLimInj();
    CondLimInj(int numPhysique, tinyxml2::XMLElement *element, std::vector<Phase*> vecPhases, int &nombreTransports, std::vector<std::string> nomTransports, std::string nomFichier = "Fichier Inconnu");
    CondLimInj(double m0, double *ak0, double *rho0, double *p0, const int nombrePhases = 1);
    CondLimInj(const CondLimInj &Source, const int lvl = 0); //Constructeur de copie (utile pour AMR)
    virtual ~CondLimInj();

    virtual void creeLimite(BordDeMaille **face);
    virtual void resolRiemannLimite(Cellule &cellGauche, const int &nombrePhases, const double &dxGauche, double &dtMax);
    virtual void resolRiemannTransportLimite(Cellule &cellGauche, const int &nombreTransports) const;

    virtual int quiSuisJe() const { return 4; };
    virtual void afficheInfos();

    //Pour methode AMR
    virtual void creerBordEnfant();  /*!< Creer un bord enfant (non initialise) */

  protected:
  private:
    int m_nombrePhase;
    double m_m0;     //!< debit massique surfacique
    double *m_ak0;
    double *m_rhok0;
    double *m_pk0;
    int m_nombreTransports;
    double *m_valeurTransport;

};

#endif // CONDLIMINJ_H
