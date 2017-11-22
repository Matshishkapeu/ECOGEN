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

#ifndef CONDLIMMURO2_H
#define CONDLIMMURO2_H

#include "CondLimMur.h"


class CondLimMurO2 : public CondLimMur
{
public:
  CondLimMurO2();
  CondLimMurO2(const CondLimMurO2& Source, const int lvl = 0); //Constructeur de copie (utile pour AMR)
  CondLimMurO2(int numPhysique);
  virtual ~CondLimMurO2();

  virtual void creeLimite(BordDeMaille **face);
  virtual void allouePentes(const int &nombrePhases, const int &nombreTransports, int &allouePenteLocal);
  virtual void calculPentes(const int &nombrePhases, const int &nombreTransports, Prim type = vecPhases);
  virtual void resolRiemann(const int &nombrePhases, const int &nombreTransports, double &dtMax, Limiteur &limiteurGlobal, Limiteur &limiteurInterface, Prim type = vecPhases);

  //Accesseurs
  virtual Phase* getPentesPhase(const int &numeroPhase) const;
  virtual Melange* getPentesMelange() const;
  virtual Transport* getPentesTransport(const int &numeroTransport) const;

  //Pour methode AMR
  virtual void creerBordEnfant();  /*!< Creer un bord enfant (non initialise) */

protected:
  int m_nombrePhases;
  Phase **m_vecPhasesPentes;         /*!< vecteur des pentes des phases */
  Melange *m_melangePentes;          /*!< vecteur des pentes de melange */
  Transport *m_vecTransportsPentes;	/*!< vecteur des pentes des transports */
private:
};

#endif // CONDLIMMURO2_H
