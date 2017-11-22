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

#ifndef PHYSADD_H
#define PHYSADD_H

class PhysAdd; //Predeclaration de la classe PhysAdd pour pouvoir inclure Cellule.h et BordDeMaille.h

#include "../Erreurs.h"
#include "../Cellule.h"
#include "../BordDeMaille.h"
#include "GrandeursPhysAdd.h"
#include "../Parallel.h"

class PhysAdd
{
  public:
    PhysAdd();
    virtual ~PhysAdd();

    virtual void ajouteGrandeurPhysAdd(Cellule *cell) { Erreurs::messageErreur("ajouteGrandeurPhysAdd non implemente pour physique additionnelle utilise"); };

    virtual std::string quiSuisJe() const { Erreurs::messageErreur("quiSuisJe non implemente pour physique additionnelle utilise"); return 0; };

    virtual double calculEnergiePhysAdd(GrandeursPhysAdd* GPA) { return 0.;  }; //!< renvoi l energie massique lie a la physique (0 si pas d energie associee)
    void calculFluxPhysAdd(BordDeMaille *bord, const int &nombrePhases);
    void calculFluxPhysAddLimite(BordDeMaille *bord, const int &nombrePhases);
    void ajoutNonConsPhysAdd(Cellule *cell, const int &nombrePhases);
    virtual void resolFluxPhysAdd(BordDeMaille *bord, const int &nombrePhases) { Erreurs::messageErreur("resolFluxPhysAdd non implemente pour physique additionnelle utilise"); };
    virtual void resolFluxPhysAddLimite(BordDeMaille *bord, const int &nombrePhases) { Erreurs::messageErreur("resolFluxPhysAddLimite non implemente pour physique additionnelle utilise"); };
    virtual void ajoutFluxPhysAdd(BordDeMaille *bord, const int &nombrePhases, const double &coefAMR);
    virtual void retireFluxPhysAdd(BordDeMaille *bord, const int &nombrePhases, const double &coefAMR);
    virtual void ajoutNonCons(Cellule *cell, const int &nombrePhases) { Erreurs::messageErreur("ajoutNonCons non implemente pour physique additionnelle utilise"); };

		virtual void reinitialiseFonctionCouleur(std::vector<Cellule *> *cellulesLvl, int &lvl) {}; //!< Utile seulement pour la capillarite

    virtual void communicationsPhysAdd(Cellule **cellules, const int &dim) { Erreurs::messageErreur("communicationsPhysAdd non implemente pour physique additionnelle utilise"); };
		virtual void communicationsPhysAddAMR(Cellule **cellules, const int &dim, const int &lvl) { Erreurs::messageErreur("communicationsPhysAddAMR non implemente pour physique additionnelle utilise"); };
    virtual int getNumTransportAssocie() const { Erreurs::messageErreur("getNumTransportAssocie non implemente pour physique additionnelle utilise"); return 0; };

  protected:
    
  private:
};

#endif // PHYSADD_H
