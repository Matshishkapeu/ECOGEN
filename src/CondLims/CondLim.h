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

#ifndef CONDLIM_H
#define CONDLIM_H

#include <iostream>
#include "../BordDeMaille.h"
#include "../Ordre2/BordDeMailleO2.h" //Ajouter pour l'AMR, a priori ne pose pas de probleme
#include "../libTierces/tinyxml2.h"
#include "../Erreurs.h"
#include "../Outils.h"

class CondLim : public BordDeMaille
{
  public:
    CondLim();
    CondLim(int numPhysique);
    CondLim(const CondLim &Source);
    virtual ~CondLim();

    virtual void creeLimite(BordDeMaille **face){ Erreurs::messageErreur("Impossible de creer la limite dans creeLimite"); };
    virtual void creeLimite(BordDeMaille **face, std::string ordreCalcul) { Erreurs::messageErreur("Impossible de creer la limite dans creeLimite"); };
    virtual void initialise(Cellule *cellGauche, Cellule *cellDroite);

    virtual void calculFlux(const int &nombrePhases, const int &nombreTransports, double &dtMax, Limiteur &limiteurGlobal, Limiteur &limiteurInterface, Prim type = vecPhases);
    virtual void calculFluxPhysAdd(const int &nombrePhases, PhysAdd &physAdd);
    virtual void resolRiemann(const int &nombrePhases, const int &nombreTransports, double &dtMax, Limiteur &limiteurGlobal, Limiteur &limiteurInterface, Prim type = vecPhases);
    virtual void ajoutFlux(const int &nombrePhases, const int &nombreTransports, const double &coefAMR, Prim type = vecPhases) {};  //Ici la fonction ne fait rien car il s agit d une limite a droite et il n y a rien a ajouter a droite.
    virtual void resolRiemannLimite(Cellule &cellGauche, const int &nombrePhases, const double &dxGauche, double &dtMax) { Erreurs::messageErreur("Attention resolRiemannLimite non prevu pour limite utilisee"); };
    virtual void resolRiemannTransportLimite(Cellule &cellGauche, const int &nombreTransports) const { Erreurs::messageErreur("Attention resolRiemannTransportLimite non prevu pour limite utilisee"); };

    virtual int quiSuisJe() const { Erreurs::messageErreur("quiSuisJe pas prevu pour la limite demandee"); return 0; };
    virtual void afficheInfos(){};

    virtual int getNumPhys() const;

    //Pour methode AMR
    virtual void calculXi(const double &critereVar, const bool &varRho, const bool &varP, const bool &varU, const bool &varAlpha) {};
    virtual void calculFluxXi(const double &dX, const double &dY, const double &dZ) {};
    virtual void raffineBordExterne(const int &nbMaillesY, const int &nbMaillesZ, const double &dXParent, const double &dYParent,
      const double &dZParent, Cellule *cellRef, const double &surfaceEnfant);
    virtual void deraffineBordExterne(Cellule *cellRef);

  protected:
    int m_numPhysique;

  private:
};

#endif // CONDLIM_H
