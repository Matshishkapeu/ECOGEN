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

#ifndef PAKCONDUCTIVITE_H
#define PAKCONDUCTIVITE_H

#include "../PAKapila.h"
#include "GPAConductivite.h"
#include "../../Eos/Eos.h"

class PAKConductivite : public PAKapila
{
  public:
    PAKConductivite();
    PAKConductivite(int& nombreGPA, Eos** eos, int &nombrePhases, std::string nomFichier = "Fichier Inconnu");
    virtual ~PAKConductivite();

    virtual void ajouteGrandeurPhysAdd(Cellule *cell);

    virtual void resolFluxPhysAdd(BordDeMaille *bord, const int &nombrePhases);
    virtual void resolFluxPhysAddLimite(BordDeMaille *bord, const int &nombrePhases);
    void resolFluxConductiviteInterne(Coord &gradTkGauche, Coord &gradTkDroite, double &alphakL, double &alphakR, int &numPhase) const;
    void resolFluxConductiviteAbs(Coord &gradTkGauche, double &alphakL, int &numPhase) const;
    void resolFluxConductiviteMur(Coord &gradTkGauche, double &alphakL, int &numPhase) const;
    void resolFluxConductiviteSortie(Coord &gradTkGauche, double &alphakL, int &numPhase) const;
    void resolFluxConductiviteInjection(Coord &gradTkGauche, double &alphakL, int &numPhase) const;
    void resolFluxConductiviteAutres(Coord &gradTkGauche, double &alphakL, int &numPhase) const;
    virtual void ajoutNonCons(Cellule *cell, const int &nombrePhases) {}; //La conductivite n a pas de terme non cons.

    virtual void communicationsPhysAdd(Cellule **cellules, const int &dim);
		virtual void communicationsPhysAddAMR(Cellule **cellules, const int &dim, const int &lvl);

    //Accesseurs
    virtual double getLambdak(int &numPhase) const;

  protected:

  private:
    double *m_lambdak;
    int m_numGPA;

    Coord m_gradTkGauche;
    Coord m_gradTkDroite;
    Coord m_normale;
    Coord m_tangente;
    Coord m_binormale;
};

#endif // PAKCONDUCTIVITE_H
