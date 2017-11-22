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

#ifndef PAKVISCOSITE_H
#define PAKVISCOSITE_H

#include "../PAKapila.h"
#include "GPAViscosite.h"
#include "../../Eos/Eos.h"

class PAKViscosite : public PAKapila
{
  public:
    PAKViscosite();
    PAKViscosite(int& nombreGPA, Eos** eos, int &nombrePhases, std::string nomFichier = "Fichier Inconnu");
    virtual ~PAKViscosite();

    virtual void ajouteGrandeurPhysAdd(Cellule *cell);

    virtual void resolFluxPhysAdd(BordDeMaille *bord, const int &nombrePhases);
    virtual void resolFluxPhysAddLimite(BordDeMaille *bord, const int &nombrePhases);
    void resolFluxViscositeInterne(Coord &vitesseGauche, Coord &vitesseDroite, Coord &gradUGauche, Coord &gradUDroite, Coord &gradVGauche, Coord &gradVDroite, Coord &gradWGauche, Coord &gradWDroite, double &muMelGauche, double &muMelDroite, double &distGauche, double &distDroite, int nombrePhases) const;
    void resolFluxViscositeAbs(Coord &vitesseGauche, Coord &gradUGauche, Coord &gradVGauche, Coord &gradWGauche, double &muMelGauche, double &distGauche, int nombrePhases) const;
    void resolFluxViscositeMur(Coord &vitesseGauche, double &muMelGauche, int nombrePhases, double &distGauche) const;
    void resolFluxViscositeAutres(Coord &vitesseGauche, Coord &gradUGauche, Coord &gradVGauche, Coord &gradWGauche, double &muMelGauche, double &distGauche, int nombrePhases) const;
    virtual void ajoutNonCons(Cellule *cell, const int &nombrePhases);

    virtual void communicationsPhysAdd(Cellule **cellules, const int &dim);
		virtual void communicationsPhysAddAMR(Cellule **cellules, const int &dim, const int &lvl);

  protected:
  
  private:
    double *m_muk;
    int m_numGPA;  //!< numero de la variable assosiee au sein de chaque cellule (m_vecGrandeursPhysAdd)

    Coord m_vitesseGauche;
    Coord m_gradUGauche;
    Coord m_gradVGauche;
    Coord m_gradWGauche;
    Coord m_vitesseDroite;
    Coord m_gradUDroite;
    Coord m_gradVDroite;
    Coord m_gradWDroite;
    Coord m_normale;
    Coord m_tangente;
    Coord m_binormale;
};

#endif // PAKVISCOSITE_H
