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

#ifndef PAKCAPILLARITE_H
#define PAKCAPILLARITE_H

#include "../PAKapila.h"
#include "GPACapillarite.h"

class PAKCapillarite : public PAKapila
{
  public:
    PAKCapillarite();
    PAKCapillarite(tinyxml2::XMLElement *element, int& nombreGPA, std::vector<std::string> nomTransports, std::string nomFichier = "Fichier Inconnu");
    virtual ~PAKCapillarite();

    virtual void ajouteGrandeurPhysAdd(Cellule *cell);

    virtual double calculEnergiePhysAdd(GrandeursPhysAdd* GPA);
    virtual void resolFluxPhysAdd(BordDeMaille *bord, const int &nombrePhases);
    virtual void resolFluxPhysAddLimite(BordDeMaille *bord, const int &nombrePhases);
    void resolFluxCapillariteInterne(Coord &vitesseGauche, Coord &vitesseDroite, Coord &gradCGauche, Coord &gradCDroite) const;
    void resolFluxCapillariteAbs(Coord &vitesseGauche, Coord &gradCGauche) const;
    void resolFluxCapillariteMur(Coord &gradCGauche) const;
    void resolFluxCapillariteSortie(Coord &vitesseGauche, Coord &gradCGauche) const;
    void resolFluxCapillariteInjection(Coord &vitesseGauche, Coord &gradCGauche) const;
    void resolFluxCapillariteAutres(Coord &vitesseGauche, Coord &gradCGauche) const;
    virtual void ajoutNonCons(Cellule *cell, const int &nombrePhases) {}; //La capillarite n a pas de terme non cons.

		virtual void reinitialiseFonctionCouleur(std::vector<Cellule *> *cellulesLvl, int &lvl);

    virtual void communicationsPhysAdd(Cellule **cellules, const int &dim);
		virtual void communicationsPhysAddAMR(Cellule **cellules, const int &dim, const int &lvl);
    virtual int getNumTransportAssocie() const;

  protected:
  
  private:
    std::string m_nomTransportAssocie;  //!< nom de la variable assosiee au sein de chaque cellule (m_vecTransport) 
    int m_numTransportAssocie;          //!< numero de la variable assosiee au sein de chaque cellule (m_vecTransport) 
    int m_numGPAGradC;                  //!< numero de la variable assosiee au sein de chaque cellule (m_vecGrandeursPhysAdd)
    double m_sigma;                     //!< tension de surface

    Coord m_vitesseGauche;
    Coord m_gradCGauche;
    Coord m_vitesseDroite;
    Coord m_gradCDroite;
    Coord m_normale;
    Coord m_tangente;
    Coord m_binormale;
};

#endif // PAKCAPILLARITE_H
