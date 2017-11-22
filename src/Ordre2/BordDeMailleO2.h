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

#ifndef BORDDEMAILLEO2_H
#define BORDDEMAILLEO2_H

#include "../BordDeMaille.h"

class BordDeMailleO2; //Predeclaration de la classe BordDeMailleO2 pour pouvoir inclure CelluleO2.h

#include "CelluleO2.h"

class BordDeMailleO2 : public BordDeMaille
{
  public:
    /** Default constructor */
    BordDeMailleO2();
    BordDeMailleO2(int lvl); //Pour AMR
    /** Default destructor */
    virtual ~BordDeMailleO2();

    virtual void allouePentes(const int &nombrePhases, const int &nombreTransports, int &allouePenteLocal);
    virtual void calculPentes(const int &nombrePhases, const int &nombreTransports, Prim type = vecPhases);
    virtual void calculFlux(const int &nombrePhases, const int &nombreTransports, double &dtMax, Limiteur &limiteurGlobal, Limiteur &limiteurInterface, Prim type = vecPhases);
    void resolRiemann(const int &nombrePhases, const int &nombreTransports, double &ondeMax, Limiteur &limiteurGlobal, Limiteur &limiteurInterface, Prim type = vecPhases); /*!< probleme de Riemann special ordre 2 */

    //Accesseurs
    virtual Phase* getPentesPhase(const int &numeroPhase) const;
    virtual Melange* getPentesMelange() const;
    virtual Transport* getPentesTransport(const int &numeroTransport) const;
    //virtual Cellule *getB(BO2 B) const;
    //virtual double getBeta(betaO2 beta) const;
    //virtual double getDistanceH(distanceHO2 dist) const;
    //virtual void setB(BO2 B, Cellule *cellule);
    //void setBeta(betaO2 beta, double &valeur);
    //virtual void setDistanceH(distanceHO2 dist, double &valeur);

    //Pour methode AMR
    virtual void creerBordEnfant();                                                                        /*!< Creer un bord enfant (non initialise) */
    virtual void creerBordEnfantInterne(const int &lvl, std::vector<BordDeMaille*> *bordsEnfantsInternes); /*!< Creer un bord enfant interne (non initialise) */

   protected:
     int m_nombrePhases;
     Phase **m_vecPhasesPentes;         /*!< vecteur des pentes des phases */
     Melange *m_melangePentes;          /*!< vecteur des pentes de melange */
     Transport *m_vecTransportsPentes;	/*!< vecteur des pentes des transports */

     //Stockage methode multipentes
     //Cellule *m_BG1M; /*!< pointeurs vers cellules Arrieres a gauche pour ordre2  */
     //Cellule *m_BG2M;
     //Cellule *m_BG3M;
     //Cellule *m_BG1P; /*!< pointeurs vers cellules Avants a gauche pour ordre2  */
     //Cellule *m_BG2P;
     //Cellule *m_BG3P;
     //Cellule *m_BD1M; /*!< pointeurs vers cellules Arrieres a droite pour ordre2  */
     //Cellule *m_BD2M;
     //Cellule *m_BD3M;
     //Cellule *m_BD1P; /*!< pointeurs vers cellules Avants a droite pour ordre2  */
     //Cellule *m_BD2P;
     //Cellule *m_BD3P;

     //double m_betaG1M;  /*!< ponderations pour ordre2 */
     //double m_betaG2M;
     //double m_betaG3M;
     //double m_betaG1P;
     //double m_betaG2P;
     //double m_betaG3P;
     //double m_betaD1M;
     //double m_betaD2M;
     //double m_betaD3M;
     //double m_betaD1P;
     //double m_betaD2P;
     //double m_betaD3P;

     //double m_distanceHGM;  /*!< distances au point geometrique pour le calcul des pentes */
     //double m_distanceHGP;
     //double m_distanceHDM;
     //double m_distanceHDP;

   private:
};

extern Phase **pentesPhasesLocal1;
extern Phase **pentesPhasesLocal2;
extern Melange *pentesMelangeLocal1;
extern Melange *pentesMelangeLocal2;
extern double *pentesTransportLocal1;
extern double *pentesTransportLocal2;

#endif // BORDDEMAILLEO2_H
