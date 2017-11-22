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

#ifndef CELLULEO2_H
#define CELLULEO2_H

#include "../Cellule.h"
#include "BordDeMailleO2.h"

class CelluleO2 : public Cellule
{
    public:
        CelluleO2();
        CelluleO2(int lvl); //Pour AMR
        virtual ~CelluleO2();
        virtual void alloue(const int &nombrePhases, const int &nombreTransports, const std::vector<PhysAdd*> &physAdd, Modele *modele);
        virtual void copiePhase(const int &numeroPhase, Phase *phase);
        virtual void calculPentesLocal(const int &nombrePhases, const int &nombreTransports, BordDeMaille &bordRef, Limiteur &limiteurGlobal, Limiteur &limiteurInterface);
        virtual void calculPentesLocalLimite(const int &nombrePhases, const int &nombreTransports, BordDeMaille &bordRef, Limiteur &limiteurGlobal, Limiteur &limiteurInterface);
        //virtual void calculMultiPente(const int &nombrePhases, BordDeMaille *bord, Limiteur *limiteurGlobal);
        virtual void sauvegardeCons(const int &nombrePhases, const int &nombreTransports);
        virtual void recuperationCons(const int &nombrePhases, const int &nombreTransports);
        virtual void predictionOrdre2(const double &dt, const int &nombrePhases, const int &nombreTransports);
        virtual void alloueEtCopiePhase(const int &numeroPhase, Phase *phase);
        virtual void calculsEtendus(const int &nombrePhases, Prim type = vecPhases);
        virtual void calculsEtendusPourCommunications(const int &nombrePhases, Prim type = vecPhases);
        virtual void projection(const Coord &normale, const Coord &tangente, const Coord &binormale, const int &nombrePhases, Prim type = vecPhases);
        virtual void projectionRepereAbsolu(const Coord &normale, const Coord &tangente, const Coord &binormale, const int &nombrePhases, Prim type = vecPhases);
        virtual void copieDansCellule(Cellule &cellSource, Prim type=vecPhases) const;

        //Accesseurs
        virtual Phase* getPhase(const int &numeroPhase, Prim type = vecPhases) const;
        virtual Phase** getPhases(Prim type = vecPhases) const;
        virtual Melange* getMelange(Prim type = vecPhases) const;
        virtual Transport& getTransport(const int &numTransport, Prim type = vecPhases) const;
        virtual Transport* getTransports(Prim type = vecPhases) const;
        virtual void setTransport(double valeur, int &numTransport, Prim type = vecPhases);

        //Pour methode AMR
        virtual void creerCelluleEnfant(const int &num, const int &lvl);                                              /*!< Creer une cellule enfant (non initialisee) */

				//Pour methodes ordre 2 parallele
				virtual void rempliTamponPentes(double *tampon, int &compteur) const;
				virtual void rempliTamponPentesAMRjeSuisCpuGauche(double *tampon, int &compteur, const int &lvl) const;
				virtual void rempliTamponPentesAMRjeSuisCpuDroite(double *tampon, int &compteur, const int &lvl) const;

    protected:
        Phase **m_vecPhasesO2;             /*!< pour stocker les valeurs predites a l ordre 2 */
        Melange *m_melangeO2;              /*!< pour stocker les valeurs predites a l ordre 2 */
        Transport *m_vecTransportsO2;		   /*!< pour stocker les valeurs predites a l ordre 2 */
        Flux *m_consSauvegarde;            /*!< Vecteur de sauvegarde des variables conservatives. De type flux car recueille la somme des flux sur l objet cellule */
        Transport *m_consTransportsSauvegarde;  /*!< Vecteur de saugevarde des grandeurs passives permettant de recueillir la somme des flux des grandeurs transportees */

    private:
};

#endif // CELLULEO2_H
