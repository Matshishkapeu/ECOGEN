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

#ifndef BORDDEMAILLE_H
#define BORDDEMAILLE_H

class BordDeMaille; //Predeclaration de la classe BordDeMaille pour pouvoir inclure Cellule.h

#include "Cellule.h"
#include "Modeles/Modele.h"
#include "Modeles/Flux.h"
#include "Maths/Coord.h"
#include "Maillages/Face.h"
#include "Maillages/FaceCartesien.h"
#include "PhysiqueAdditionnelle/PhysAdd.h"

enum BO2 { BG1M, BG2M, BG3M, BG1P, BG2P, BG3P, BD1M, BD2M, BD3M, BD1P, BD2P, BD3P };
enum betaO2 { betaG1M, betaG2M, betaG3M, betaG1P, betaG2P, betaG3P, betaD1M, betaD2M, betaD3M, betaD1P, betaD2P, betaD3P };
enum distanceHO2 { distanceHGM, distanceHGP, distanceHDM, distanceHDP };

class BordDeMaille
{
  public:
    /** Default constructor */
    BordDeMaille();
    BordDeMaille(int lvl); //Pour AMR
    /** Default destructor */
    virtual ~BordDeMaille();

    void setFace(Face *face);

    virtual void calculFlux(const int &nombrePhases, const int &nombreTransports, double &dtMax, Limiteur &limiteurGlobal, Limiteur &limiteurInterface, Prim type = vecPhases);
    virtual void calculFluxPhysAdd(const int &nombrePhases, PhysAdd &physAdd);
    virtual void resolRiemann(const int &nombrePhases, const int &nombreTransports, double &ondeMax, Limiteur &limiteurGlobal, Limiteur &limiteurInterface, Prim type = vecPhases);
    virtual void initialise(Cellule *cellGauche, Cellule *cellDroite);
    void initialiseGauche(Cellule *cellGauche);
    virtual void initialiseDroite(Cellule *cellDroite);
    virtual void ajoutFlux(const int &nombrePhases, const int &nombreTransports, const double &coefAMR);
    void retireFlux(const int &nombrePhases, const int &nombreTransports, const double &coefAMR);
    double distance(Cellule *c);

    void EffetsSurface1D(const int &nombrePhases);

    void associeModele(Modele *mod);

		virtual int quiSuisJe() const {	return 0; };

    //Inutilise pour Bord de Maille ordre 1
    virtual void allouePentes(const int &nombrePhases, const int &nombreTransports, int &allouePenteLocal) {};   /*!< Ne fait rien pour des Bord de Maille ordre 1 */
    virtual void calculPentes(const int &nombrePhases, const int &nombreTransports, Prim type = vecPhases) {};   /*!< Ne fait rien pour des Bord de Maille ordre 1 */
    virtual Phase* getPentesPhase(const int &numeroPhase) const { return 0; };                                   /*!< Ne fait rien pour des Bord de Maille ordre 1 */
    virtual Melange* getPentesMelange() const { return 0; };                                                     /*!< Ne fait rien pour des Bord de Maille ordre 1 */
    virtual Transport* getPentesTransport(const int &numeroTransport) const { return 0; };                       /*!< Ne fait rien pour des Bord de Maille ordre 1 */
    //virtual Cellule *getB(BO2 B) const { return 0; };                                                          /*!< Ne fait rien pour des Bord de Maille ordre 1 */
    //virtual double getBeta(betaO2 beta) const { return 0.; };                                                  /*!< Ne fait rien pour des Bord de Maille ordre 1 */
    //virtual double getDistanceH(distanceHO2 dist) const { return 0.; };                                        /*!< Ne fait rien pour des Bord de Maille ordre 1 */
    //virtual void setB(BO2 B, Cellule *cellule) {};                                                             /*!< Ne fait rien pour des Bord de Maille ordre 1 */
    //virtual void setBeta(betaO2 beta, double &valeur) {};                                                      /*!< Ne fait rien pour des Bord de Maille ordre 1 */
    //virtual void setDistanceH(distanceHO2 dist, double &valeur) {};                                            /*!< Ne fait rien pour des Bord de Maille ordre 1 */


    //Accesseurs
    Face *getFace();                                            /*!< Attention, getFace() non const */
    Modele *getMod() const;
    Cellule *getCellGauche() const;
    Cellule *getCellDroite() const;
    virtual int getNumPhys() const { return -1; };
    //virtual double getDebit(int numPhase) const { Erreurs::messageErreur("getDebits non prevu pour BordDeMaille"); return 0.; }

    //Pour methode AMR
    virtual void calculXi(const double &critereVar, const bool &varRho, const bool &varP, const bool &varU, const bool &varAlpha);  /*!< Calcul de la variable Xi pour critere de (de)raffinement a priori */
    void calculCritereAMR(const double &critereVar, std::string nomVariable, int num = 0);                                          /*!< Calcul de xi via le critere de variation */
    virtual void calculFluxXi(const double &dXLocal, const double &dYLocal, const double &dZLocal);                                 /*!< Calcul des flux de Xi (diffusion) pour smoothing */
    virtual void creerBordEnfant();                                                                        /*!< Creer un bord enfant (non initialise) */
    virtual void creerBordEnfantInterne(const int &lvl, std::vector<BordDeMaille*> *bordsEnfantsInternes); /*!< Creer un bord enfant interne (non initialise) */
    void creerFaceEnfant(BordDeMaille *bordParent);             /*!< Creer une face enfant (non initialise) */
    virtual void raffineBordExterne(const int &nbMaillesY, const int &nbMaillesZ, const double &dXParent, const double &dYParent,
      const double &dZParent, Cellule *cellRef, const double &surfaceEnfant);      /*!< Raffinement du bord externe en creant si besoin des bords enfants + liaisons cellules/bords */
		void raffineBordExterneGhost(const int &nbMaillesY, const int &nbMaillesZ, const double &dXParent, const double &dYParent,
			const double &dZParent, Cellule *cellRef, const double &surfaceEnfant);      /*!< Raffinement du bord externe pour les cellules fantomes en creant si besoin des bords enfants + liaisons cellules/bords */
    virtual void deraffineBordExterne(Cellule *cellRef);        /*!< Deraffinement du bord externe en supprimant si besoin ses bords enfants + liaisons cellules/bords */
    void finaliseFace();                                        /*!< Supprime la face correspondante au bord */
    void deraffineBordsEnfants();                               /*!< Supprime les bords enfants */
    void constructionTableauBordsExternesLvl(std::vector<BordDeMaille *> *bordsLvl); /*!< Construction du nouveau tableau de bords du niveau (lvl + 1), bords externes ajoutes ici */
    bool getSplit() const;                                      /*!< Renvoie si oui ou non le bord est splitte */
    int getLvl() const;                                         /*!< Renvoie le niveau du bord */
    int getNombreBordsEnfants() const;                          /*!< Renvoie le nombre de bords enfants de ce bord*/
    BordDeMaille *getBordEnfant(const int &numEnfant);          /*!< Renvoie le bord enfant correspondant au numero */

   protected:
    Cellule *m_cellGauche;
    Cellule *m_cellDroite;
    Modele* m_mod;
    Face *m_face;
    
    //Attributs pour methode AMR
    int m_lvl;                                             /*!< Niveau dans l arbre AMR du bord */
    std::vector<BordDeMaille*> m_bordsEnfants;             /*!< Tableau de bords enfants (taille : 1 en 1D, 2 en 2D et 4 en 3D) */

  private:
};

//Utile pour la resolution des problemes de Riemann
extern Cellule *cellGauche;
extern Cellule *cellDroite;
extern Coord normale;
extern Coord tangente;
extern Coord binormale;

#endif // BORDDEMAILLE_H
