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

#ifndef CELLULE_H
#define CELLULE_H

#include <vector>
#include <fstream>
#include "Modeles/Phase.h"
#include "Maths/Coord.h"
#include "Transport/Transport.h"

class Cellule; //Predeclaration de la classe Cellule pour pouvoir inclure les .h ci-dessous

#include "Modeles/Melange.h"
#include "PhysiqueAdditionnelle/GrandeursPhysAdd.h"
#include "BordDeMaille.h"
#include "Modeles/Modele.h"
#include "Modeles/Flux.h"
#include "Maillages/Element.h"
#include "Geometries/DomaineGeometrique.h"


class Cellule
{
    public:
        Cellule();
        Cellule(int lvl); //Pour AMR
        virtual ~Cellule();
        void ajouteBord(BordDeMaille *bord);
        void supprimeBord(BordDeMaille *bord);
        int getBordsSize() const;
        BordDeMaille* getBord(int &b);
        virtual void alloue(const int &nombrePhases, const int &nombreTransports, const std::vector<PhysAdd*> &physAdd, Modele *modele);
        void alloueEosYk(const int &nombrePhases, Modele *modele);
        void rempli(std::vector<DomaineGeometrique*> &domaines);
        virtual void alloueEtCopiePhase(const int &numeroPhase, Phase *phase);
        virtual void copiePhase(const int &numeroPhase, Phase *phase);
        void copieMelange(Melange *melange);
        void miseAZeroCons(const int &nombrePhases, const int &nombreTransports);
        void miseAZeroConsGenerale(const int &nombrePhases, const int &nombreTransports);
        void miseAZeroFluxTemp(const int &nombrePhases);
        void evolutionTemporelle(const double &dt, const int &nombrePhases, const int &nombreTransports);
        void construitPrim(const int &nombrePhases);
        void construitCons(const int &nombrePhases);
        void relaxPressions(const int &nombrePhases);
        void relaxPTMu(const int &nombrePhases);
        void correctionEnergie(const int &nombrePhases);
        void integreTermeSource(const double &dt, const int &nombrePhases) {};
        void ecritPhasesMelange(const int &nombrePhases, const int &nombreTransports, std::ofstream &fluxFichier) const;
        virtual void calculsEtendus(const int &nombrePhases, Prim type = vecPhases);
        void calculsEtendusPourRiemann(const int &nombrePhases);
        virtual void calculsEtendusPourCommunications(const int &nombrePhases, Prim type = vecPhases);
        virtual void projection(const Coord &normale, const Coord &tangente, const Coord &binormale, const int &nombrePhases, Prim type = vecPhases);
        virtual void projectionRepereAbsolu(const Coord &normale, const Coord &tangente, const Coord &binormale, const int &nombrePhases, Prim type = vecPhases);
        virtual void copieDansCellule(Cellule &cellSource, Prim type=vecPhases) const { Erreurs::messageErreur("methode copie non dispo pour cellule"); };
        void copieVec(Phase **vecPhases, Melange *melange, Transport *vecTransports);
        //void ecritureCoupe1Dde2D(std::ofstream &fluxFichier, std::string variableConstanteCoupe, const double &valeurCoupe, const double &dL);    /*!< Ecriture de la coupe 1D de la simulation 2D  */
        //void ecritureCoupe1Dde3D(std::ofstream &fluxFichier, std::string variableConstanteCoupe1, std::string variableConstanteCoupe2, const double &valeurCoupe1, const double &valeurCoupe2, const double &dL1, const double &dL2);                                            /*!< Ecriture de la coupe 1D de la simulation 3D  */

        //Pour physiques additionnelles
        void preparePhysAdd();
        double selectionneScalaire(std::string nomVariable, int num=0) const;
        void setScalaire(std::string nomVariable, const double &valeur, int num = 0, int indice = -1);
        Coord selectionneVecteur(std::string nomVecteur, int num=0, int indice=-1) const;
        void setVecteur(std::string nomVecteur, const Coord &valeur, int num=0, int indice=-1);
        
        Coord calculGradient(std::string nomVariable, int num=-1);
       
        GrandeursPhysAdd* getGPA(int &numGPA) const; //!< Permet de recuperer une grandeur physique additionelle

        Coord getGradTk(int &numPhase, int &numPhysAdd) const;
        void setGradTk(int &numPhase, int &numPhysAdd, double *tampon, int &compteur);
        void ajoutNonConsPhysAdd(const int &nombrePhases, PhysAdd &physAdd);

				void reinitialiseFonctionCouleur(); //!< Re-initialise la fonction couleur (transport) avec alpha
        
        //Accesseurs
        virtual Phase* getPhase(const int &numeroPhase, Prim type=vecPhases) const;
        virtual Phase** getPhases(Prim type=vecPhases) const;
        virtual Melange* getMelange(Prim type = vecPhases) const;
        Flux* getCons() const;
        void setCons(Flux *cons);
        Coord getPosition() const;
        void setElement(Element *element, const int &numCellule);
        Element* getElement();
        virtual void setTransport(double valeur, int &numTransport, Prim type = vecPhases);
        virtual Transport& getTransport(const int &numTransport, Prim type = vecPhases) const;
        virtual Transport* getTransports(Prim type = vecPhases) const;
        Transport* getConsTransport(const int &numTransport) const;
        void setConsTransport(double valeur, const int &numTransport);
        std::vector<GrandeursPhysAdd*>& getVecGrandeursPhysAdd();
        int getNombrePhases() const;
        int getNombreTransports() const;
        double getXi() const { return m_xi; };
        double getGradient();

        //Inutilise pour cellules ordre 1
        virtual void calculPentesLocal(const int &nombrePhases, const int &nombreTransports, BordDeMaille &bord, Limiteur &limiteurGlobal, Limiteur &limiteurInterface) {};  /*!< Ne fait rien pour des cellules ordre 1 */
        virtual void calculPentesLocalLimite(const int &nombrePhases, const int &nombreTransports, BordDeMaille &bord, Limiteur &limiteurGlobal, Limiteur &limiteurInterface) {};  /*!< Ne fait rien pour des cellules ordre 1 */
        virtual void calculMultiPente(const int &nombrePhases, BordDeMaille *bord, Limiteur *limiteurGlobal) {};                      /*!< Ne fait rien pour des cellules ordre 1 */
        virtual Phase* getPentes(const int &numeroPhase) const { return 0; };                                                         /*!< Ne fait rien pour des cellules ordre 1 */
        virtual Transport* getPentesTransport(const int &numeroTransport) const { return 0; };                                        /*!< Ne fait rien pour des cellules ordre 1 */
        virtual void sauvegardeCons(const int &nombrePhases, const int &nombreTransports) {};                                         /*!< Ne fait rien pour des cellules ordre 1 */
        virtual void recuperationCons(const int &nombrePhases, const int &nombreTransports) {};                                       /*!< Ne fait rien pour des cellules ordre 1 */
        virtual void predictionOrdre2(const double &dt, const int &nombrePhases, const int &nombreTransports) {};                     /*!< Ne fait rien pour des cellules ordre 1 */

        //BordDeMaille *getBordDeMaille(); //FP//TODO// A FAIRE...

        void afficheInfos() const;

        //Distance de la cellule a une autre cellule ou un bord de maille
        double distance(Cellule *c);        /*!< Distance totale  */
        double distanceX(Cellule *c);       /*!< Distance selon x */
        double distanceY(Cellule *c);       /*!< Distance selon y */
        double distanceZ(Cellule *c);       /*!< Distance selon z */
        double distance(BordDeMaille *b);   /*!< Distance totale  */
        double distanceX(BordDeMaille *b);  /*!< Distance selon x */
        double distanceY(BordDeMaille *b);  /*!< Distance selon y */
        double distanceZ(BordDeMaille *b);  /*!< Distance selon z */

        bool traverseObjet(const ObjetGeometrique &objet) const;

        //Ecriture
        void ecritureGnuplotAMR(std::ofstream &fluxFichier, const int &dim, ObjetGeometrique *objet = 0);

        //Pour methode AMR
        void miseAZeroXi();                                              /*!< Mise a zero de m_consXi */
        void miseAZeroConsXi();                                          /*!< Mise a zero de m_consXi */
        void evolutionTemporelleXi(const double &dtDiff);                /*!< Evolution temporelle de Xi pour smoothing */
        void choixRaffine(const double &xiSplit, const int &nbMaillesY, const int &nbMaillesZ,
          const double &dX, const double &dY, const double &dZ, const std::vector<PhysAdd*> &physAdd, Modele *modele, int &nbMaillesTotalAMR); /*!< Choix du raffinement ou non de la cellule parent */
        void choixDeraffine(const double &xiJoin, int &nbMaillesTotalAMR);                                                                     /*!< Choix du deraffinement ou non de la cellule parent */
        void raffineCelluleEtBords(const int &nbMaillesY, const int &nbMaillesZ, const double &dX, const double &dY, const double &dZ,
          const std::vector<PhysAdd*> &physAdd, Modele *modele);         /*!< Raffinement de la cellule parent en creant des cellules enfants */
        virtual void creerCelluleEnfant(const int &num, const int &lvl); /*!< Creer une cellule enfant (non initialisee) */
        void deraffineCelluleEtBords();                                  /*!< Deraffinement de la cellule parent en supprimant les cellules enfants */
        void moyenneEnfantsDansParent();                                 /*!< Moyenne des variables des cellules enfants dans la cellule parent, utile pour le calcul de xi. */
        bool lvlVoisinTropGrand();                                       /*!< Regarde si un ou des voisins a/ont un niveau AMR trop eleve */
        bool lvlVoisinTropPetit();                                       /*!< Regarde si au moins un des voisins a un niveau AMR trop petit pour pouvoir raffiner */
        void constructionTableauxCellulesLvlEtBordsInternesLvl(std::vector<Cellule *> *cellulesLvl, std::vector<BordDeMaille *> *bordsLvl);   /*!< Construction des nouveaux tableaux de cellules et de bords du niveau (lvl + 1), juste bords internes ajoutes ici */
        int getLvl();                                                    /*!< Niveau dans l arbre AMR de la cellule */
        bool getSplit();                                                 /*!< Renvoie si oui ou non la cellule est splittee */
        double getXi();                                                  /*!< Renvoie la valeur de Xi de la cellule */
        void setXi(double valeur);                                       /*!< Ecrit la valeur de Xi de la cellule */
        void ajoutFluxXi(double valeur);                                 /*!< Ajout du flux de xi sur la cellule */
        void retireFluxXi(double valeur);                                /*!< Retrait du flux de xi sur la cellule */
        int getNombreCellulesEnfants();                                  /*!< Niveau dans l arbre AMR de la cellule */
        Cellule* getCelluleEnfant(const int &num);                       /*!< Renvoie la cellule enfant correspondant au numero donne */

        //Pour parallele non-AMR
        void rempliTamponPrimitives(double *tampon, int &compteur, Prim type = vecPhases) const;
        void recupereTamponPrimitives(double *tampon, int &compteur, Eos **eos, Prim type = vecPhases);
        virtual void rempliTamponPentes(double *tampon, int &compteur) const {};                                     /*!< Ne fait rien pour des cellules ordre 1 */
        virtual void recupereTamponPentes(double *tampon, int &compteur) {};                                         /*!< Ne fait rien pour des cellules ordre 1 */
        void rempliTamponScalaire(double *tampon, int &compteur, std::string nomVariable) const;
        void recupereTamponScalaire(double *tampon, int &compteur, std::string nomVariable);
        void rempliTamponVecteur(double *tampon, int &compteur, const int &dim, std::string nomVecteur, int num = 0, int indice = -1) const;
        void recupereTamponVecteur(double *tampon, int &compteur, const int &dim, std::string nomVecteur, int num = 0, int indice = -1);
        void rempliTamponTransports(double *tampon, int &compteur) const;
        void recupereTamponTransports(double *tampon, int &compteur);
        virtual bool estCelluleO2Ghost() const { return false; };

        //Pour parallele AMR
        void rempliTamponPrimitivesAMRjeSuisCpuGauche(double *tampon, int &compteur, const int &lvl, Prim type = vecPhases) const;
        void rempliTamponPrimitivesAMRjeSuisCpuDroite(double *tampon, int &compteur, const int &lvl, Prim type = vecPhases) const;
        void recupereTamponPrimitivesAMR(double *tampon, int &compteur, const int &lvl, Eos **eos, Prim type = vecPhases);
        virtual void rempliTamponPentesAMRjeSuisCpuGauche(double *tampon, int &compteur, const int &lvl) const {};   /*!< Ne fait rien pour des cellules ordre 1 */
        virtual void rempliTamponPentesAMRjeSuisCpuDroite(double *tampon, int &compteur, const int &lvl) const {};   /*!< Ne fait rien pour des cellules ordre 1 */
        virtual void recupereTamponPentesAMR(double *tampon, int &compteur, const int &lvl) {};                      /*!< Ne fait rien pour des cellules ordre 1 */
        void rempliTamponScalaireAMRjeSuisCpuGauche(double *tampon, int &compteur, const int &lvl, std::string nomVariable) const;
        void rempliTamponScalaireAMRjeSuisCpuDroite(double *tampon, int &compteur, const int &lvl, std::string nomVariable) const;
        void recupereTamponScalaireAMR(double *tampon, int &compteur, const int &lvl, std::string nomVariable);
        void rempliTamponVecteurAMRjeSuisCpuGauche(double *tampon, int &compteur, const int &lvl, const int &dim, std::string nomVecteur, int num = 0, int indice = -1) const;
        void rempliTamponVecteurAMRjeSuisCpuDroite(double *tampon, int &compteur, const int &lvl, const int &dim, std::string nomVecteur, int num = 0, int indice = -1) const;
        void recupereTamponVecteurAMR(double *tampon, int &compteur, const int &lvl, const int &dim, std::string nomVecteur, int num = 0, int indice = -1);
        void rempliTamponTransportsAMRjeSuisCpuGauche(double *tampon, int &compteur, const int &lvl) const;
        void rempliTamponTransportsAMRjeSuisCpuDroite(double *tampon, int &compteur, const int &lvl) const;
        void recupereTamponTransportsAMR(double *tampon, int &compteur, const int &lvl);
        void choixRaffineDeraffineGhost(const int &nbMaillesY, const int &nbMaillesZ, const double &dX, const double &dY, const double &dZ,
          const std::vector<PhysAdd*> &physAdd, Modele *modele, std::vector<Cellule *> *cellulesLvlGhost);           /*!< Choix du raffinement, deraffinement ou non de la cellule fantome parent + Mise a jour du tableau de cellules fantomes de niveau lvl + 1 */
        void raffineCelluleEtBordsGhost(const int &nbMaillesY, const int &nbMaillesZ, const double &dX, const double &dY, const double &dZ,
          const std::vector<PhysAdd*> &physAdd, Modele *modele);                            /*!< Raffinement de la cellule fantome parent en creant des cellules fantomes enfants */
        void deraffineCelluleEtBordsGhost();                                                /*!< Deraffinement de la cellule fantome parent en supprimant les cellules fantomes enfants */
        void rempliTamponXiJeSuisCpuGauche(double *tampon, int &compteur, const int &lvl) const;
        void rempliTamponXiJeSuisCpuDroite(double *tampon, int &compteur, const int &lvl) const;
        void recupereTamponXi(double *tampon, int &compteur, const int &lvl);
        void rempliTamponSplitJeSuisCpuGauche(bool *tampon, int &compteur, const int &lvl) const;
        void rempliTamponSplitJeSuisCpuDroite(bool *tampon, int &compteur, const int &lvl) const;
        void recupereTamponSplit(bool *tampon, int &compteur, const int &lvl);
        void rempliNombreElementsAEnvoyerAVoisinJeSuisCpuGauche(int &nombreNombreElementsAEnvoyerAVoisin, const int &lvl);
        void rempliNombreElementsAEnvoyerAVoisinJeSuisCpuDroite(int &nombreNombreElementsAEnvoyerAVoisin, const int &lvl);

    protected:
      int m_nombrePhases;
      int m_nombreTransports;
      Phase **m_vecPhases;
      Melange *m_melange;
      Transport *m_vecTransports;                            /*!< Vecteur de grandeurs passives transportees par l ecoulement */
      Flux *m_cons;                                          /*!< Vecteur de variables conservatives. De type flux car recueille la somme des flux sur l objet cellule */
      Transport *m_consTransports;                           /*!< Vecteur de grandeurs passives permettant de recueillir la somme des flux des grandeurs transportees */
      Element *m_element;                                    /*!< Pointeur vers element de maillage correspondant */
      std::vector<BordDeMaille*> m_bords;                    /*!< Tableau de Pointeurs vers les Bords de la Cellule */
      std::vector<GrandeursPhysAdd*> m_vecGrandeursPhysAdd;  /*!< Vecteur de Pointeurs vers les Grandeurs des Physiques Additionnelles de la Cellule */
     
      //Attributs pour methode AMR
      int m_lvl;                                             /*!< Niveau dans l arbre AMR de la cellule */
      double m_xi;                                           /*!< Variable servant de critere pour le raffinement et deraffinement des cellules */
      double m_consXi;                                       /*!< Variable tampon pour stocker la somme des flux de Xi */
			bool m_split;                                          /*!< Variable pour savoir si la cellule a des enfants (split) ou non */
      std::vector<Cellule*> m_cellulesEnfants;               /*!< Tableau de pointeurs vers les cellules enfants */
      std::vector<BordDeMaille*> m_bordsEnfantsInternes;     /*!< Tableau de pointeurs vers les bords enfants internes de la cellule */

    private:
};

#endif // CELLULE_H
