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

#include "BordDeMaille.h"
#include <iostream>

using namespace std;

//Utile pour la resolution des problemes de Riemann
Cellule *cellGauche;
Cellule *cellDroite;
Coord normale;
Coord tangente;
Coord binormale;

//***********************************************************************

BordDeMaille::BordDeMaille() : m_mod(0), m_cellGauche(0), m_cellDroite(0), m_face(0), m_bordsEnfants(0)
{
  m_lvl = 0;
}

//***********************************************************************

BordDeMaille::BordDeMaille(int lvl) : m_mod(0), m_cellGauche(0), m_cellDroite(0), m_face(0), m_bordsEnfants(0)
{
  m_lvl = lvl;
}

//***********************************************************************

BordDeMaille::~BordDeMaille()
{
  for (unsigned int i = 0; i < m_bordsEnfants.size(); i++) {
    m_bordsEnfants[i]->finaliseFace();
    delete m_bordsEnfants[i];
  }
  m_bordsEnfants.clear();
}

//***********************************************************************

void BordDeMaille::initialise(Cellule *cellGauche, Cellule *cellDroite)
{
  m_cellGauche = cellGauche;
  m_cellDroite = cellDroite;
}

//***********************************************************************

void BordDeMaille::initialiseGauche(Cellule *cellGauche)
{
  m_cellGauche = cellGauche;
}

//***********************************************************************

void BordDeMaille::initialiseDroite(Cellule *cellDroite)
{
  m_cellDroite = cellDroite;
}

//***********************************************************************

void BordDeMaille::setFace(Face *face)
{
  m_face = face;
}

//***********************************************************************

void BordDeMaille::calculFlux(const int &nombrePhases, const int &nombreTransports, double &dtMax, Limiteur &limiteurGlobal, Limiteur &limiteurInterface, Prim type)
{
  this->resolRiemann(nombrePhases, nombreTransports, dtMax, limiteurGlobal, limiteurInterface, type);

  if (m_cellGauche->getLvl() == m_cellDroite->getLvl()) {     //CoefAMR = 1 pour les deux
    this->ajoutFlux(nombrePhases, nombreTransports, 1.);      //Ajout du flux sur maille droite
    this->retireFlux(nombrePhases, nombreTransports, 1.);     //Retrait du flux sur maille gauche
  }
  else if (m_cellGauche->getLvl() > m_cellDroite->getLvl()) { //CoefAMR = 1 pour la gauche et 0.5 pour la droite
    this->ajoutFlux(nombrePhases, nombreTransports, 0.5);     //Ajout du flux sur maille droite
    this->retireFlux(nombrePhases, nombreTransports, 1.);     //Retrait du flux sur maille gauche
  }
  else {                                                      //CoefAMR = 0.5 pour la gauche et 1 pour la droite
    this->ajoutFlux(nombrePhases, nombreTransports, 1.);      //Ajout du flux sur maille droite
    this->retireFlux(nombrePhases, nombreTransports, 0.5);    //Retrait du flux sur maille gauche
  }
}

//***********************************************************************

void BordDeMaille::calculFluxPhysAdd(const int &nombrePhases, PhysAdd &physAdd)
{
  physAdd.calculFluxPhysAdd(this, nombrePhases);
}

//***********************************************************************

void BordDeMaille::resolRiemann(const int &nombrePhases, const int &nombreTransports, double &dtMax, Limiteur &limiteurGlobal, Limiteur &limiteurInterface, Prim type)
{
  cellGauche->copieVec(m_cellGauche->getPhases(), m_cellGauche->getMelange(), m_cellGauche->getTransports());
  cellDroite->copieVec(m_cellDroite->getPhases(), m_cellDroite->getMelange(), m_cellDroite->getTransports());

  //Projection des vitesses sur repere attache a la face
  normale = m_face->getNormale();
  tangente = m_face->getTangente();
  binormale = m_face->getBinormale();
  cellGauche->projection(normale, tangente, binormale, nombrePhases);
  cellDroite->projection(normale, tangente, binormale, nombrePhases);

  //Calcul des variables etendus (Phases, Melange, PhysAdd)
  cellGauche->calculsEtendusPourRiemann(nombrePhases);
  cellDroite->calculsEtendusPourRiemann(nombrePhases);

  //Probleme de Riemann
  double dxGauche(m_cellGauche->getElement()->getLCFL());
  double dxDroite(m_cellDroite->getElement()->getLCFL());
  dxGauche = dxGauche*pow(2., (double)m_lvl);
  dxDroite = dxDroite*pow(2., (double)m_lvl);
  m_mod->resolRiemannInterne(*cellGauche, *cellDroite, nombrePhases, dxGauche, dxDroite, dtMax);
  //Traitement des fonctions de transport (m_Sm connu : doit etre place apres l appel au Solveur de Riemann)
  if (nombreTransports > 0) { m_mod->resolRiemannTransportInterne(*cellGauche, *cellDroite, nombreTransports); }

  //Projection du flux sur le repere absolu
  m_mod->projectionRepereAbsolu(normale, tangente, binormale);
}

//***********************************************************************

void BordDeMaille::ajoutFlux(const int &nombrePhases, const int &nombreTransports, const double &coefAMR)
{
  double volume(m_cellDroite->getElement()->getVolume());
  double surface(m_face->getSurface());
  double coefA(surface / volume); //pas de "pas de temps"
  coefA = coefA*coefAMR;
  m_cellDroite->getCons()->ajoutFlux(coefA, nombrePhases);
  m_cellDroite->getCons()->ajoutNonCons(coefA, m_cellDroite, nombrePhases);
  double sM(m_mod->getSM());
  for (int k = 0; k < nombreTransports; k++) {
    m_cellDroite->getConsTransport(k)->ajoutFlux(coefA, k);
    m_cellDroite->getConsTransport(k)->ajoutNonCons(coefA, m_cellDroite->getTransport(k).getValeur(), sM);
  }
}

//***********************************************************************

void BordDeMaille::retireFlux(const int &nombrePhases, const int &nombreTransports, const double &coefAMR)
{
  double volume(m_cellGauche->getElement()->getVolume());
  double surface(m_face->getSurface());
  double coefA(surface / volume); //pas de "pas de temps"
  coefA = coefA*coefAMR;
  m_cellGauche->getCons()->retireFlux(coefA, nombrePhases);
  m_cellGauche->getCons()->retireNonCons(coefA, m_cellGauche, nombrePhases);
  double sM(m_mod->getSM());
  for (int k = 0; k < nombreTransports; k++) {
    m_cellGauche->getConsTransport(k)->retireFlux(coefA, k);
    m_cellGauche->getConsTransport(k)->retireNonCons(coefA, m_cellGauche->getTransport(k).getValeur(), sM);
  }
}

//***********************************************************************

double BordDeMaille::distance(Cellule *c)
{
  return m_face->distance(c->getElement());
}

//***********************************************************************
void BordDeMaille::EffetsSurface1D(const int &nombrePhases)
{
  Coord normale = m_face->getNormale();
  double surface = m_face->getSurface();
  if (m_cellDroite != NULL){ m_cellDroite->getCons()->ajoutTuyere1D(normale, surface, m_cellDroite, nombrePhases); }
  m_cellGauche->getCons()->retireTuyere1D(normale, surface, m_cellGauche, nombrePhases);
}

//***********************************************************************

void BordDeMaille::associeModele(Modele *modele)
{
  m_mod = modele;
}

//***********************************************************************

Face *BordDeMaille::getFace()
{
  return m_face;
}

//***********************************************************************

Modele *BordDeMaille::getMod() const
{
  return m_mod;
}

//***********************************************************************

Cellule *BordDeMaille::getCellGauche() const
{
  return m_cellGauche;
}

//***********************************************************************

Cellule *BordDeMaille::getCellDroite() const
{
  return m_cellDroite;
}

//****************************************************************************
//******************************Methode AMR***********************************
//****************************************************************************

void BordDeMaille::calculXi(const double &critereVar, const bool &varRho, const bool &varP, const bool &varU, const bool &varAlpha)
{
  if (varRho) { this->calculCritereAMR(critereVar, "RHO"); }
  if (varP) {
    if (m_cellGauche->getXi() < 0.99 || m_cellDroite->getXi() < 0.99) { this->calculCritereAMR(critereVar, "P"); }
  }
  if (varU) {
    if (m_cellGauche->getXi() < 0.99 || m_cellDroite->getXi() < 0.99) { this->calculCritereAMR(critereVar, "u"); }
  }
  if (varAlpha) {
    if (m_cellGauche->getXi() < 0.99 || m_cellDroite->getXi() < 0.99) { this->calculCritereAMR(critereVar, "ALPHA", 1); }
  }
}

//***********************************************************************

void BordDeMaille::calculCritereAMR(const double &critereVar, string nomVariable, int num)
{
  double valeurMin, variation, cd, cg;
  // Recuperation des valeurs de la variable en question a gauche et a droite
  cg = m_cellGauche->selectionneScalaire(nomVariable, num);
  cd = m_cellDroite->selectionneScalaire(nomVariable, num);

  // Valeur de la variation
  valeurMin = min(fabs(cd), fabs(cg));
	if (valeurMin < 1.e-4) { valeurMin = 1.e-4; } //Utile pour alpha (quasi-seulement) ou vitesse
  variation = fabs(cd - cg) / valeurMin;

  //Mise a jour de xi si la variation est superieure au critere
  if (variation >= critereVar) {
    m_cellGauche->setXi(1.);
    m_cellDroite->setXi(1.);
  }
}

//***********************************************************************

void BordDeMaille::calculFluxXi(const double &dXLocal, const double &dYLocal, const double &dZLocal)
{
  //Gradient de Xi sur la normale a la face
  double gradBord;
  double xig = m_cellGauche->getXi();
  double xid = m_cellDroite->getXi();
  double distance(m_cellGauche->distance(m_cellDroite));
  gradBord = (xid - xig) / distance;

  //Incrementation de consXi en fonction du bord en question
  //Ajout du flux sur maille droite
  //Retrait du flux sur maille gauche
  Coord normale = m_face->getNormale();
  double fluxXi(0.);
  if (fabs(normale.getX()) > 0.9) { fluxXi = - gradBord / dXLocal; }
  else if (fabs(normale.getY()) > 0.9) { fluxXi = - gradBord / dYLocal; }
  else { fluxXi = - gradBord / dZLocal; }
  m_cellDroite->ajoutFluxXi(fluxXi);
  m_cellGauche->retireFluxXi(fluxXi);
}

//***********************************************************************

void BordDeMaille::creerBordEnfant()
{
  m_bordsEnfants.push_back(new BordDeMaille(m_lvl + 1));
}

//***********************************************************************

void BordDeMaille::creerBordEnfantInterne(const int &lvl, vector<BordDeMaille*> *bordsEnfantsInternes)
{
  (*bordsEnfantsInternes).push_back(new BordDeMaille(lvl + 1));
}

//***********************************************************************

void BordDeMaille::creerFaceEnfant(BordDeMaille *bordParent)
{
  m_face = bordParent->m_face->creerNouvelleFace();
}

//***********************************************************************

void BordDeMaille::raffineBordExterne(const int &nbMaillesY, const int &nbMaillesZ, const double &dXParent, const double &dYParent,
  const double &dZParent, Cellule *cellRef, const double &surfaceEnfant)
{
  //La creation des bords enfants n'est pas systematique, on regarde d'abord si ces bords enfants ne sont pas deja crees.
  //Dans tous les cas on re-attribut les liaisons cellules/bords.

  double epsilon(1.e-6);
  int allouePenteLocal(1);

  if (nbMaillesZ == 1) {
    if (nbMaillesY == 1) {

      //--------------------------------------------------
      //--------------------- Cas 1D ---------------------
      //--------------------------------------------------

      //Bord pas encore split -> Creation des bords enfants
      //---------------------------------------------------
      if (m_lvl == cellRef->getLvl()) {

        this->creerBordEnfant();
        m_bordsEnfants[0]->m_face = m_face->creerNouvelleFace();
        m_bordsEnfants[0]->m_face->initialiseAutres(surfaceEnfant, m_face->getNormale(), m_face->getTangente(), m_face->getBinormale());
        m_bordsEnfants[0]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY(), m_face->getPos().getZ());
        if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
          //Bord numero 1 (gauche)
          m_bordsEnfants[0]->initialiseGauche(m_cellGauche);
          m_bordsEnfants[0]->initialiseDroite(cellRef->getCelluleEnfant(0));
          m_cellGauche->ajouteBord(m_bordsEnfants[0]);
          cellRef->getCelluleEnfant(0)->ajouteBord(m_bordsEnfants[0]);
        }
        else {
          //Bord numero 2 (droite)
          m_bordsEnfants[0]->initialiseGauche(cellRef->getCelluleEnfant(1));
          m_bordsEnfants[0]->initialiseDroite(m_cellDroite);
          cellRef->getCelluleEnfant(1)->ajouteBord(m_bordsEnfants[0]);
          m_cellDroite->ajouteBord(m_bordsEnfants[0]);
        }
        m_bordsEnfants[0]->associeModele(m_mod);
        m_bordsEnfants[0]->allouePentes(cellRef->getNombrePhases(), cellRef->getNombreTransports(), allouePenteLocal);
      }

      //Bord deja split -> on met seulement a jour les liaisons cellules/bords
      //----------------------------------------------------------------------
      else {
        if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
          //Bord numero 1 (gauche)
          m_cellDroite = cellRef->getCelluleEnfant(0);
          cellRef->getCelluleEnfant(0)->ajouteBord(this);
        }
        else {
          //Bord numero 2 (droite)
          m_cellGauche = cellRef->getCelluleEnfant(1);
          cellRef->getCelluleEnfant(1)->ajouteBord(this);
        }
      }
    }
    else {

      //--------------------------------------------------
      //--------------------- Cas 2D ---------------------
      //--------------------------------------------------

      //Bord pas encore split -> Creation des bords enfants
      //---------------------------------------------------
      if (m_lvl == cellRef->getLvl()) {

        //Creation des bords et faces enfants avec premiere initialisation
        //----------------------------------------------------------------
        for (int i = 0; i < 2; i++) {
          this->creerBordEnfant();
          m_bordsEnfants[i]->m_face = m_face->creerNouvelleFace();
          m_bordsEnfants[i]->m_face->initialiseAutres(surfaceEnfant, m_face->getNormale(), m_face->getTangente(), m_face->getBinormale());
        }

        //Face selon X
        //------------
        if (fabs(m_face->getNormale().getX()) > epsilon) {
          //Cote gauche
          if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
            for (int i = 0; i < 2; i++) {
              //Premiere face
              if (i == 0) {
                m_bordsEnfants[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ());
                m_bordsEnfants[i]->initialiseGauche(m_cellGauche);
                m_bordsEnfants[i]->initialiseDroite(cellRef->getCelluleEnfant(0));
                m_cellGauche->ajouteBord(m_bordsEnfants[i]);
                cellRef->getCelluleEnfant(0)->ajouteBord(m_bordsEnfants[i]);
              }
              //Deuxieme face
              else {
                m_bordsEnfants[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ());
                m_bordsEnfants[i]->initialiseGauche(m_cellGauche);
                m_bordsEnfants[i]->initialiseDroite(cellRef->getCelluleEnfant(2));
                m_cellGauche->ajouteBord(m_bordsEnfants[i]);
                cellRef->getCelluleEnfant(2)->ajouteBord(m_bordsEnfants[i]);
              }
            }
          }
          //Cote droite
          else {
            for (int i = 0; i < 2; i++) {
              //Premiere face
              if (i == 0) {
                m_bordsEnfants[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ());
                m_bordsEnfants[i]->initialiseGauche(cellRef->getCelluleEnfant(1));
                m_bordsEnfants[i]->initialiseDroite(m_cellDroite);
                cellRef->getCelluleEnfant(1)->ajouteBord(m_bordsEnfants[i]);
                m_cellDroite->ajouteBord(m_bordsEnfants[i]);
              }
              //Deuxieme face
              else {
                m_bordsEnfants[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ());
                m_bordsEnfants[i]->initialiseGauche(cellRef->getCelluleEnfant(3));
                m_bordsEnfants[i]->initialiseDroite(m_cellDroite);
                cellRef->getCelluleEnfant(3)->ajouteBord(m_bordsEnfants[i]);
                m_cellDroite->ajouteBord(m_bordsEnfants[i]);
              }
            }
          }
        }

        //Face selon Y
        //------------
        else {
          //Cote bas
          if (m_face->getPos().getY() < cellRef->getElement()->getPosition().getY()) {
            for (int i = 0; i < 2; i++) {
              //Premiere face
              if (i == 0) {
                m_bordsEnfants[i]->m_face->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ());
                m_bordsEnfants[i]->initialiseGauche(m_cellGauche);
                m_bordsEnfants[i]->initialiseDroite(cellRef->getCelluleEnfant(0));
                m_cellGauche->ajouteBord(m_bordsEnfants[i]);
                cellRef->getCelluleEnfant(0)->ajouteBord(m_bordsEnfants[i]);
              }
              //Deuxieme face
              else {
                m_bordsEnfants[i]->m_face->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ());
                m_bordsEnfants[i]->initialiseGauche(m_cellGauche);
                m_bordsEnfants[i]->initialiseDroite(cellRef->getCelluleEnfant(1));
                m_cellGauche->ajouteBord(m_bordsEnfants[i]);
                cellRef->getCelluleEnfant(1)->ajouteBord(m_bordsEnfants[i]);
              }
            }
          }
          //Cote haut
          else {
            for (int i = 0; i < 2; i++) {
              //Premiere face
              if (i == 0) {
                m_bordsEnfants[i]->m_face->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ());
                m_bordsEnfants[i]->initialiseGauche(cellRef->getCelluleEnfant(2));
                m_bordsEnfants[i]->initialiseDroite(m_cellDroite);
                cellRef->getCelluleEnfant(2)->ajouteBord(m_bordsEnfants[i]);
                m_cellDroite->ajouteBord(m_bordsEnfants[i]);
              }
              //Deuxieme face
              else {
                m_bordsEnfants[i]->m_face->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ());
                m_bordsEnfants[i]->initialiseGauche(cellRef->getCelluleEnfant(3));
                m_bordsEnfants[i]->initialiseDroite(m_cellDroite);
                cellRef->getCelluleEnfant(3)->ajouteBord(m_bordsEnfants[i]);
                m_cellDroite->ajouteBord(m_bordsEnfants[i]);
              }
            }
          }
        }

        //Association du modele et des pentes
        //-----------------------------------
        for (int i = 0; i < 2; i++) {
          m_bordsEnfants[i]->associeModele(m_mod);
          m_bordsEnfants[i]->allouePentes(cellRef->getNombrePhases(), cellRef->getNombreTransports(), allouePenteLocal);
        }

      }

      //Bord deja split -> on met seulement a jour les liaisons cellules/bords
      //----------------------------------------------------------------------
      else {

        //Face selon X
        //------------
        if (fabs(m_face->getNormale().getX()) > epsilon) {
          //Cote gauche
          if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
            //Premiere face
            if (m_face->getPos().getY() < cellRef->getElement()->getPosition().getY()) {
              m_cellDroite = cellRef->getCelluleEnfant(0);
              cellRef->getCelluleEnfant(0)->ajouteBord(this);
            }
            //Deuxieme face
            else {
              m_cellDroite = cellRef->getCelluleEnfant(2);
              cellRef->getCelluleEnfant(2)->ajouteBord(this);
            }
          }
          //Cote droite
          else {
            //Premiere face
            if (m_face->getPos().getY() < cellRef->getElement()->getPosition().getY()) {
              m_cellGauche = cellRef->getCelluleEnfant(1);
              cellRef->getCelluleEnfant(1)->ajouteBord(this);
            }
            //Deuxieme face
            else {
              m_cellGauche = cellRef->getCelluleEnfant(3);
              cellRef->getCelluleEnfant(3)->ajouteBord(this);
            }
          }
        }

        //Face selon Y
        //------------
        else {
          //Cote bas
          if (m_face->getPos().getY() < cellRef->getElement()->getPosition().getY()) {
            //Premiere face
            if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
              m_cellDroite = cellRef->getCelluleEnfant(0);
              cellRef->getCelluleEnfant(0)->ajouteBord(this);
            }
            //Deuxieme face
            else {
              m_cellDroite = cellRef->getCelluleEnfant(1);
              cellRef->getCelluleEnfant(1)->ajouteBord(this);
            }
          }
          //Cote haut
          else {
            //Premiere face
            if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
              m_cellGauche = cellRef->getCelluleEnfant(2);
              cellRef->getCelluleEnfant(2)->ajouteBord(this);
            }
            //Deuxieme face
            else {
              m_cellGauche = cellRef->getCelluleEnfant(3);
              cellRef->getCelluleEnfant(3)->ajouteBord(this);
            }
          }
        }
      }
    }
  }
  else {

    //--------------------------------------------------
    //--------------------- Cas 3D ---------------------
    //--------------------------------------------------

    //Bord pas encore split -> Creation des bords enfants
    //---------------------------------------------------
    if (m_lvl == cellRef->getLvl()) {

      //Creation des bords et faces enfants avec premiere initialisation
      //----------------------------------------------------------------
      for (int i = 0; i < 4; i++) {
        this->creerBordEnfant();
        m_bordsEnfants[i]->m_face = m_face->creerNouvelleFace();
        m_bordsEnfants[i]->m_face->initialiseAutres(surfaceEnfant, m_face->getNormale(), m_face->getTangente(), m_face->getBinormale());
      }

      //Face selon X
      //------------
      if (fabs(m_face->getNormale().getX()) > epsilon) {
        //Cote gauche
        if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
          for (int i = 0; i < 4; i++) {
            //Premiere face
            if (i == 0) {
              m_bordsEnfants[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ() - 0.25*dZParent);
              m_bordsEnfants[i]->initialiseGauche(m_cellGauche);
              m_bordsEnfants[i]->initialiseDroite(cellRef->getCelluleEnfant(0));
              m_cellGauche->ajouteBord(m_bordsEnfants[i]);
              cellRef->getCelluleEnfant(0)->ajouteBord(m_bordsEnfants[i]);
            }
            //Deuxieme face
            else if (i == 1) {
              m_bordsEnfants[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ() + 0.25*dZParent);
              m_bordsEnfants[i]->initialiseGauche(m_cellGauche);
              m_bordsEnfants[i]->initialiseDroite(cellRef->getCelluleEnfant(4));
              m_cellGauche->ajouteBord(m_bordsEnfants[i]);
              cellRef->getCelluleEnfant(4)->ajouteBord(m_bordsEnfants[i]);
            }
            //Troisieme face
            else if (i == 2) {
              m_bordsEnfants[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ() - 0.25*dZParent);
              m_bordsEnfants[i]->initialiseGauche(m_cellGauche);
              m_bordsEnfants[i]->initialiseDroite(cellRef->getCelluleEnfant(2));
              m_cellGauche->ajouteBord(m_bordsEnfants[i]);
              cellRef->getCelluleEnfant(2)->ajouteBord(m_bordsEnfants[i]);
            }
            //Quatrieme face
            else {
              m_bordsEnfants[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ() + 0.25*dZParent);
              m_bordsEnfants[i]->initialiseGauche(m_cellGauche);
              m_bordsEnfants[i]->initialiseDroite(cellRef->getCelluleEnfant(6));
              m_cellGauche->ajouteBord(m_bordsEnfants[i]);
              cellRef->getCelluleEnfant(6)->ajouteBord(m_bordsEnfants[i]);
            }
          }
        }
        //Cote droite
        else {
          for (int i = 0; i < 4; i++) {
            //Premiere face
            if (i == 0) {
              m_bordsEnfants[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ() - 0.25*dZParent);
              m_bordsEnfants[i]->initialiseGauche(cellRef->getCelluleEnfant(1));
              m_bordsEnfants[i]->initialiseDroite(m_cellDroite);
              cellRef->getCelluleEnfant(1)->ajouteBord(m_bordsEnfants[i]);
              m_cellDroite->ajouteBord(m_bordsEnfants[i]);
            }
            //Deuxieme face
            else if (i == 1) {
              m_bordsEnfants[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ() + 0.25*dZParent);
              m_bordsEnfants[i]->initialiseGauche(cellRef->getCelluleEnfant(5));
              m_bordsEnfants[i]->initialiseDroite(m_cellDroite);
              cellRef->getCelluleEnfant(5)->ajouteBord(m_bordsEnfants[i]);
              m_cellDroite->ajouteBord(m_bordsEnfants[i]);
            }
            //Troisieme face
            else if (i == 2) {
              m_bordsEnfants[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ() - 0.25*dZParent);
              m_bordsEnfants[i]->initialiseGauche(cellRef->getCelluleEnfant(3));
              m_bordsEnfants[i]->initialiseDroite(m_cellDroite);
              cellRef->getCelluleEnfant(3)->ajouteBord(m_bordsEnfants[i]);
              m_cellDroite->ajouteBord(m_bordsEnfants[i]);
            }
            //Quatrieme face
            else {
              m_bordsEnfants[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ() + 0.25*dZParent);
              m_bordsEnfants[i]->initialiseGauche(cellRef->getCelluleEnfant(7));
              m_bordsEnfants[i]->initialiseDroite(m_cellDroite);
              cellRef->getCelluleEnfant(7)->ajouteBord(m_bordsEnfants[i]);
              m_cellDroite->ajouteBord(m_bordsEnfants[i]);
            }
          }
        }
      }

      //Face selon Y
      //------------
      else if (fabs(m_face->getNormale().getY()) > epsilon) {
        //Cote bas
        if (m_face->getPos().getY() < cellRef->getElement()->getPosition().getY()) {
          for (int i = 0; i < 4; i++) {
            //Premiere face
            if (i == 0) {
              m_bordsEnfants[i]->m_face->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() - 0.25*dZParent);
              m_bordsEnfants[i]->initialiseGauche(m_cellGauche);
              m_bordsEnfants[i]->initialiseDroite(cellRef->getCelluleEnfant(0));
              m_cellGauche->ajouteBord(m_bordsEnfants[i]);
              cellRef->getCelluleEnfant(0)->ajouteBord(m_bordsEnfants[i]);
            }
            //Deuxieme face
            else if (i == 1) {
              m_bordsEnfants[i]->m_face->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() - 0.25*dZParent);
              m_bordsEnfants[i]->initialiseGauche(m_cellGauche);
              m_bordsEnfants[i]->initialiseDroite(cellRef->getCelluleEnfant(1));
              m_cellGauche->ajouteBord(m_bordsEnfants[i]);
              cellRef->getCelluleEnfant(1)->ajouteBord(m_bordsEnfants[i]);
            }
            //Troisieme face
            else if (i == 2) {
              m_bordsEnfants[i]->m_face->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() + 0.25*dZParent);
              m_bordsEnfants[i]->initialiseGauche(m_cellGauche);
              m_bordsEnfants[i]->initialiseDroite(cellRef->getCelluleEnfant(4));
              m_cellGauche->ajouteBord(m_bordsEnfants[i]);
              cellRef->getCelluleEnfant(4)->ajouteBord(m_bordsEnfants[i]);
            }
            //Quatrieme face
            else {
              m_bordsEnfants[i]->m_face->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() + 0.25*dZParent);
              m_bordsEnfants[i]->initialiseGauche(m_cellGauche);
              m_bordsEnfants[i]->initialiseDroite(cellRef->getCelluleEnfant(5));
              m_cellGauche->ajouteBord(m_bordsEnfants[i]);
              cellRef->getCelluleEnfant(5)->ajouteBord(m_bordsEnfants[i]);
            }
          }
        }
        //Cote haut
        else {
          for (int i = 0; i < 4; i++) {
            //Premiere face
            if (i == 0) {
              m_bordsEnfants[i]->m_face->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() - 0.25*dZParent);
              m_bordsEnfants[i]->initialiseGauche(cellRef->getCelluleEnfant(2));
              m_bordsEnfants[i]->initialiseDroite(m_cellDroite);
              cellRef->getCelluleEnfant(2)->ajouteBord(m_bordsEnfants[i]);
              m_cellDroite->ajouteBord(m_bordsEnfants[i]);
            }
            //Deuxieme face
            else if (i == 1) {
              m_bordsEnfants[i]->m_face->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() - 0.25*dZParent);
              m_bordsEnfants[i]->initialiseGauche(cellRef->getCelluleEnfant(3));
              m_bordsEnfants[i]->initialiseDroite(m_cellDroite);
              cellRef->getCelluleEnfant(3)->ajouteBord(m_bordsEnfants[i]);
              m_cellDroite->ajouteBord(m_bordsEnfants[i]);
            }
            //Troisieme face
            else if (i == 2) {
              m_bordsEnfants[i]->m_face->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() + 0.25*dZParent);
              m_bordsEnfants[i]->initialiseGauche(cellRef->getCelluleEnfant(6));
              m_bordsEnfants[i]->initialiseDroite(m_cellDroite);
              cellRef->getCelluleEnfant(6)->ajouteBord(m_bordsEnfants[i]);
              m_cellDroite->ajouteBord(m_bordsEnfants[i]);
            }
            //Quatrieme face
            else {
              m_bordsEnfants[i]->m_face->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() + 0.25*dZParent);
              m_bordsEnfants[i]->initialiseGauche(cellRef->getCelluleEnfant(7));
              m_bordsEnfants[i]->initialiseDroite(m_cellDroite);
              cellRef->getCelluleEnfant(7)->ajouteBord(m_bordsEnfants[i]);
              m_cellDroite->ajouteBord(m_bordsEnfants[i]);
            }
          }
        }
      }

      //Face selon Z
      //------------
      else {
        //Cote devant
        if (m_face->getPos().getZ() < cellRef->getElement()->getPosition().getZ()) {
          for (int i = 0; i < 4; i++) {
            //Premiere face
            if (i == 0) {
              m_bordsEnfants[i]->m_face->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ());
              m_bordsEnfants[i]->initialiseGauche(m_cellGauche);
              m_bordsEnfants[i]->initialiseDroite(cellRef->getCelluleEnfant(0));
              m_cellGauche->ajouteBord(m_bordsEnfants[i]);
              cellRef->getCelluleEnfant(0)->ajouteBord(m_bordsEnfants[i]);
            }
            //Deuxieme face
            else if (i == 1) {
              m_bordsEnfants[i]->m_face->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ());
              m_bordsEnfants[i]->initialiseGauche(m_cellGauche);
              m_bordsEnfants[i]->initialiseDroite(cellRef->getCelluleEnfant(1));
              m_cellGauche->ajouteBord(m_bordsEnfants[i]);
              cellRef->getCelluleEnfant(1)->ajouteBord(m_bordsEnfants[i]);
            }
            //Troisieme face
            else if (i == 2) {
              m_bordsEnfants[i]->m_face->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ());
              m_bordsEnfants[i]->initialiseGauche(m_cellGauche);
              m_bordsEnfants[i]->initialiseDroite(cellRef->getCelluleEnfant(2));
              m_cellGauche->ajouteBord(m_bordsEnfants[i]);
              cellRef->getCelluleEnfant(2)->ajouteBord(m_bordsEnfants[i]);
            }
            //Quatrieme face
            else {
              m_bordsEnfants[i]->m_face->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ());
              m_bordsEnfants[i]->initialiseGauche(m_cellGauche);
              m_bordsEnfants[i]->initialiseDroite(cellRef->getCelluleEnfant(3));
              m_cellGauche->ajouteBord(m_bordsEnfants[i]);
              cellRef->getCelluleEnfant(3)->ajouteBord(m_bordsEnfants[i]);
            }
          }
        }
        //Cote derriere
        else {
          for (int i = 0; i < 4; i++) {
            //Premiere face
            if (i == 0) {
              m_bordsEnfants[i]->m_face->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ());
              m_bordsEnfants[i]->initialiseGauche(cellRef->getCelluleEnfant(4));
              m_bordsEnfants[i]->initialiseDroite(m_cellDroite);
              cellRef->getCelluleEnfant(4)->ajouteBord(m_bordsEnfants[i]);
              m_cellDroite->ajouteBord(m_bordsEnfants[i]);
            }
            //Deuxieme face
            else if (i == 1) {
              m_bordsEnfants[i]->m_face->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ());
              m_bordsEnfants[i]->initialiseGauche(cellRef->getCelluleEnfant(5));
              m_bordsEnfants[i]->initialiseDroite(m_cellDroite);
              cellRef->getCelluleEnfant(5)->ajouteBord(m_bordsEnfants[i]);
              m_cellDroite->ajouteBord(m_bordsEnfants[i]);
            }
            //Troisieme face
            else if (i == 2) {
              m_bordsEnfants[i]->m_face->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ());
              m_bordsEnfants[i]->initialiseGauche(cellRef->getCelluleEnfant(6));
              m_bordsEnfants[i]->initialiseDroite(m_cellDroite);
              cellRef->getCelluleEnfant(6)->ajouteBord(m_bordsEnfants[i]);
              m_cellDroite->ajouteBord(m_bordsEnfants[i]);
            }
            //Quatrieme face
            else {
              m_bordsEnfants[i]->m_face->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ());
              m_bordsEnfants[i]->initialiseGauche(cellRef->getCelluleEnfant(7));
              m_bordsEnfants[i]->initialiseDroite(m_cellDroite);
              cellRef->getCelluleEnfant(7)->ajouteBord(m_bordsEnfants[i]);
              m_cellDroite->ajouteBord(m_bordsEnfants[i]);
            }
          }
        }
      }

      //Association du modele et des pentes
      //-----------------------------------
      for (int i = 0; i < 4; i++) {
        m_bordsEnfants[i]->associeModele(m_mod);
        m_bordsEnfants[i]->allouePentes(cellRef->getNombrePhases(), cellRef->getNombreTransports(), allouePenteLocal);
      }

    }

    //Bord deja split -> on met seulement a jour les liaisons cellules/bords
    //----------------------------------------------------------------------
    else {

      //Face selon X
      //------------
      if (fabs(fabs(m_face->getNormale().getX()) - 1.) < epsilon) {
        //Cote gauche
        if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
          if (m_face->getPos().getY() < cellRef->getElement()->getPosition().getY()) {
            if (m_face->getPos().getZ() < cellRef->getElement()->getPosition().getZ()) {
              //Premiere face
              m_cellDroite = cellRef->getCelluleEnfant(0);
              cellRef->getCelluleEnfant(0)->ajouteBord(this);
            }
            else {
              //Deuxieme face
              m_cellDroite = cellRef->getCelluleEnfant(4);
              cellRef->getCelluleEnfant(4)->ajouteBord(this);
            }
          }
          else {
            if (m_face->getPos().getZ() < cellRef->getElement()->getPosition().getZ()) {
              //Troisieme face
              m_cellDroite = cellRef->getCelluleEnfant(2);
              cellRef->getCelluleEnfant(2)->ajouteBord(this);
            }
            else {
              //Quatrieme face
              m_cellDroite = cellRef->getCelluleEnfant(6);
              cellRef->getCelluleEnfant(6)->ajouteBord(this);
            }
          }
        }
        //Cote droite
        else {
          if (m_face->getPos().getY() < cellRef->getElement()->getPosition().getY()) {
            if (m_face->getPos().getZ() < cellRef->getElement()->getPosition().getZ()) {
              //Premiere face
              m_cellGauche = cellRef->getCelluleEnfant(1);
              cellRef->getCelluleEnfant(1)->ajouteBord(this);
            }
            else {
              //Deuxieme face
              m_cellGauche = cellRef->getCelluleEnfant(5);
              cellRef->getCelluleEnfant(5)->ajouteBord(this);
            }
          }
          else {
            if (m_face->getPos().getZ() < cellRef->getElement()->getPosition().getZ()) {
              //Troisieme face
              m_cellGauche = cellRef->getCelluleEnfant(3);
              cellRef->getCelluleEnfant(3)->ajouteBord(this);
            }
            else {
              //Quatrieme face
              m_cellGauche = cellRef->getCelluleEnfant(7);
              cellRef->getCelluleEnfant(7)->ajouteBord(this);
            }
          }
        }
      }

      //Face selon Y
      //------------
      else if (fabs(fabs(m_face->getNormale().getY()) - 1.) < epsilon) {
        //Cote bas
        if (m_face->getPos().getY() < cellRef->getElement()->getPosition().getY()) {
          if (m_face->getPos().getZ() < cellRef->getElement()->getPosition().getZ()) {
            if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
              //Premiere face
              m_cellDroite = cellRef->getCelluleEnfant(0);
              cellRef->getCelluleEnfant(0)->ajouteBord(this);
            }
            else {
              //Deuxieme face
              m_cellDroite = cellRef->getCelluleEnfant(1);
              cellRef->getCelluleEnfant(1)->ajouteBord(this);
            }
          }
          else {
            if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
              //Troisieme face
              m_cellDroite = cellRef->getCelluleEnfant(4);
              cellRef->getCelluleEnfant(4)->ajouteBord(this);
            }
            else {
              //Quatrieme face
              m_cellDroite = cellRef->getCelluleEnfant(5);
              cellRef->getCelluleEnfant(5)->ajouteBord(this);
            }
          }
        }
        //Cote haut
        else {
          if (m_face->getPos().getZ() < cellRef->getElement()->getPosition().getZ()) {
            if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
              //Premiere face
              m_cellGauche = cellRef->getCelluleEnfant(2);
              cellRef->getCelluleEnfant(2)->ajouteBord(this);
            }
            else {
              //Deuxieme face
              m_cellGauche = cellRef->getCelluleEnfant(3);
              cellRef->getCelluleEnfant(3)->ajouteBord(this);
            }
          }
          else {
            if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
              //Troisieme face
              m_cellGauche = cellRef->getCelluleEnfant(6);
              cellRef->getCelluleEnfant(6)->ajouteBord(this);
            }
            else {
              //Quatrieme face
              m_cellGauche = cellRef->getCelluleEnfant(7);
              cellRef->getCelluleEnfant(7)->ajouteBord(this);
            }
          }
        }
      }

      //Face selon Z
      //------------
      else {
        //Cote devant
        if (m_face->getPos().getZ() < cellRef->getElement()->getPosition().getZ()) {
          if (m_face->getPos().getY() < cellRef->getElement()->getPosition().getY()) {
            if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
              //Premiere face
              m_cellDroite = cellRef->getCelluleEnfant(0);
              cellRef->getCelluleEnfant(0)->ajouteBord(this);
            }
            else {
              //Deuxieme face
              m_cellDroite = cellRef->getCelluleEnfant(1);
              cellRef->getCelluleEnfant(1)->ajouteBord(this);
            }
          }
          else {
            if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
              //Troisieme face
              m_cellDroite = cellRef->getCelluleEnfant(2);
              cellRef->getCelluleEnfant(2)->ajouteBord(this);
            }
            else {
              //Quatrieme face
              m_cellDroite = cellRef->getCelluleEnfant(3);
              cellRef->getCelluleEnfant(3)->ajouteBord(this);
            }
          }
        }
        //Cote derriere
        else {
          if (m_face->getPos().getY() < cellRef->getElement()->getPosition().getY()) {
            if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
              //Premiere face
              m_cellGauche = cellRef->getCelluleEnfant(4);
              cellRef->getCelluleEnfant(4)->ajouteBord(this);
            }
            else {
              //Deuxieme face
              m_cellGauche = cellRef->getCelluleEnfant(5);
              cellRef->getCelluleEnfant(5)->ajouteBord(this);
            }
          }
          else {
            if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
              //Troisieme face
              m_cellGauche = cellRef->getCelluleEnfant(6);
              cellRef->getCelluleEnfant(6)->ajouteBord(this);
            }
            else {
              //Quatrieme face
              m_cellGauche = cellRef->getCelluleEnfant(7);
              cellRef->getCelluleEnfant(7)->ajouteBord(this);
            }
          }
        }
      }
    }
  }
}


//***********************************************************************

void BordDeMaille::raffineBordExterneGhost(const int &nbMaillesY, const int &nbMaillesZ, const double &dXParent, const double &dYParent,
	const double &dZParent, Cellule *cellRef, const double &surfaceEnfant)
{
	//La creation des bords enfants n'est pas systematique, on regarde d'abord si ces bords enfants ne sont pas deja crees.
	//Dans tous les cas on re-attribut les liaisons cellules/bords.

	int allouePenteLocal(1);

	if (nbMaillesZ == 1) {
		if (nbMaillesY == 1) {

			//--------------------------------------------------
			//--------------------- Cas 1D ---------------------
			//--------------------------------------------------

			//Bord pas encore split -> Creation des bords enfants
			//---------------------------------------------------
			if (m_lvl == cellRef->getLvl()) {

				this->creerBordEnfant();
				m_bordsEnfants[0]->m_face = m_face->creerNouvelleFace();
				m_bordsEnfants[0]->m_face->initialiseAutres(surfaceEnfant, m_face->getNormale(), m_face->getTangente(), m_face->getBinormale());
				m_bordsEnfants[0]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY(), m_face->getPos().getZ());
				if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
					//Bord numero 1 (gauche)
					m_bordsEnfants[0]->initialiseGauche(m_cellGauche);
					m_bordsEnfants[0]->initialiseDroite(cellRef->getCelluleEnfant(0));
					m_cellGauche->ajouteBord(m_bordsEnfants[0]);
					cellRef->getCelluleEnfant(0)->ajouteBord(m_bordsEnfants[0]);
				}
				else {
					//Bord numero 2 (droite)
					m_bordsEnfants[0]->initialiseGauche(cellRef->getCelluleEnfant(0));
					m_bordsEnfants[0]->initialiseDroite(m_cellDroite);
					cellRef->getCelluleEnfant(0)->ajouteBord(m_bordsEnfants[0]);
					m_cellDroite->ajouteBord(m_bordsEnfants[0]);
				}
				m_bordsEnfants[0]->associeModele(m_mod);
				m_bordsEnfants[0]->allouePentes(cellRef->getNombrePhases(), cellRef->getNombreTransports(), allouePenteLocal);
			}

			//Bord deja split -> on met seulement a jour les liaisons cellules/bords
			//----------------------------------------------------------------------
			else {
				if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
					//Bord numero 1 (gauche)
					m_cellDroite = cellRef->getCelluleEnfant(0);
					cellRef->getCelluleEnfant(0)->ajouteBord(this);
				}
				else {
					//Bord numero 2 (droite)
					m_cellGauche = cellRef->getCelluleEnfant(0);
					cellRef->getCelluleEnfant(0)->ajouteBord(this);
				}
			}
		}
		else {

			//--------------------------------------------------
			//--------------------- Cas 2D ---------------------
			//--------------------------------------------------

			//Bord pas encore split -> Creation des bords enfants
			//---------------------------------------------------
			if (m_lvl == cellRef->getLvl()) {

				//Creation des bords et faces enfants avec premiere initialisation
				//----------------------------------------------------------------
				for (int i = 0; i < 2; i++) {
					this->creerBordEnfant();
					m_bordsEnfants[i]->m_face = m_face->creerNouvelleFace();
					m_bordsEnfants[i]->m_face->initialiseAutres(surfaceEnfant, m_face->getNormale(), m_face->getTangente(), m_face->getBinormale());
				}

				//Face selon X obligatoire
				//------------------------
				//Cote gauche
				if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
					for (int i = 0; i < 2; i++) {
						//Premiere face
						if (i == 0) {
							m_bordsEnfants[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ());
							m_bordsEnfants[i]->initialiseGauche(m_cellGauche);
							m_bordsEnfants[i]->initialiseDroite(cellRef->getCelluleEnfant(0));
							m_cellGauche->ajouteBord(m_bordsEnfants[i]);
							cellRef->getCelluleEnfant(0)->ajouteBord(m_bordsEnfants[i]);
						}
						//Deuxieme face
						else {
							m_bordsEnfants[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ());
							m_bordsEnfants[i]->initialiseGauche(m_cellGauche);
							m_bordsEnfants[i]->initialiseDroite(cellRef->getCelluleEnfant(1));
							m_cellGauche->ajouteBord(m_bordsEnfants[i]);
							cellRef->getCelluleEnfant(1)->ajouteBord(m_bordsEnfants[i]);
						}
					}
				}
				//Cote droite
				else {
					for (int i = 0; i < 2; i++) {
						//Premiere face
						if (i == 0) {
							m_bordsEnfants[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ());
							m_bordsEnfants[i]->initialiseGauche(cellRef->getCelluleEnfant(0));
							m_bordsEnfants[i]->initialiseDroite(m_cellDroite);
							cellRef->getCelluleEnfant(0)->ajouteBord(m_bordsEnfants[i]);
							m_cellDroite->ajouteBord(m_bordsEnfants[i]);
						}
						//Deuxieme face
						else {
							m_bordsEnfants[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ());
							m_bordsEnfants[i]->initialiseGauche(cellRef->getCelluleEnfant(1));
							m_bordsEnfants[i]->initialiseDroite(m_cellDroite);
							cellRef->getCelluleEnfant(1)->ajouteBord(m_bordsEnfants[i]);
							m_cellDroite->ajouteBord(m_bordsEnfants[i]);
						}
					}
				}

				//Association du modele et des pentes
				//-----------------------------------
				for (int i = 0; i < 2; i++) {
					m_bordsEnfants[i]->associeModele(m_mod);
					m_bordsEnfants[i]->allouePentes(cellRef->getNombrePhases(), cellRef->getNombreTransports(), allouePenteLocal);
				}

			}

			//Bord deja split -> on met seulement a jour les liaisons cellules/bords
			//----------------------------------------------------------------------
			else {

				//Face selon X obligatoire
				//------------------------
				//Cote gauche
				if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
					//Premiere face
					if (m_face->getPos().getY() < cellRef->getElement()->getPosition().getY()) {
						m_cellDroite = cellRef->getCelluleEnfant(0);
						cellRef->getCelluleEnfant(0)->ajouteBord(this);
					}
					//Deuxieme face
					else {
						m_cellDroite = cellRef->getCelluleEnfant(1);
						cellRef->getCelluleEnfant(1)->ajouteBord(this);
					}
				}
				//Cote droite
				else {
					//Premiere face
					if (m_face->getPos().getY() < cellRef->getElement()->getPosition().getY()) {
						m_cellGauche = cellRef->getCelluleEnfant(0);
						cellRef->getCelluleEnfant(0)->ajouteBord(this);
					}
					//Deuxieme face
					else {
						m_cellGauche = cellRef->getCelluleEnfant(1);
						cellRef->getCelluleEnfant(1)->ajouteBord(this);
					}
				}
			}
		}
	}
	else {

		//--------------------------------------------------
		//--------------------- Cas 3D ---------------------
		//--------------------------------------------------

		//Bord pas encore split -> Creation des bords enfants
		//---------------------------------------------------
		if (m_lvl == cellRef->getLvl()) {

			//Creation des bords et faces enfants avec premiere initialisation
			//----------------------------------------------------------------
			for (int i = 0; i < 4; i++) {
				this->creerBordEnfant();
				m_bordsEnfants[i]->m_face = m_face->creerNouvelleFace();
				m_bordsEnfants[i]->m_face->initialiseAutres(surfaceEnfant, m_face->getNormale(), m_face->getTangente(), m_face->getBinormale());
			}

			//Face selon X obligatoire
			//------------------------
			//Cote gauche
			if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
				for (int i = 0; i < 4; i++) {
					//Premiere face
					if (i == 0) {
						m_bordsEnfants[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ() - 0.25*dZParent);
						m_bordsEnfants[i]->initialiseGauche(m_cellGauche);
						m_bordsEnfants[i]->initialiseDroite(cellRef->getCelluleEnfant(0));
						m_cellGauche->ajouteBord(m_bordsEnfants[i]);
						cellRef->getCelluleEnfant(0)->ajouteBord(m_bordsEnfants[i]);
					}
					//Deuxieme face
					else if (i == 1) {
						m_bordsEnfants[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ() + 0.25*dZParent);
						m_bordsEnfants[i]->initialiseGauche(m_cellGauche);
						m_bordsEnfants[i]->initialiseDroite(cellRef->getCelluleEnfant(2));
						m_cellGauche->ajouteBord(m_bordsEnfants[i]);
						cellRef->getCelluleEnfant(2)->ajouteBord(m_bordsEnfants[i]);
					}
					//Troisieme face
					else if (i == 2) {
						m_bordsEnfants[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ() - 0.25*dZParent);
						m_bordsEnfants[i]->initialiseGauche(m_cellGauche);
						m_bordsEnfants[i]->initialiseDroite(cellRef->getCelluleEnfant(1));
						m_cellGauche->ajouteBord(m_bordsEnfants[i]);
						cellRef->getCelluleEnfant(1)->ajouteBord(m_bordsEnfants[i]);
					}
					//Quatrieme face
					else {
						m_bordsEnfants[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ() + 0.25*dZParent);
						m_bordsEnfants[i]->initialiseGauche(m_cellGauche);
						m_bordsEnfants[i]->initialiseDroite(cellRef->getCelluleEnfant(3));
						m_cellGauche->ajouteBord(m_bordsEnfants[i]);
						cellRef->getCelluleEnfant(3)->ajouteBord(m_bordsEnfants[i]);
					}
				}
			}
			//Cote droite
			else {
				for (int i = 0; i < 4; i++) {
					//Premiere face
					if (i == 0) {
						m_bordsEnfants[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ() - 0.25*dZParent);
						m_bordsEnfants[i]->initialiseGauche(cellRef->getCelluleEnfant(0));
						m_bordsEnfants[i]->initialiseDroite(m_cellDroite);
						cellRef->getCelluleEnfant(0)->ajouteBord(m_bordsEnfants[i]);
						m_cellDroite->ajouteBord(m_bordsEnfants[i]);
					}
					//Deuxieme face
					else if (i == 1) {
						m_bordsEnfants[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ() + 0.25*dZParent);
						m_bordsEnfants[i]->initialiseGauche(cellRef->getCelluleEnfant(2));
						m_bordsEnfants[i]->initialiseDroite(m_cellDroite);
						cellRef->getCelluleEnfant(2)->ajouteBord(m_bordsEnfants[i]);
						m_cellDroite->ajouteBord(m_bordsEnfants[i]);
					}
					//Troisieme face
					else if (i == 2) {
						m_bordsEnfants[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ() - 0.25*dZParent);
						m_bordsEnfants[i]->initialiseGauche(cellRef->getCelluleEnfant(1));
						m_bordsEnfants[i]->initialiseDroite(m_cellDroite);
						cellRef->getCelluleEnfant(1)->ajouteBord(m_bordsEnfants[i]);
						m_cellDroite->ajouteBord(m_bordsEnfants[i]);
					}
					//Quatrieme face
					else {
						m_bordsEnfants[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ() + 0.25*dZParent);
						m_bordsEnfants[i]->initialiseGauche(cellRef->getCelluleEnfant(3));
						m_bordsEnfants[i]->initialiseDroite(m_cellDroite);
						cellRef->getCelluleEnfant(3)->ajouteBord(m_bordsEnfants[i]);
						m_cellDroite->ajouteBord(m_bordsEnfants[i]);
					}
				}
			}

			//Association du modele et des pentes
			//-----------------------------------
			for (int i = 0; i < 4; i++) {
				m_bordsEnfants[i]->associeModele(m_mod);
				m_bordsEnfants[i]->allouePentes(cellRef->getNombrePhases(), cellRef->getNombreTransports(), allouePenteLocal);
			}

		}

		//Bord deja split -> on met seulement a jour les liaisons cellules/bords
		//----------------------------------------------------------------------
		else {

			//Face selon X obligatoire
			//------------------------
			//Cote gauche
			if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
				if (m_face->getPos().getY() < cellRef->getElement()->getPosition().getY()) {
					if (m_face->getPos().getZ() < cellRef->getElement()->getPosition().getZ()) {
						//Premiere face
						m_cellDroite = cellRef->getCelluleEnfant(0);
						cellRef->getCelluleEnfant(0)->ajouteBord(this);
					}
					else {
						//Deuxieme face
						m_cellDroite = cellRef->getCelluleEnfant(2);
						cellRef->getCelluleEnfant(2)->ajouteBord(this);
					}
				}
				else {
					if (m_face->getPos().getZ() < cellRef->getElement()->getPosition().getZ()) {
						//Troisieme face
						m_cellDroite = cellRef->getCelluleEnfant(1);
						cellRef->getCelluleEnfant(1)->ajouteBord(this);
					}
					else {
						//Quatrieme face
						m_cellDroite = cellRef->getCelluleEnfant(3);
						cellRef->getCelluleEnfant(3)->ajouteBord(this);
					}
				}
			}
			//Cote droite
			else {
				if (m_face->getPos().getY() < cellRef->getElement()->getPosition().getY()) {
					if (m_face->getPos().getZ() < cellRef->getElement()->getPosition().getZ()) {
						//Premiere face
						m_cellGauche = cellRef->getCelluleEnfant(0);
						cellRef->getCelluleEnfant(0)->ajouteBord(this);
					}
					else {
						//Deuxieme face
						m_cellGauche = cellRef->getCelluleEnfant(2);
						cellRef->getCelluleEnfant(2)->ajouteBord(this);
					}
				}
				else {
					if (m_face->getPos().getZ() < cellRef->getElement()->getPosition().getZ()) {
						//Troisieme face
						m_cellGauche = cellRef->getCelluleEnfant(1);
						cellRef->getCelluleEnfant(1)->ajouteBord(this);
					}
					else {
						//Quatrieme face
						m_cellGauche = cellRef->getCelluleEnfant(3);
						cellRef->getCelluleEnfant(3)->ajouteBord(this);
					}
				}
			}
		}
	}
}

//***********************************************************************

void BordDeMaille::deraffineBordExterne(Cellule *cellRef)
{
  //On parcourt seulement les bords parents pour regarder si la cellule voisine de celle de reference a des enfants,
  //si oui (enfants), on ne peut pas deraffiner le bord, on reaffecte donc les liaisons cellules/bords des bords enfants,
  //si non (pas enfants), on peut deraffiner le bord et mettre a jour les cellules voisines.
  //Plus, si je suis un bord parent, qui a donc des enfants, mais que la cellule de reference ne les connait pas encore, on les ajoute a ses bords.

  //Parcours les bords parents
  if (cellRef->getLvl() == m_lvl) {
    //cellRef est la cellule gauche
    if (m_cellGauche == cellRef) {
      //Si la cellule voisine (droite) a des enfants, on reaffecte les liaisons cellules/bords des bords enfants
      if (m_cellDroite->getSplit()) {
        for (unsigned int bordEnfant = 0; bordEnfant < m_bordsEnfants.size(); bordEnfant++) {
          m_bordsEnfants[bordEnfant]->initialiseGauche(cellRef);
        }
      }
      //La cellule voisine (droite) n'a pas d'enfants, on deraffine le bord et met a jour les cellules voisines
      else {
        //Il faut aussi enlever ses bords enfants des bords de la cellule droite et de mes bords
        for (unsigned int bordEnfant = 0; bordEnfant < m_bordsEnfants.size(); bordEnfant++) {
          m_cellDroite->supprimeBord(m_bordsEnfants[bordEnfant]);
          cellRef->supprimeBord(m_bordsEnfants[bordEnfant]);
        }
        this->deraffineBordsEnfants();
      }
    }
    //cellRef est la cellule droite
    else {
      //Si la cellule voisine (gauche) a des enfants, on reaffecte les liaisons cellules/bords des bords enfants
      if (m_cellGauche->getSplit()) {
        for (unsigned int bordEnfant = 0; bordEnfant < m_bordsEnfants.size(); bordEnfant++) {
          m_bordsEnfants[bordEnfant]->initialiseDroite(cellRef);
        }
      }
      //La cellule voisine (gauche) n'a pas d'enfants, on deraffine le bord et met a jour les cellules voisines
      else {
        //Il faut aussi enlever ses bords enfants des bords de la cellule gauche et de mes bords
        for (unsigned int bordEnfant = 0; bordEnfant < m_bordsEnfants.size(); bordEnfant++) {
          m_cellGauche->supprimeBord(m_bordsEnfants[bordEnfant]);
          cellRef->supprimeBord(m_bordsEnfants[bordEnfant]);
        }
        this->deraffineBordsEnfants();
      }
    }

    //Si je suis un bord parent, qui a donc des enfants, mais que la cellule de reference ne les connait pas encore, on les ajoute a ses bords.
    for (unsigned int bordEnfant = 0; bordEnfant < m_bordsEnfants.size(); bordEnfant++) {
      bool ajoutEnfantAuxBordsCelluleRef(true);
      for (int i = 0; i < cellRef->getBordsSize(); i++) {
        if (cellRef->getBord(i) == m_bordsEnfants[bordEnfant]) { ajoutEnfantAuxBordsCelluleRef = false; break; }
      }
      if (ajoutEnfantAuxBordsCelluleRef) { cellRef->ajouteBord(m_bordsEnfants[bordEnfant]); }
    }
  }
}

//***********************************************************************

void BordDeMaille::finaliseFace()
{
  delete m_face;
}

//***********************************************************************

void BordDeMaille::deraffineBordsEnfants()
{
  for (unsigned int i = 0; i < m_bordsEnfants.size(); i++) {
    m_bordsEnfants[i]->finaliseFace();
    delete m_bordsEnfants[i];
  }
  m_bordsEnfants.clear();
}

//***********************************************************************

void BordDeMaille::constructionTableauBordsExternesLvl(vector<BordDeMaille *> *bordsLvl)
{
  for (unsigned int i = 0; i < m_bordsEnfants.size(); i++) {
    bordsLvl[m_lvl + 1].push_back(m_bordsEnfants[i]);
  }
}

//***********************************************************************

bool BordDeMaille::getSplit() const
{
  bool split = false;
  if (m_bordsEnfants.size() > 0) { split = true; }
  return split;
}

//***********************************************************************

int BordDeMaille::getLvl() const
{
  return m_lvl;
}

//***********************************************************************

int BordDeMaille::getNombreBordsEnfants() const
{
  return m_bordsEnfants.size();
}

//***********************************************************************

BordDeMaille* BordDeMaille::getBordEnfant(const int &numEnfant)
{
  return m_bordsEnfants[numEnfant];
}

//***********************************************************************