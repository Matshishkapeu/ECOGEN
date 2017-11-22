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

#include "CondLim.h"
#include <iostream>

using namespace std;

//***********************************************************************

CondLim::CondLim(){}

//***********************************************************************

CondLim::CondLim(int numPhysique) : m_numPhysique(numPhysique)
{}

//***********************************************************************

CondLim::CondLim(const CondLim &Source) : m_numPhysique(Source.m_numPhysique)
{}

//***********************************************************************

CondLim::~CondLim() {}

//***********************************************************************

void CondLim::initialise(Cellule *cellGauche, Cellule *cellDroite)
{
  m_cellGauche = cellGauche;
  m_cellDroite = NULL;
}

//***********************************************************************

void CondLim::calculFlux(const int &nombrePhases, const int &nombreTransports, double &dtMax, Limiteur &limiteurGlobal, Limiteur &limiteurInterface, Prim type)
{
  this->resolRiemann(nombrePhases, nombreTransports, dtMax, limiteurGlobal, limiteurInterface, type);
  this->retireFlux(nombrePhases, nombreTransports, 1.); //Retrait du flux sur maille gauche
}

//***********************************************************************

void CondLim::calculFluxPhysAdd(const int &nombrePhases, PhysAdd &physAdd)
{
  physAdd.calculFluxPhysAddLimite(this, nombrePhases);
}

//***********************************************************************

void CondLim::resolRiemann(const int &nombrePhases, const int &nombreTransports, double &dtMax, Limiteur &limiteurGlobal, Limiteur &limiteurInterface, Prim type)
{
  cellGauche->copieVec(m_cellGauche->getPhases(type), m_cellGauche->getMelange(type), m_cellGauche->getTransports(type));

  //Projection des vitesses sur repere attache a la face
  Coord normale = m_face->getNormale();
  Coord tangente = m_face->getTangente();
  Coord binormale = m_face->getBinormale();
  cellGauche->projection(normale, tangente, binormale, nombrePhases);

  //Calcul des variables etendus (Phases, Melange, PhysAdd)
  cellGauche->calculsEtendusPourRiemann(nombrePhases);

  //Probleme de Riemann
  double dxGauche(m_cellGauche->getElement()->getLCFL());
  dxGauche = dxGauche*pow(2., (double)m_lvl);
  this->resolRiemannLimite(*cellGauche, nombrePhases, dxGauche, dtMax);
  //Traitement des fonctions de transport (m_Sm connu : doit etre place apres l appel au Solveur de Riemann)
  if (nombreTransports > 0) { this->resolRiemannTransportLimite(*cellGauche, nombreTransports); }

  //Projection du flux sur le repere absolu
  m_mod->projectionRepereAbsolu(normale, tangente, binormale);
}

//****************************************************************************

int CondLim::getNumPhys() const
{
  return m_numPhysique;
}

//****************************************************************************
//******************************Methode AMR***********************************
//****************************************************************************

void CondLim::raffineBordExterne(const int &nbMaillesY, const int &nbMaillesZ, const double &dXParent, const double &dYParent,
  const double &dZParent, Cellule *cellRef, const double &surfaceEnfant)
{
  //Le bord est une CL -> Creation des bords enfants

  double epsilon(1.e-6);
  int allouePenteLocal = 1;

  if (nbMaillesZ == 1) {
    if (nbMaillesY == 1) {

      //--------------------------------------------------
      //--------------------- Cas 1D ---------------------
      //--------------------------------------------------

      this->creerBordEnfant();
      m_bordsEnfants[0]->creerFaceEnfant(this);
      m_bordsEnfants[0]->getFace()->initialiseAutres(surfaceEnfant, m_face->getNormale(), m_face->getTangente(), m_face->getBinormale());
      m_bordsEnfants[0]->getFace()->setPos(m_face->getPos().getX(), m_face->getPos().getY(), m_face->getPos().getZ());
      if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
        //Bord numero 1 (gauche)
        m_bordsEnfants[0]->initialiseGauche(cellRef->getCelluleEnfant(0));
        cellRef->getCelluleEnfant(0)->ajouteBord(m_bordsEnfants[0]);
      }
      else {
        //Bord numero 2 (droite)
        m_bordsEnfants[0]->initialiseGauche(cellRef->getCelluleEnfant(1));
        cellRef->getCelluleEnfant(1)->ajouteBord(m_bordsEnfants[0]);
      }
      m_bordsEnfants[0]->associeModele(m_mod);
      m_bordsEnfants[0]->allouePentes(cellRef->getNombrePhases(), cellRef->getNombreTransports(), allouePenteLocal);
    }
    else {

      //--------------------------------------------------
      //--------------------- Cas 2D ---------------------
      //--------------------------------------------------

      //Creation des bords et faces enfants avec premiere initialisation
      //----------------------------------------------------------------
      for (int i = 0; i < 2; i++) {
        this->creerBordEnfant();
        m_bordsEnfants[i]->creerFaceEnfant(this);
        m_bordsEnfants[i]->getFace()->initialiseAutres(surfaceEnfant, m_face->getNormale(), m_face->getTangente(), m_face->getBinormale());
      }

      //Face selon X
      //------------
      if (fabs(m_face->getNormale().getX()) > epsilon) {
        //Cote gauche
        if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
          for (int i = 0; i < 2; i++) {
            //Premiere face
            if (i == 0) {
              m_bordsEnfants[i]->getFace()->setPos(m_face->getPos().getX(), m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ());
              m_bordsEnfants[i]->initialiseGauche(cellRef->getCelluleEnfant(0));
              cellRef->getCelluleEnfant(0)->ajouteBord(m_bordsEnfants[i]);
            }
            //Deuxieme face
            else {
              m_bordsEnfants[i]->getFace()->setPos(m_face->getPos().getX(), m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ());
              m_bordsEnfants[i]->initialiseGauche(cellRef->getCelluleEnfant(2));
              cellRef->getCelluleEnfant(2)->ajouteBord(m_bordsEnfants[i]);
            }
          }
        }
        //Cote droite
        else {
          for (int i = 0; i < 2; i++) {
            //Premiere face
            if (i == 0) {
              m_bordsEnfants[i]->getFace()->setPos(m_face->getPos().getX(), m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ());
              m_bordsEnfants[i]->initialiseGauche(cellRef->getCelluleEnfant(1));
              cellRef->getCelluleEnfant(1)->ajouteBord(m_bordsEnfants[i]);
            }
            //Deuxieme face
            else {
              m_bordsEnfants[i]->getFace()->setPos(m_face->getPos().getX(), m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ());
              m_bordsEnfants[i]->initialiseGauche(cellRef->getCelluleEnfant(3));
              cellRef->getCelluleEnfant(3)->ajouteBord(m_bordsEnfants[i]);
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
              m_bordsEnfants[i]->getFace()->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ());
              m_bordsEnfants[i]->initialiseGauche(cellRef->getCelluleEnfant(0));
              cellRef->getCelluleEnfant(0)->ajouteBord(m_bordsEnfants[i]);
            }
            //Deuxieme face
            else {
              m_bordsEnfants[i]->getFace()->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ());
              m_bordsEnfants[i]->initialiseGauche(cellRef->getCelluleEnfant(1));
              cellRef->getCelluleEnfant(1)->ajouteBord(m_bordsEnfants[i]);
            }
          }
        }
        //Cote haut
        else {
          for (int i = 0; i < 2; i++) {
            //Premiere face
            if (i == 0) {
              m_bordsEnfants[i]->getFace()->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ());
              m_bordsEnfants[i]->initialiseGauche(cellRef->getCelluleEnfant(2));
              cellRef->getCelluleEnfant(2)->ajouteBord(m_bordsEnfants[i]);
            }
            //Deuxieme face
            else {
              m_bordsEnfants[i]->getFace()->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ());
              m_bordsEnfants[i]->initialiseGauche(cellRef->getCelluleEnfant(3));
              cellRef->getCelluleEnfant(3)->ajouteBord(m_bordsEnfants[i]);
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
  }
  else {

    //--------------------------------------------------
    //--------------------- Cas 3D ---------------------
    //--------------------------------------------------

    //Creation des bords et faces enfants avec premiere initialisation
    //----------------------------------------------------------------
    for (int i = 0; i < 4; i++) {
      this->creerBordEnfant();
      m_bordsEnfants[i]->creerFaceEnfant(this);
      m_bordsEnfants[i]->getFace()->initialiseAutres(surfaceEnfant, m_face->getNormale(), m_face->getTangente(), m_face->getBinormale());
    }

    //Face selon X
    //------------
    if (fabs(m_face->getNormale().getX()) > epsilon) {
      //Cote gauche
      if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
        for (int i = 0; i < 4; i++) {
          //Premiere face
          if (i == 0) {
            m_bordsEnfants[i]->getFace()->setPos(m_face->getPos().getX(), m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ() - 0.25*dZParent);
            m_bordsEnfants[i]->initialiseGauche(cellRef->getCelluleEnfant(0));
            cellRef->getCelluleEnfant(0)->ajouteBord(m_bordsEnfants[i]);
          }
          //Deuxieme face
          else if (i == 1) {
            m_bordsEnfants[i]->getFace()->setPos(m_face->getPos().getX(), m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ() + 0.25*dZParent);
            m_bordsEnfants[i]->initialiseGauche(cellRef->getCelluleEnfant(4));
            cellRef->getCelluleEnfant(4)->ajouteBord(m_bordsEnfants[i]);
          }
          //Troisieme face
          else if (i == 2) {
            m_bordsEnfants[i]->getFace()->setPos(m_face->getPos().getX(), m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ() - 0.25*dZParent);
            m_bordsEnfants[i]->initialiseGauche(cellRef->getCelluleEnfant(2));
            cellRef->getCelluleEnfant(2)->ajouteBord(m_bordsEnfants[i]);
          }
          //Quatrieme face
          else {
            m_bordsEnfants[i]->getFace()->setPos(m_face->getPos().getX(), m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ() + 0.25*dZParent);
            m_bordsEnfants[i]->initialiseGauche(cellRef->getCelluleEnfant(6));
            cellRef->getCelluleEnfant(6)->ajouteBord(m_bordsEnfants[i]);
          }
        }
      }
      //Cote droite
      else {
        for (int i = 0; i < 4; i++) {
          //Premiere face
          if (i == 0) {
            m_bordsEnfants[i]->getFace()->setPos(m_face->getPos().getX(), m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ() - 0.25*dZParent);
            m_bordsEnfants[i]->initialiseGauche(cellRef->getCelluleEnfant(1));
            cellRef->getCelluleEnfant(1)->ajouteBord(m_bordsEnfants[i]);
          }
          //Deuxieme face
          else if (i == 1) {
            m_bordsEnfants[i]->getFace()->setPos(m_face->getPos().getX(), m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ() + 0.25*dZParent);
            m_bordsEnfants[i]->initialiseGauche(cellRef->getCelluleEnfant(5));
            cellRef->getCelluleEnfant(5)->ajouteBord(m_bordsEnfants[i]);
          }
          //Troisieme face
          else if (i == 2) {
            m_bordsEnfants[i]->getFace()->setPos(m_face->getPos().getX(), m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ() - 0.25*dZParent);
            m_bordsEnfants[i]->initialiseGauche(cellRef->getCelluleEnfant(3));
            cellRef->getCelluleEnfant(3)->ajouteBord(m_bordsEnfants[i]);
          }
          //Quatrieme face
          else {
            m_bordsEnfants[i]->getFace()->setPos(m_face->getPos().getX(), m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ() + 0.25*dZParent);
            m_bordsEnfants[i]->initialiseGauche(cellRef->getCelluleEnfant(7));
            cellRef->getCelluleEnfant(7)->ajouteBord(m_bordsEnfants[i]);
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
            m_bordsEnfants[i]->getFace()->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() - 0.25*dZParent);
            m_bordsEnfants[i]->initialiseGauche(cellRef->getCelluleEnfant(0));
            cellRef->getCelluleEnfant(0)->ajouteBord(m_bordsEnfants[i]);
          }
          //Deuxieme face
          else if (i == 1) {
            m_bordsEnfants[i]->getFace()->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() - 0.25*dZParent);
            m_bordsEnfants[i]->initialiseGauche(cellRef->getCelluleEnfant(1));
            cellRef->getCelluleEnfant(1)->ajouteBord(m_bordsEnfants[i]);
          }
          //Troisieme face
          else if (i == 2) {
            m_bordsEnfants[i]->getFace()->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() + 0.25*dZParent);
            m_bordsEnfants[i]->initialiseGauche(cellRef->getCelluleEnfant(4));
            cellRef->getCelluleEnfant(4)->ajouteBord(m_bordsEnfants[i]);
          }
          //Quatrieme face
          else {
            m_bordsEnfants[i]->getFace()->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() + 0.25*dZParent);
            m_bordsEnfants[i]->initialiseGauche(cellRef->getCelluleEnfant(5));
            cellRef->getCelluleEnfant(5)->ajouteBord(m_bordsEnfants[i]);
          }
        }
      }
      //Cote haut
      else {
        for (int i = 0; i < 4; i++) {
          //Premiere face
          if (i == 0) {
            m_bordsEnfants[i]->getFace()->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() - 0.25*dZParent);
            m_bordsEnfants[i]->initialiseGauche(cellRef->getCelluleEnfant(2));
            cellRef->getCelluleEnfant(2)->ajouteBord(m_bordsEnfants[i]);
          }
          //Deuxieme face
          else if (i == 1) {
            m_bordsEnfants[i]->getFace()->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() - 0.25*dZParent);
            m_bordsEnfants[i]->initialiseGauche(cellRef->getCelluleEnfant(3));
            cellRef->getCelluleEnfant(3)->ajouteBord(m_bordsEnfants[i]);
          }
          //Troisieme face
          else if (i == 2) {
            m_bordsEnfants[i]->getFace()->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() + 0.25*dZParent);
            m_bordsEnfants[i]->initialiseGauche(cellRef->getCelluleEnfant(6));
            cellRef->getCelluleEnfant(6)->ajouteBord(m_bordsEnfants[i]);
          }
          //Quatrieme face
          else {
            m_bordsEnfants[i]->getFace()->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() + 0.25*dZParent);
            m_bordsEnfants[i]->initialiseGauche(cellRef->getCelluleEnfant(7));
            cellRef->getCelluleEnfant(7)->ajouteBord(m_bordsEnfants[i]);
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
            m_bordsEnfants[i]->getFace()->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ());
            m_bordsEnfants[i]->initialiseGauche(cellRef->getCelluleEnfant(0));
            cellRef->getCelluleEnfant(0)->ajouteBord(m_bordsEnfants[i]);
          }
          //Deuxieme face
          else if (i == 1) {
            m_bordsEnfants[i]->getFace()->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ());
            m_bordsEnfants[i]->initialiseGauche(cellRef->getCelluleEnfant(1));
            cellRef->getCelluleEnfant(1)->ajouteBord(m_bordsEnfants[i]);
          }
          //Troisieme face
          else if (i == 2) {
            m_bordsEnfants[i]->getFace()->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ());
            m_bordsEnfants[i]->initialiseGauche(cellRef->getCelluleEnfant(2));
            cellRef->getCelluleEnfant(2)->ajouteBord(m_bordsEnfants[i]);
          }
          //Quatrieme face
          else {
            m_bordsEnfants[i]->getFace()->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ());
            m_bordsEnfants[i]->initialiseGauche(cellRef->getCelluleEnfant(3));
            cellRef->getCelluleEnfant(3)->ajouteBord(m_bordsEnfants[i]);
          }
        }
      }
      //Cote derriere
      else {
        for (int i = 0; i < 4; i++) {
          //Premiere face
          if (i == 0) {
            m_bordsEnfants[i]->getFace()->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ());
            m_bordsEnfants[i]->initialiseGauche(cellRef->getCelluleEnfant(4));
            cellRef->getCelluleEnfant(4)->ajouteBord(m_bordsEnfants[i]);
          }
          //Deuxieme face
          else if (i == 1) {
            m_bordsEnfants[i]->getFace()->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ());
            m_bordsEnfants[i]->initialiseGauche(cellRef->getCelluleEnfant(5));
            cellRef->getCelluleEnfant(5)->ajouteBord(m_bordsEnfants[i]);
          }
          //Troisieme face
          else if (i == 2) {
            m_bordsEnfants[i]->getFace()->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ());
            m_bordsEnfants[i]->initialiseGauche(cellRef->getCelluleEnfant(6));
            cellRef->getCelluleEnfant(6)->ajouteBord(m_bordsEnfants[i]);
          }
          //Quatrieme face
          else {
            m_bordsEnfants[i]->getFace()->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ());
            m_bordsEnfants[i]->initialiseGauche(cellRef->getCelluleEnfant(7));
            cellRef->getCelluleEnfant(7)->ajouteBord(m_bordsEnfants[i]);
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
}

//***********************************************************************

void CondLim::deraffineBordExterne(Cellule *cellRef)
{
  //Dans le cas des CL, on parcourt toujours les bords parents mais on peut directement les deraffiner et mettre a jour la cellule gauche.

  //Parcours les bords parents
  if (cellRef->getLvl() == m_lvl) {
    //cellRef est forcement la cellule gauche, on deraffine le bord parent et on met a jour la cellule gauche (de reference)
    for (unsigned int bordEnfant = 0; bordEnfant < m_bordsEnfants.size(); bordEnfant++) {
      cellRef->supprimeBord(m_bordsEnfants[bordEnfant]);
    }
    this->deraffineBordsEnfants();
  }
}

//****************************************************************************