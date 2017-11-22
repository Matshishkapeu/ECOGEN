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

#include <iostream>
#include <algorithm>
#include <sstream>

#include "MaillageCartesienAMR.h"

using namespace std;

//***********************************************************************

MaillageCartesienAMR::MaillageCartesienAMR(double lX, int nombreMaillesX,
  double lY, int nombreMaillesY,
  double lZ, int nombreMaillesZ,
	int lvlMax, double critereVar, bool varRho, bool varP, bool varU, bool varAlpha, double xiSplit, double xiJoin) :
  MaillageCartesien(lX, nombreMaillesX, lY, nombreMaillesY, lZ, nombreMaillesZ),
  m_lvlMax(lvlMax), m_critereVar(critereVar), m_varRho(varRho), m_varP(varP), m_varU(varU), m_varAlpha(varAlpha), m_xiSplit(xiSplit), m_xiJoin(xiJoin)
{
  m_type = UNS;
}

//***********************************************************************

MaillageCartesienAMR::~MaillageCartesienAMR(){
  delete[] m_cellulesLvlGhost;
}

//***********************************************************************

void MaillageCartesienAMR::genereTableauxCellulesBordsLvl(Cellule **cellules, BordDeMaille **bord, vector<Cellule *> **cellulesLvl,
	vector<BordDeMaille *> **bordsLvl)
{
	(*cellulesLvl) = new vector<Cellule *>[m_lvlMax + 1];
	for (int i = 0; i < m_nombreCellulesCalcul; i++) { (*cellulesLvl)[0].push_back(cellules[i]); }

	(*bordsLvl) = new vector<BordDeMaille *>[m_lvlMax + 1];
	for (int i = 0; i < m_nombreFacesTotal; i++) { (*bordsLvl)[0].push_back(bord[i]); }

	m_cellules = cellulesLvl;

	if (Ncpu > 1) {
		//Genere les tableaux de cellules fantomes par niveau
		m_cellulesLvlGhost = new vector<Cellule *>[m_lvlMax + 1];
		for (int i = m_nombreCellulesCalcul; i < m_nombreCellulesTotal; i++) {
			m_cellulesLvlGhost[0].push_back(cellules[i]);
		}
	}
}

//***********************************************************************

void MaillageCartesienAMR::procedureRaffinementInitialisation(vector<Cellule *> *cellulesLvl, vector<BordDeMaille *> *bordsLvl,
  const vector<PhysAdd*> &physAdd, Modele *modele, int &nbMaillesTotalAMR, vector<DomaineGeometrique*> &domaines,	Cellule **cellules, Eos **eos)
{
  nbMaillesTotalAMR = m_nombreCellulesCalcul;
  m_L = min(m_dX, m_dY);
  m_L = min(m_L, m_dZ);

  for (int iterInit = 0; iterInit < 2; iterInit++) {
    for (int lvl = 0; lvl < m_lvlMax; lvl++) {
			if (Ncpu > 1) { Calcul_Parallele.communicationsPrimitivesAMR(cellules, eos, lvl); }
      this->procedureRaffinement(cellulesLvl, bordsLvl, lvl, physAdd, modele, nbMaillesTotalAMR, cellules, eos);
      for (unsigned int i = 0; i < cellulesLvl[lvl + 1].size(); i++) {
        cellulesLvl[lvl + 1][i]->rempli(domaines);
      }
      for (unsigned int i = 0; i < cellulesLvl[lvl + 1].size(); i++) {
        cellulesLvl[lvl + 1][i]->calculsEtendus(m_nombrePhases);
      }
      for (unsigned int i = 0; i < cellulesLvl[lvl].size(); i++) {
        cellulesLvl[lvl][i]->moyenneEnfantsDansParent();
      }
    }
  }
  for (int lvl = 0; lvl <= m_lvlMax; lvl++) {
		if (Ncpu > 1) { Calcul_Parallele.communicationsPrimitivesAMR(cellules, eos, lvl); }
    for (unsigned int i = 0; i < cellulesLvl[lvl].size(); i++) {
      if (!cellulesLvl[lvl][i]->getSplit()) { cellulesLvl[lvl][i]->calculsEtendus(m_nombrePhases); }
    }
  }
}

//***********************************************************************

void MaillageCartesienAMR::procedureRaffinement(vector<Cellule *> *cellulesLvl, vector<BordDeMaille *> *bordsLvl, const int &lvl,
  const vector<PhysAdd*> &physAdd, Modele *modele, int &nbMaillesTotalAMR, Cellule **cellules, Eos **eos)
{
  //1) Calcul de Xi dans chaque cellule de niveau lvl
  //-------------------------------------------------
  for (unsigned int i = 0; i < cellulesLvl[lvl].size(); i++) { cellulesLvl[lvl][i]->miseAZeroXi(); }
  //bool varP2(true);
  //if (lvl >= 3) { varP2 = false; }
  //for (unsigned int i = 0; i < bordsLvl[lvl].size(); i++) { bordsLvl[lvl][i]->calculXi(m_critereVar, m_varRho, varP2, m_varU, m_varAlpha); }
  for (unsigned int i = 0; i < bordsLvl[lvl].size(); i++) { bordsLvl[lvl][i]->calculXi(m_critereVar, m_varRho, m_varP, m_varU, m_varAlpha); }
  if (Ncpu > 1) { Calcul_Parallele.communicationsXi(cellules, lvl); }

  //2) Smoothing de Xi
  //------------------
  for (int iterDiff = 0; iterDiff < 4; iterDiff++) { //Nombre d iterations arbitraire
		//Mise a zero cons xi
    for (unsigned int i = 0; i < cellulesLvl[lvl].size(); i++) { cellulesLvl[lvl][i]->miseAZeroConsXi(); }

    //Calcul des flux
    double dXLocal, dYLocal, dZLocal;
    dXLocal = m_dX / pow(2., (double)lvl);
    dYLocal = m_dY / pow(2., (double)lvl);
    dZLocal = m_dZ / pow(2., (double)lvl);
    for (unsigned int i = 0; i < bordsLvl[lvl].size(); i++) { bordsLvl[lvl][i]->calculFluxXi(dXLocal, dYLocal, dZLocal); }

    //Evolution temporelle
    double dLLocal = m_L / pow(2., (double)lvl);
    double cflDiff = 0.4;
    double dtDiff = cflDiff*0.5*pow(dLLocal, 2.);
    for (unsigned int i = 0; i < cellulesLvl[lvl].size(); i++) { cellulesLvl[lvl][i]->evolutionTemporelleXi(dtDiff); }
		if (Ncpu > 1) { Calcul_Parallele.communicationsXi(cellules, lvl); }
  }

	if (lvl < m_lvlMax) {
    int lvlPlus1 = lvl + 1;
    //3) Raffinement des cellules et bords
    //------------------------------------
    for (unsigned int i = 0; i < cellulesLvl[lvl].size(); i++) { cellulesLvl[lvl][i]->choixRaffine(m_xiSplit, m_nombreMaillesY, m_nombreMaillesZ, m_dX, m_dY, m_dZ, physAdd, modele, nbMaillesTotalAMR); }

    //4) Deraffinement des cellules et bords
    //--------------------------------------
    bool deraffine = false;
    for (unsigned int i = 0; i < cellulesLvl[lvl].size(); i++) { cellulesLvl[lvl][i]->choixDeraffine(m_xiJoin, nbMaillesTotalAMR); }

    if (Ncpu > 1) {
      //5) Raffinement et deraffinement des cellules fantomes
      //-----------------------------------------------------
      //Communication split + Raffinement et deraffinement des cellules fantomes + Reconstruction du tableau de cellules fantomes de niveau lvl + 1
      Calcul_Parallele.communicationsSplit(cellules, lvl);
      m_cellulesLvlGhost[lvlPlus1].clear();
      for (unsigned int i = 0; i < m_cellulesLvlGhost[lvl].size(); i++) { m_cellulesLvlGhost[lvl][i]->choixRaffineDeraffineGhost(m_nombreMaillesY, m_nombreMaillesZ, m_dX, m_dY, m_dZ, physAdd, modele, m_cellulesLvlGhost); }
      //Communications primitives pour mettre a jour les cellules deraffinees
      Calcul_Parallele.communicationsPrimitivesAMR(cellules, eos, lvl);

      //6) Mise a jour des communications persistantes au niveau lvl + 1
      //----------------------------------------------------------------
      Calcul_Parallele.communicationsNombreCellulesGhost(cellules, lvlPlus1);	//Communication des nombres d'elements a envoyer et a recevoir de chaque cote de la limite parallele
      Calcul_Parallele.miseAJourCommunicationsPersistantesLvl(lvlPlus1, m_geometrie);
    }

    //7) Reconstruction des tableaux de cellules et bords lvl + 1
    //-----------------------------------------------------------
    cellulesLvl[lvlPlus1].clear();
    bordsLvl[lvlPlus1].clear();
    for (unsigned int i = 0; i < cellulesLvl[lvl].size(); i++) { cellulesLvl[lvl][i]->constructionTableauxCellulesLvlEtBordsInternesLvl(cellulesLvl, bordsLvl); }
    for (unsigned int i = 0; i < bordsLvl[lvl].size(); i++) { bordsLvl[lvl][i]->constructionTableauBordsExternesLvl(bordsLvl); }
  }
}

//***********************************************************************

string MaillageCartesienAMR::quiSuisJe() const
{
  return "CARTESIEN_AMR";
}

//**************************************************************************
//******************************** ECRITURE ********************************
//**************************************************************************

void MaillageCartesienAMR::ecritEntetePiece(std::ofstream &fluxFichier, std::vector<Cellule *> *cellulesLvl, int lvl) const
{
  int nombreMailles = 0, nombrePointsParMaille = 4;
  for (unsigned int i = 0; i < cellulesLvl[lvl].size(); i++) {
    if (!cellulesLvl[lvl][i]->getSplit()) { nombreMailles += 1; }
  }
  if (m_nombreMaillesZ > 1) { nombrePointsParMaille = 8; }

  fluxFichier << "    <Piece NumberOfPoints=\"" << nombrePointsParMaille*nombreMailles << "\" NumberOfCells=\"" << nombreMailles << "\">" << endl;
}

//***********************************************************************

void MaillageCartesienAMR::recupereNoeuds(std::vector<double> &jeuDonnees, int lvl) const
{
  int dimZ = 0;
  if (m_nombreMaillesZ > 1) dimZ = 1;

  double dXsur2, dYsur2, dZsur2;
  dXsur2 = m_dX / pow(2., (double)lvl) / 2.;
  dYsur2 = m_dY / pow(2., (double)lvl) / 2.;
  dZsur2 = m_dZ / pow(2., (double)lvl) / 2.;

  for (unsigned int i = 0; i < (*m_cellules)[lvl].size(); i++) {
    if (!(*m_cellules)[lvl][i]->getSplit()) {

      //Point 0
      jeuDonnees.push_back((*m_cellules)[lvl][i]->getElement()->getPosition().getX() - dXsur2);
      jeuDonnees.push_back((*m_cellules)[lvl][i]->getElement()->getPosition().getY() - dYsur2);
      jeuDonnees.push_back((*m_cellules)[lvl][i]->getElement()->getPosition().getZ() - dZsur2*dimZ);
      //Point 1
      jeuDonnees.push_back((*m_cellules)[lvl][i]->getElement()->getPosition().getX() + dXsur2);
      jeuDonnees.push_back((*m_cellules)[lvl][i]->getElement()->getPosition().getY() - dYsur2);
      jeuDonnees.push_back((*m_cellules)[lvl][i]->getElement()->getPosition().getZ() - dZsur2*dimZ);
      //Point 2
      jeuDonnees.push_back((*m_cellules)[lvl][i]->getElement()->getPosition().getX() + dXsur2);
      jeuDonnees.push_back((*m_cellules)[lvl][i]->getElement()->getPosition().getY() + dYsur2);
      jeuDonnees.push_back((*m_cellules)[lvl][i]->getElement()->getPosition().getZ() - dZsur2*dimZ);
      //Point 3
      jeuDonnees.push_back((*m_cellules)[lvl][i]->getElement()->getPosition().getX() - dXsur2);
      jeuDonnees.push_back((*m_cellules)[lvl][i]->getElement()->getPosition().getY() + dYsur2);
      jeuDonnees.push_back((*m_cellules)[lvl][i]->getElement()->getPosition().getZ() - dZsur2*dimZ);

      if (dimZ > 0.99) {
        //Point 4
        jeuDonnees.push_back((*m_cellules)[lvl][i]->getElement()->getPosition().getX() - dXsur2);
        jeuDonnees.push_back((*m_cellules)[lvl][i]->getElement()->getPosition().getY() - dYsur2);
        jeuDonnees.push_back((*m_cellules)[lvl][i]->getElement()->getPosition().getZ() + dZsur2);
        //Point 5
        jeuDonnees.push_back((*m_cellules)[lvl][i]->getElement()->getPosition().getX() + dXsur2);
        jeuDonnees.push_back((*m_cellules)[lvl][i]->getElement()->getPosition().getY() - dYsur2);
        jeuDonnees.push_back((*m_cellules)[lvl][i]->getElement()->getPosition().getZ() + dZsur2);
        //Point 6
        jeuDonnees.push_back((*m_cellules)[lvl][i]->getElement()->getPosition().getX() + dXsur2);
        jeuDonnees.push_back((*m_cellules)[lvl][i]->getElement()->getPosition().getY() + dYsur2);
        jeuDonnees.push_back((*m_cellules)[lvl][i]->getElement()->getPosition().getZ() + dZsur2);
        //Point 7
        jeuDonnees.push_back((*m_cellules)[lvl][i]->getElement()->getPosition().getX() - dXsur2);
        jeuDonnees.push_back((*m_cellules)[lvl][i]->getElement()->getPosition().getY() + dYsur2);
        jeuDonnees.push_back((*m_cellules)[lvl][i]->getElement()->getPosition().getZ() + dZsur2);

      }
    } //Fin cellule non split
  } //Fin Cellules
}

//***********************************************************************

void MaillageCartesienAMR::recupereConnectivite(std::vector<double> &jeuDonnees, int lvl) const
{
  int dimZ(0);
  if (m_nombreMaillesZ > 1) dimZ = 1;
  int nombrePointsParMaille(4);
  if (m_nombreMaillesZ > 1) { nombrePointsParMaille = 8; }

  if (dimZ < 0.99) {
    int numCell(0);
    for (unsigned int i = 0; i < (*m_cellules)[lvl].size(); i++) {
      if (!(*m_cellules)[lvl][i]->getSplit()) {
        jeuDonnees.push_back(numCell*nombrePointsParMaille);
        jeuDonnees.push_back(numCell*nombrePointsParMaille+1);
        jeuDonnees.push_back(numCell*nombrePointsParMaille+2);
        jeuDonnees.push_back(numCell*nombrePointsParMaille+3);
        numCell++;
      }
    }
  }
  else {
    int numCell(0);
    for (unsigned int i = 0; i < (*m_cellules)[lvl].size(); i++) {
      if (!(*m_cellules)[lvl][i]->getSplit()) {
        jeuDonnees.push_back(numCell*nombrePointsParMaille);
        jeuDonnees.push_back(numCell*nombrePointsParMaille + 1);
        jeuDonnees.push_back(numCell*nombrePointsParMaille + 2);
        jeuDonnees.push_back(numCell*nombrePointsParMaille + 3);
        jeuDonnees.push_back(numCell*nombrePointsParMaille + 4);
        jeuDonnees.push_back(numCell*nombrePointsParMaille + 5);
        jeuDonnees.push_back(numCell*nombrePointsParMaille + 6);
        jeuDonnees.push_back(numCell*nombrePointsParMaille + 7);
        numCell++;
      }
    }
  }
}

//***********************************************************************

void MaillageCartesienAMR::recupereOffsets(std::vector<double> &jeuDonnees, int lvl) const
{
  int nombrePointsParMaille(4);
  if (m_nombreMaillesZ > 1) { nombrePointsParMaille = 8; }
  int numCell(0);
  for (unsigned int i = 0; i < (*m_cellules)[lvl].size(); i++) {
    if (!(*m_cellules)[lvl][i]->getSplit()) {
      jeuDonnees.push_back((numCell + 1)*nombrePointsParMaille);
      numCell++;
    }
  }
}

//****************************************************************************

void MaillageCartesienAMR::recupereTypeCell(std::vector<double> &jeuDonnees, int lvl) const
{
  int type(9);
  if (m_nombreMaillesZ > 1) { type = 12; }
  int numCell(0);
  for (unsigned int i = 0; i < (*m_cellules)[lvl].size(); i++) {
    if (!(*m_cellules)[lvl][i]->getSplit()) {
      jeuDonnees.push_back(type);
      numCell++;
    }
  }
}

//***********************************************************************

void MaillageCartesienAMR::recupereDonnees(vector<Cellule *> *cellulesLvl, std::vector<double> &jeuDonnees, const int var, int phase, int lvl) const
{
  for (unsigned int i = 0; i < cellulesLvl[lvl].size(); i++) {
    if (!cellulesLvl[lvl][i]->getSplit()) {
      if (var > 0) { //On veut recuperer les donnees scalaires
        if (phase >= 0) { jeuDonnees.push_back(cellulesLvl[lvl][i]->getPhase(phase)->renvoieScalaire(var)); } //Donnees de phases
        else if (phase == -1) { jeuDonnees.push_back(cellulesLvl[lvl][i]->getMelange()->renvoieScalaire(var)); }               //Donnees de melange
        else if (phase == -2) { jeuDonnees.push_back(cellulesLvl[lvl][i]->getTransport(var - 1).getValeur()); }
        else if (phase == -3) { jeuDonnees.push_back(cellulesLvl[lvl][i]->getXi()); }
        else if (phase == -4) { jeuDonnees.push_back(cellulesLvl[lvl][i]->getGradient()); }
        else { Erreurs::messageErreur("MaillageCartesienAMR::recupereDonnees : numero de phase inconnu : ", phase); }
      }
      else { //On veut recuperer les donnees vectorielles
        if (phase >= 0) { //Donnees de phases
          jeuDonnees.push_back(cellulesLvl[lvl][i]->getPhase(phase)->renvoieVecteur(-var)->getX());
          jeuDonnees.push_back(cellulesLvl[lvl][i]->getPhase(phase)->renvoieVecteur(-var)->getY());
          jeuDonnees.push_back(cellulesLvl[lvl][i]->getPhase(phase)->renvoieVecteur(-var)->getZ());
        }
        else if(phase == -1){  //Donnees de melange
          jeuDonnees.push_back(cellulesLvl[lvl][i]->getMelange()->renvoieVecteur(-var)->getX());
          jeuDonnees.push_back(cellulesLvl[lvl][i]->getMelange()->renvoieVecteur(-var)->getY());
          jeuDonnees.push_back(cellulesLvl[lvl][i]->getMelange()->renvoieVecteur(-var)->getZ());
        }
        else { Erreurs::messageErreur("MaillageCartesienAMR::recupereDonnees : numero de phase inconnu : ", phase); }
      } //Fin vecteur
    } //Fin split
  } //fin lvl
}

//****************************************************************************
//****************************** Parallele ***********************************
//****************************************************************************

void MaillageCartesienAMR::initialiseCommunicationsPersistantes(const int nombrePhases, const int nombreTransports, Cellule **cellules, string ordreCalcul)
{
	m_nombrePhases = nombrePhases;
	m_nombreTransports = nombreTransports;
	int nombreVariablesPhaseATransmettre = cellules[0]->getPhase(0)->nombreVariablesATransmettre();
	nombreVariablesPhaseATransmettre *= m_nombrePhases;
	int nombreVariablesMelangeATransmettre = cellules[0]->getMelange()->nombreVariablesATransmettre();
	int m_nombreVariablesPrimitives = nombreVariablesPhaseATransmettre + nombreVariablesMelangeATransmettre + m_nombreTransports;
  int m_nombreVariablesPentes(0);
  if (ordreCalcul == "ORDRE2") {
    int nombrePentesPhaseATransmettre = cellules[0]->getPhase(0)->nombrePentesATransmettre();
    nombrePentesPhaseATransmettre *= m_nombrePhases;
    int nombrePentesMelangeATransmettre = cellules[0]->getMelange()->nombrePentesATransmettre();
    m_nombreVariablesPentes = nombrePentesPhaseATransmettre + nombrePentesMelangeATransmettre + m_nombreTransports;
  }
	Calcul_Parallele.initialiseCommunicationsPersistantesAMR(m_nombreVariablesPrimitives, m_nombreVariablesPentes, m_nombreTransports, m_geometrie, m_lvlMax);
}

//***********************************************************************

void MaillageCartesienAMR::communicationsPrimitives(Cellule **cellules, Eos **eos, const int &lvl, Prim type)
{
	Calcul_Parallele.communicationsPrimitivesAMR(cellules, eos, lvl, type);
}

//***********************************************************************

void MaillageCartesienAMR::communicationsPentes(Cellule **cellules, const int &lvl)
{
	Calcul_Parallele.communicationsPentesAMR(cellules, lvl);
}

//***********************************************************************

void MaillageCartesienAMR::communicationsScalaire(Cellule **cellules, string nomScalaire, const int &lvl)
{
	Calcul_Parallele.communicationsScalaireAMR(cellules, nomScalaire, lvl);
}

//***********************************************************************

void MaillageCartesienAMR::communicationsVecteur(Cellule **cellules, string nomVecteur, const int &dim, const int &lvl, int num, int indice)
{
	Calcul_Parallele.communicationsVecteurAMR(cellules, nomVecteur, m_geometrie, lvl, num, indice);
}

//***********************************************************************

void MaillageCartesienAMR::communicationsPhysAdd(const vector<PhysAdd*> &physAdd, Cellule **cellules, const int &lvl)
{
	for (unsigned int pa = 0; pa < physAdd.size(); pa++) { physAdd[pa]->communicationsPhysAddAMR(cellules, m_geometrie, lvl); }
}

//***********************************************************************

void MaillageCartesienAMR::communicationsTransports(Cellule **cellules, const int &lvl)
{
  Calcul_Parallele.communicationsTransportsAMR(cellules, lvl);
}

//***********************************************************************

void MaillageCartesienAMR::finaliseParallele(const int &lvlMax)
{
	Calcul_Parallele.finaliseAMR(lvlMax);
}

//***********************************************************************
