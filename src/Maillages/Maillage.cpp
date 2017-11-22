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

#include "Maillage.h"

using namespace std;

//***********************************************************************

Maillage::Maillage() :
  m_numFichier(0)
{}

//***********************************************************************

Maillage::~Maillage(){}

//***********************************************************************

int Maillage::getNombreCellules() const
{
  return m_nombreCellulesCalcul;
}

//***********************************************************************

int Maillage::getNombreCellulesTotal() const
{
  return m_nombreCellulesTotal;
}

//***********************************************************************

int Maillage::getNombreFaces() const
{
  return m_nombreFacesTotal;
}

//***********************************************************************

int Maillage::getNumFichier() const
{
  return m_numFichier;
}

//***********************************************************************

void Maillage::ecritSolutionGnuplot(std::vector<Cellule *> *cellulesLvl, std::ofstream &fluxFichier, ObjetGeometrique *objet) const
{
  for (unsigned int c = 0; c < cellulesLvl[0].size(); c++) {
    cellulesLvl[0][c]->ecritureGnuplotAMR(fluxFichier, m_geometrie, objet);
  }
}

//****************************************************************************
//***************************** Methode AMR **********************************
//****************************************************************************

void Maillage::genereTableauxCellulesBordsLvl(Cellule **cellules, BordDeMaille **bord, vector<Cellule *> **cellulesLvl,
  vector<BordDeMaille *> **bordsLvl)
{
  (*cellulesLvl) = new vector<Cellule *>[1];
  for (int i = 0; i < m_nombreCellulesCalcul; i++) { (*cellulesLvl)[0].push_back(cellules[i]); }

  (*bordsLvl) = new vector<BordDeMaille *>[1];
  for (int i = 0; i < m_nombreFacesTotal; i++) { (*bordsLvl)[0].push_back(bord[i]); }

  //m_cellules = cellulesLvl; //KS//FP// Utile !??? Si vraiment utile, que dans MaillageCartesienAMR
}


//****************************************************************************
//****************************** Parallele ***********************************
//****************************************************************************

void Maillage::initialiseCommunicationsPersistantes(const int nombrePhases, const int nombreTransports, Cellule **cellules, string ordreCalcul)
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
	Calcul_Parallele.initialiseCommunicationsPersistantes(m_nombreVariablesPrimitives, m_nombreVariablesPentes, m_nombreTransports, m_geometrie);
}

//***********************************************************************

void Maillage::communicationsPrimitives(Cellule **cellules, Eos **eos, const int &lvl, Prim type)
{
	Calcul_Parallele.communicationsPrimitives(cellules, eos, type);
}

//***********************************************************************

void Maillage::communicationsPentes(Cellule **cellules, const int &lvl)
{
	Calcul_Parallele.communicationsPentes(cellules);
}

//***********************************************************************

void Maillage::communicationsScalaire(Cellule **cellules, string nomScalaire, const int &lvl)
{
	Calcul_Parallele.communicationsScalaire(cellules, nomScalaire);
}

//***********************************************************************

void Maillage::communicationsVecteur(Cellule **cellules, string nomVecteur, const int &dim, const int &lvl, int num, int indice)
{
	Calcul_Parallele.communicationsVecteur(cellules, nomVecteur, m_geometrie, num, indice);
}

//***********************************************************************

void Maillage::communicationsPhysAdd(const vector<PhysAdd*> &physAdd, Cellule **cellules, const int &lvl)
{
	for (unsigned int pa = 0; pa < physAdd.size(); pa++) { physAdd[pa]->communicationsPhysAdd(cellules, m_geometrie); }
}

//***********************************************************************

void Maillage::communicationsTransports(Cellule **cellules, const int &lvl)
{
  Calcul_Parallele.communicationsTransports(cellules);
}


//***********************************************************************

void Maillage::finaliseParallele(const int &lvlMax)
{
	Calcul_Parallele.finalise(lvlMax);
}

//***********************************************************************