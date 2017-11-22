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

#include "CelluleO2Ghost.h"
#include <iostream>

using namespace std;

//***********************************************************************

CelluleO2Ghost::CelluleO2Ghost() : CelluleO2(), m_vecPhasesPentesGhost(0), m_melangePentesGhost(0), m_vecTransportsPentesGhost(0) {}

//***********************************************************************

CelluleO2Ghost::CelluleO2Ghost(int lvl) : CelluleO2(lvl), m_vecPhasesPentesGhost(0), m_melangePentesGhost(0), m_vecTransportsPentesGhost(0) {}

//***********************************************************************

CelluleO2Ghost::~CelluleO2Ghost()
{
	for (int k = 0; k < m_nombrePhases; k++) {
		delete m_vecPhasesPentesGhost[k];
	}
	delete[] m_vecPhasesPentesGhost;
	delete m_melangePentesGhost;
	delete[] m_vecTransportsPentesGhost;
}

//***********************************************************************

void CelluleO2Ghost::alloue(const int &nombrePhases, const int &nombreTransports, const std::vector<PhysAdd*> &physAdd, Modele *modele)
{
	m_nombrePhases = nombrePhases;
	m_nombreTransports = nombreTransports;
	m_vecPhases = new Phase*[nombrePhases];
	m_vecPhasesO2 = new Phase*[nombrePhases];
	for (int k = 0; k < nombrePhases; k++) {
		modele->allouePhase(&m_vecPhases[k]);
		modele->allouePhase(&m_vecPhasesO2[k]);
	}
	modele->alloueMelange(&m_melange);
	modele->alloueMelange(&m_melangeO2);
	modele->alloueCons(&m_cons, nombrePhases);
	modele->alloueCons(&m_consSauvegarde, nombrePhases);
	if (nombreTransports > 0) {
		m_vecTransports = new Transport[nombreTransports];
		m_consTransports = new Transport[nombreTransports];
		m_consTransportsSauvegarde = new Transport[nombreTransports];
		m_vecTransportsO2 = new Transport[nombreTransports];
	}
	for (unsigned int k = 0; k < physAdd.size(); k++) {
		physAdd[k]->ajouteGrandeurPhysAdd(this);
	}

	//Allocation des pentes fantomes, specifique aux limites paralleles
	m_vecPhasesPentesGhost = new Phase*[nombrePhases];
	for (int k = 0; k < nombrePhases; k++) {
		modele->allouePhase(&m_vecPhasesPentesGhost[k]);
		m_vecPhasesPentesGhost[k]->miseAZero();
	}
	modele->alloueMelange(&m_melangePentesGhost);
	m_melangePentesGhost->miseAZero();
	m_vecTransportsPentesGhost = new double[nombreTransports];
	for (int k = 0; k < nombreTransports; k++) {
		m_vecTransportsPentesGhost[k] = 0.;
	}
}

//***********************************************************************

void CelluleO2Ghost::calculPentesLocal(const int &nombrePhases, const int &nombreTransports, BordDeMaille &bordRef, Limiteur &limiteurGlobal, Limiteur &limiteurInterface)
{
	//Mise a zero des pentes locales
	//------------------------------
	double sommeCoeff(0.);
	for (int k = 0; k < nombrePhases; k++) {
		pentesPhasesLocal1[k]->miseAZero();
	}
	pentesMelangeLocal1->miseAZero();
	for (int k = 0; k < nombreTransports; k++) {
		pentesTransportLocal1[k] = 0.;
	}

	//Boucle sur les bords pour la détermination de la pente du cote de bordRef
	//-------------------------------------------------------------------------
	for (unsigned int b = 0; b < m_bords.size(); b++) {
		if (!m_bords[b]->getSplit()) {
			for (int k = 0; k < nombrePhases; k++) { pentesPhasesLocal1[k]->multiplieEtAjoute(*m_bords[b]->getPentesPhase(k), 1.); }
			pentesMelangeLocal1->multiplieEtAjoute(*m_bords[b]->getPentesMelange(), 1.);
			for (int k = 0; k < nombreTransports; k++) { pentesTransportLocal1[k] += m_bords[b]->getPentesTransport(k)->getValeur(); }
			sommeCoeff += 1.;
		}
	}

	//Normalisation des pentes
	//------------------------
	if (sommeCoeff > 1.e-5) {
		for (int k = 0; k < nombrePhases; k++) { pentesPhasesLocal1[k]->diviser(sommeCoeff); }
		pentesMelangeLocal1->diviser(sommeCoeff);
		for (int k = 0; k < nombreTransports; k++) { pentesTransportLocal1[k] /= sommeCoeff; }
	}

	//Limitations des pentes
	//----------------------
	for (int k = 0; k < nombrePhases; k++) {
		pentesPhasesLocal1[k]->limitePentes(*pentesPhasesLocal1[k], *m_vecPhasesPentesGhost[k], limiteurGlobal, limiteurInterface);
	}
	pentesMelangeLocal1->limitePentes(*pentesMelangeLocal1, *m_melangePentesGhost, limiteurGlobal);
	for (int k = 0; k < nombreTransports; k++) {
		pentesTransportLocal1[k] = limiteurInterface.limitePente(pentesTransportLocal1[k], m_vecTransportsPentesGhost[k]);
	}
}

//***********************************************************************

void CelluleO2Ghost::creerCelluleEnfant(const int &num, const int &lvl)
{
	m_cellulesEnfants.push_back(new CelluleO2Ghost(lvl + 1));
}

//***********************************************************************

void CelluleO2Ghost::recupereTamponPentes(double *tampon, int &compteur)
{
	for (int k = 0; k < m_nombrePhases; k++) {
		m_vecPhasesPentesGhost[k]->recupereTamponPentes(tampon, compteur);
	}
	m_melangePentesGhost->recupereTamponPentes(tampon, compteur);
	for (int k = 0; k < m_nombreTransports; k++) {
		m_vecTransportsPentesGhost[k] = tampon[++compteur];
	}
}

//***********************************************************************

void CelluleO2Ghost::recupereTamponPentesAMR(double *tampon, int &compteur, const int &lvl)
{
	if (m_lvl == lvl) {
		for (int k = 0; k < m_nombrePhases; k++) {
			m_vecPhasesPentesGhost[k]->recupereTamponPentes(tampon, compteur);
		}
		m_melangePentesGhost->recupereTamponPentes(tampon, compteur);
		for (int k = 0; k < m_nombreTransports; k++) {
			m_vecTransportsPentesGhost[k] = tampon[++compteur];
		}
	}
	else {
		for (unsigned int i = 0; i < m_cellulesEnfants.size(); i++) {
			m_cellulesEnfants[i]->recupereTamponPentesAMR(tampon, compteur, lvl);
		}
	}
}

//***********************************************************************