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

#include "CelluleO2.h"
#include <iostream>
#include "../Modeles/Phase.h"

using namespace std;

//***********************************************************************

CelluleO2::CelluleO2() : Cellule(), m_vecPhasesO2(0), m_melangeO2(0), m_vecTransportsO2(0), m_consSauvegarde(0), m_consTransportsSauvegarde(0) {}

//***********************************************************************

CelluleO2::CelluleO2(int lvl) : Cellule(lvl), m_vecPhasesO2(0), m_melangeO2(0), m_vecTransportsO2(0), m_consSauvegarde(0), m_consTransportsSauvegarde(0) {}

//***********************************************************************

CelluleO2::~CelluleO2()
{
  for (int k = 0; k < m_nombrePhases; k++) {
    delete m_vecPhasesO2[k];
  }
  delete[] m_vecPhasesO2;
  if (m_vecTransportsO2 != 0) delete[] m_vecTransportsO2;
  delete m_melangeO2;
  delete m_consSauvegarde;
	if (m_consTransportsSauvegarde != 0) delete[] m_consTransportsSauvegarde;
}

//***********************************************************************

void CelluleO2::alloue(const int &nombrePhases, const int &nombreTransports, const std::vector<PhysAdd*> &physAdd, Modele *modele)
{
  m_nombrePhases = nombrePhases;
  m_nombreTransports = nombreTransports;
  m_vecPhases = new Phase*[nombrePhases];
  m_vecPhasesO2 = new Phase*[nombrePhases];
  for (int k = 0; k < nombrePhases; k++){
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
}

//***********************************************************************

void CelluleO2::copiePhase(const int &numeroPhase, Phase *phase)
{
  m_vecPhases[numeroPhase]->copiePhase(*phase);
  m_vecPhasesO2[numeroPhase]->copiePhase(*phase);
}

//***********************************************************************

void CelluleO2::calculPentesLocal(const int &nombrePhases, const int &nombreTransports, BordDeMaille &bordRef, Limiteur &limiteurGlobal, Limiteur &limiteurInterface)
{
	//Solution pour multiD cartesien (peut etre une ebauche pour le NS, a voir...)

	//Mise a zero des pentes locales
	//------------------------------
  double coeff(0.), posBordRef(0.);
	double sommeCoeff(0.), sommeCoeff2(0.);
	for (int k = 0; k < nombrePhases; k++) {
		pentesPhasesLocal1[k]->miseAZero();
		pentesPhasesLocal2[k]->miseAZero();
	}
	pentesMelangeLocal1->miseAZero();
	pentesMelangeLocal2->miseAZero();
	for (int k = 0; k < nombreTransports; k++) {
		pentesTransportLocal1[k] = 0.;
		pentesTransportLocal2[k] = 0.;
	}
	//Boucle sur les bords pour la détermination des pentes de chaque cote de la cellule
	//----------------------------------------------------------------------------------
	for (unsigned int b = 0; b < m_bords.size(); b++) {
    //Calcul de la pente a gauche et a droite de la cellule (AMR) en se basant sur celle de reference (bord a gauche ou a droite, inconnu)
    if (m_bords[b]->getPentesPhase(0) != 0) { //Bord de type BordDeMailleO2 ou CondLimMurO2
      if (m_bords[b] == &bordRef) {
				for (int k = 0; k < nombrePhases; k++) { pentesPhasesLocal1[k]->multiplieEtAjoute(*m_bords[b]->getPentesPhase(k), 1.); }
				pentesMelangeLocal1->multiplieEtAjoute(*m_bords[b]->getPentesMelange(), 1.);
				for (int k = 0; k < nombreTransports; k++) { pentesTransportLocal1[k] += m_bords[b]->getPentesTransport(k)->getValeur(); }
        sommeCoeff += 1.;
      }
      else {
				if (!m_bords[b]->getSplit()) {
				//Produit scalaire des normales avec celle de reference
					coeff = fabs(m_bords[b]->getFace()->getNormale().scalaire(bordRef.getFace()->getNormale()));
					if (coeff > 1.e-6) {
						//Face majoritement selon X
						if (fabs(bordRef.getFace()->getNormale().getX()) > 0.5) {
							posBordRef = bordRef.getFace()->getPos().getX();
							//Cote bordRef
							if (fabs(posBordRef - m_bords[b]->getFace()->getPos().getX()) <= fabs(posBordRef - m_element->getPosition().getX())) {
								for (int k = 0; k < nombrePhases; k++) { pentesPhasesLocal1[k]->multiplieEtAjoute(*m_bords[b]->getPentesPhase(k), coeff); }
								pentesMelangeLocal1->multiplieEtAjoute(*m_bords[b]->getPentesMelange(), coeff);
								for (int k = 0; k < nombreTransports; k++) { pentesTransportLocal1[k] += coeff*(m_bords[b]->getPentesTransport(k)->getValeur()); }
								sommeCoeff += coeff;
							}
							//Autre cote
							else {
								for (int k = 0; k < nombrePhases; k++) { pentesPhasesLocal2[k]->multiplieEtAjoute(*m_bords[b]->getPentesPhase(k), coeff); }
								pentesMelangeLocal2->multiplieEtAjoute(*m_bords[b]->getPentesMelange(), coeff);
								for (int k = 0; k < nombreTransports; k++) { pentesTransportLocal2[k] += coeff*(m_bords[b]->getPentesTransport(k)->getValeur()); }
								sommeCoeff2 += coeff;
							}
						}
						//Face majoritement selon Y
						else if (fabs(bordRef.getFace()->getNormale().getY()) > 0.5) {
							posBordRef = bordRef.getFace()->getPos().getY();
							//Cote bordRef
							if (fabs(posBordRef - m_bords[b]->getFace()->getPos().getY()) <= fabs(posBordRef - m_element->getPosition().getY())) {
								for (int k = 0; k < nombrePhases; k++) { pentesPhasesLocal1[k]->multiplieEtAjoute(*m_bords[b]->getPentesPhase(k), coeff); }
								pentesMelangeLocal1->multiplieEtAjoute(*m_bords[b]->getPentesMelange(), coeff);
								for (int k = 0; k < nombreTransports; k++) { pentesTransportLocal1[k] += coeff*(m_bords[b]->getPentesTransport(k)->getValeur()); }
								sommeCoeff += coeff;
							}
							//Autre cote
							else {
								for (int k = 0; k < nombrePhases; k++) { pentesPhasesLocal2[k]->multiplieEtAjoute(*m_bords[b]->getPentesPhase(k), coeff); }
								pentesMelangeLocal2->multiplieEtAjoute(*m_bords[b]->getPentesMelange(), coeff);
								for (int k = 0; k < nombreTransports; k++) { pentesTransportLocal2[k] += coeff*(m_bords[b]->getPentesTransport(k)->getValeur()); }
								sommeCoeff2 += coeff;
							}
						}
						//Face majoritement selon Z
						else {
							posBordRef = bordRef.getFace()->getPos().getZ();
							//Cote bordRef
							if (fabs(posBordRef - m_bords[b]->getFace()->getPos().getZ()) <= fabs(posBordRef - m_element->getPosition().getZ())) {
								for (int k = 0; k < nombrePhases; k++) { pentesPhasesLocal1[k]->multiplieEtAjoute(*m_bords[b]->getPentesPhase(k), coeff); }
								pentesMelangeLocal1->multiplieEtAjoute(*m_bords[b]->getPentesMelange(), coeff);
								for (int k = 0; k < nombreTransports; k++) { pentesTransportLocal1[k] += coeff*(m_bords[b]->getPentesTransport(k)->getValeur()); }
								sommeCoeff += coeff;
							}
							//Autre cote
							else {
								for (int k = 0; k < nombrePhases; k++) { pentesPhasesLocal2[k]->multiplieEtAjoute(*m_bords[b]->getPentesPhase(k), coeff); }
								pentesMelangeLocal2->multiplieEtAjoute(*m_bords[b]->getPentesMelange(), coeff);
								for (int k = 0; k < nombreTransports; k++) { pentesTransportLocal2[k] += coeff*(m_bords[b]->getPentesTransport(k)->getValeur()); }
								sommeCoeff2 += coeff;
							}
						}
					}
				}
      }
    }
	} //fin boucle sur les bords

	//Normalisation des pentes
	//------------------------
	if (sommeCoeff > 1.e-8) {
		for (int k = 0; k < nombrePhases; k++) { pentesPhasesLocal1[k]->diviser(sommeCoeff);	}
		pentesMelangeLocal1->diviser(sommeCoeff);
		for (int k = 0; k < nombreTransports; k++) { pentesTransportLocal1[k] /= sommeCoeff; }
	}
	if (sommeCoeff2 > 1.e-8) {
		for (int k = 0; k < nombrePhases; k++) { pentesPhasesLocal2[k]->diviser(sommeCoeff2);	}
		pentesMelangeLocal2->diviser(sommeCoeff2);
		for (int k = 0; k < nombreTransports; k++) { pentesTransportLocal2[k] /= sommeCoeff2; }
	}

	//Limitations des pentes
	//----------------------
  for (int k = 0; k < nombrePhases; k++) {
    pentesPhasesLocal1[k]->limitePentes(*pentesPhasesLocal1[k], *pentesPhasesLocal2[k], limiteurGlobal, limiteurInterface);
  }
  pentesMelangeLocal1->limitePentes(*pentesMelangeLocal1, *pentesMelangeLocal2, limiteurGlobal);
  for (int k = 0; k < nombreTransports; k++) {
    pentesTransportLocal1[k] = limiteurInterface.limitePente(pentesTransportLocal1[k], pentesTransportLocal2[k]);
  }
}

//***********************************************************************

void CelluleO2::calculPentesLocalLimite(const int &nombrePhases, const int &nombreTransports, BordDeMaille &bordRef, Limiteur &limiteurGlobal, Limiteur &limiteurInterface)
{
  //Solution pour multiD cartesien (peut etre une ebauche pour le NS, a voir...)

  //Mise a zero des pentes locales
  //------------------------------
  double coeff(0.), posBordRef(0.);
  double sommeCoeff2(0.);
  for (int k = 0; k < nombrePhases; k++) {
    pentesPhasesLocal1[k]->miseAZero();
    pentesPhasesLocal2[k]->miseAZero();
  }
  pentesMelangeLocal1->miseAZero();
  pentesMelangeLocal2->miseAZero();
  for (int k = 0; k < nombreTransports; k++) {
    pentesTransportLocal1[k] = 0.;
    pentesTransportLocal2[k] = 0.;
  }

  //Recupere la pente cote CL
  //-------------------------
  for (int k = 0; k < nombrePhases; k++) { pentesPhasesLocal1[k]->multiplieEtAjoute(*bordRef.getPentesPhase(k), 1.); }
  pentesMelangeLocal1->multiplieEtAjoute(*bordRef.getPentesMelange(), 1.);
  for (int k = 0; k < nombreTransports; k++) { pentesTransportLocal1[k] += bordRef.getPentesTransport(k)->getValeur(); }

  //Boucle sur les bords pour la détermination de la pente cote oppose a la CL
  //--------------------------------------------------------------------------
  for (unsigned int b = 0; b < m_bords.size(); b++) {
    //Calcul de la pente a gauche et a droite de la cellule (AMR) en se basant sur celle de reference (bord a gauche ou a droite, inconnu)
    if (m_bords[b]->getPentesPhase(0) != 0) { //Bord de type BordDeMaille/O2
      if (m_bords[b] != &bordRef) {
        if (!m_bords[b]->getSplit()) {
          //Produit scalaire des normales avec celle de reference
          coeff = fabs(m_bords[b]->getFace()->getNormale().scalaire(bordRef.getFace()->getNormale()));
          if (coeff > 1.e-6) {
            //Face majoritement selon X
            if (fabs(bordRef.getFace()->getNormale().getX()) > 0.5) {
              posBordRef = bordRef.getFace()->getPos().getX();
              //Autre cote
              if (fabs(posBordRef - m_bords[b]->getFace()->getPos().getX()) >= fabs(posBordRef - m_element->getPosition().getX())) {
                for (int k = 0; k < nombrePhases; k++) { pentesPhasesLocal2[k]->multiplieEtAjoute(*m_bords[b]->getPentesPhase(k), coeff); }
                pentesMelangeLocal2->multiplieEtAjoute(*m_bords[b]->getPentesMelange(), coeff);
                for (int k = 0; k < nombreTransports; k++) { pentesTransportLocal2[k] += coeff*(m_bords[b]->getPentesTransport(k)->getValeur()); }
                sommeCoeff2 += coeff;
              }
            }
            //Face majoritement selon Y
            else if (fabs(bordRef.getFace()->getNormale().getY()) > 0.5) {
              posBordRef = bordRef.getFace()->getPos().getY();
              //Autre cote
              if (fabs(posBordRef - m_bords[b]->getFace()->getPos().getY()) >= fabs(posBordRef - m_element->getPosition().getY())) {
                for (int k = 0; k < nombrePhases; k++) { pentesPhasesLocal2[k]->multiplieEtAjoute(*m_bords[b]->getPentesPhase(k), coeff); }
                pentesMelangeLocal2->multiplieEtAjoute(*m_bords[b]->getPentesMelange(), coeff);
                for (int k = 0; k < nombreTransports; k++) { pentesTransportLocal2[k] += coeff*(m_bords[b]->getPentesTransport(k)->getValeur()); }
                sommeCoeff2 += coeff;
              }
            }
            //Face majoritement selon Z
            else {
              posBordRef = bordRef.getFace()->getPos().getZ();
              //Autre cote
              if (fabs(posBordRef - m_bords[b]->getFace()->getPos().getZ()) >= fabs(posBordRef - m_element->getPosition().getZ())) {
                for (int k = 0; k < nombrePhases; k++) { pentesPhasesLocal2[k]->multiplieEtAjoute(*m_bords[b]->getPentesPhase(k), coeff); }
                pentesMelangeLocal2->multiplieEtAjoute(*m_bords[b]->getPentesMelange(), coeff);
                for (int k = 0; k < nombreTransports; k++) { pentesTransportLocal2[k] += coeff*(m_bords[b]->getPentesTransport(k)->getValeur()); }
                sommeCoeff2 += coeff;
              }
            }
          }
        }
      }
    }
  } //fin boucle sur les bords

  //Normalisation de la pente
  //-------------------------
  if (sommeCoeff2 > 1.e-8) {
    for (int k = 0; k < nombrePhases; k++) { pentesPhasesLocal2[k]->diviser(sommeCoeff2); }
    pentesMelangeLocal2->diviser(sommeCoeff2);
    for (int k = 0; k < nombreTransports; k++) { pentesTransportLocal2[k] /= sommeCoeff2; }
  }

  //Limitations des pentes
  //----------------------
  for (int k = 0; k < nombrePhases; k++) {
    pentesPhasesLocal1[k]->limitePentes(*pentesPhasesLocal1[k], *pentesPhasesLocal2[k], limiteurGlobal, limiteurInterface);
  }
  pentesMelangeLocal1->limitePentes(*pentesMelangeLocal1, *pentesMelangeLocal2, limiteurGlobal);
  for (int k = 0; k < nombreTransports; k++) {
    pentesTransportLocal1[k] = limiteurInterface.limitePente(pentesTransportLocal1[k], pentesTransportLocal2[k]);
  }
}

//***********************************************************************

//void CelluleO2::calculMultiPente(const int &nombrePhases, BordDeMaille *bord, Limiteur *limiteurGlobal)
//{
//  //Calcul des pentes (Le Touze et al. 2014)
//  //pentesPlus correspond aux pentes aval
//  //pentesMoins correspond aux pentes amont
//  Phase *pentesPlus; Phase *pentesMoins;
//  Phase *phaseHM=0, *phaseHP=0;
//  double distanceP, distanceM;
//
//  for (int k = 0; k < nombrePhases; k++) {
//    this->getPhase(k)->initialisePhase(&pentesPlus);
//    this->getPhase(k)->initialisePhase(&pentesMoins);
//    //Si cellule de gauche
//    if (bord->getCellGauche() == this) {
//      if (bord->getB(BG1M) != 0) {
//        this->getPhase(k)->initialisePhase(&phaseHM);
//        double beta1 = bord->getBeta(betaG1M);
//        double beta2 = bord->getBeta(betaG2M);
//        *phaseHM = beta1 * *(bord->getB(BG1M)->getPhase(k)) + beta2 * *(bord->getB(BG2M)->getPhase(k));
//        distanceM = bord->getDistanceH(distanceHGM);
//      }
//      if (bord->getB(BG1P) != 0) {
//        this->getPhase(k)->initialisePhase(&phaseHP);
//        double beta1 = bord->getBeta(betaG1P);
//        double beta2 = bord->getBeta(betaG2P);
//        *phaseHP = beta1 * *(bord->getB(BG1P)->getPhase(k)) + beta2 * *(bord->getB(BG2P)->getPhase(k));
//        distanceP = bord->getDistanceH(distanceHGP);
//      }
//    }
//    //Si cellule de droite
//    else {
//      if (bord->getB(BD1M) != 0) {
//        this->getPhase(k)->initialisePhase(&phaseHM);
//        double beta1 = bord->getBeta(betaD1M);
//        double beta2 = bord->getBeta(betaD2M);
//        *phaseHM = beta1 * *(bord->getB(BD1M)->getPhase(k)) + beta2 * *(bord->getB(BD2M)->getPhase(k));
//        distanceM = bord->getDistanceH(distanceHDM);
//      }
//      if (bord->getB(BD1P) != 0) {
//        this->getPhase(k)->initialisePhase(&phaseHP);
//        double beta1 = bord->getBeta(betaD1P);
//        double beta2 = bord->getBeta(betaD2P);
//        *phaseHP = beta1 * *(bord->getB(BD1P)->getPhase(k)) + beta2 * *(bord->getB(BD2P)->getPhase(k));
//        distanceP = bord->getDistanceH(distanceHDP);
//      }
//    }
//    if (phaseHP != 0) { pentesPlus->calculPentes((this->getPhase(k)), phaseHP, distanceP); }
//    else { pentesPlus->miseAZero(); }
//    if (phaseHM != 0) { pentesMoins->calculPentes(phaseHM, (this->getPhase(k)), distanceM); }
//    else { pentesMoins->miseAZero(); }
//
//    //pentesPlus->miseAZero(); //FP//DEV// test a enlever
//    //pentesMoins->miseAZero();
//
//    //Une fois les pentes + et - calcules, on limite et on stock dans m_vecPhasesPentes
//    //m_vecPhasesPentes[k]->limitePentes(pentesPlus, pentesMoins, limiteurGlobal);
//    delete pentesMoins; delete pentesPlus;
//    if (phaseHM != 0) delete phaseHM;
//    if (phaseHP != 0) delete phaseHP;
//
//  }
//}

//***********************************************************************

void CelluleO2::sauvegardeCons(const int &nombrePhases, const int &nombreTransports)
{
  m_consSauvegarde->setCons(m_cons, nombrePhases);
  for (int k = 0; k < nombreTransports; k++) { m_consTransportsSauvegarde[k].setValeur(m_consTransports[k].getValeur()); }
}

//***********************************************************************

void CelluleO2::recuperationCons(const int &nombrePhases, const int &nombreTransports)
{
  m_cons->setCons(m_consSauvegarde, nombrePhases);
  for (int k = 0; k < nombreTransports; k++) { m_consTransports[k].setValeur(m_consTransportsSauvegarde[k].getValeur()); }
}

//***********************************************************************

void CelluleO2::predictionOrdre2(const double &dt, const int &nombrePhases, const int &nombreTransports)
{
  m_cons->miseEnTampon(*this, nombrePhases);                        //On determine Un grace au vecteur primitif des phases et de melange que l on stocke dans fluxTempXXX
  m_cons->multiplie(0.5*dt, nombrePhases);                          //On multiplie m_cons (bilan des flux) par dt/2
  m_cons->ajoutFlux(1., nombrePhases);                              //On y ajoute fluxTempXXX (Un) -> on obtient Un+1/2 dans m_cons
  m_cons->construitPrim(m_vecPhasesO2, m_melangeO2, nombrePhases);  //On peut reconstruire m_vecPhasesO2 et m_melangeO2 a partir de m_cons
  for (int k = 0; k < nombreTransports; k++) {
    m_consTransports[k].multiplie(0.5*dt);
    m_vecTransportsO2[k].setValeur(m_vecTransports[k].getValeur());
    m_vecTransportsO2[k].ajoute(m_consTransports[k].getValeur());
  }

  // Relaxations et correction des energies //FP//Q// Si je le vire, ca doit marcher quand meme
  //FP//TODO// Faire appel aux modeles (il n'y a pas toujours de la relaxation etc.)
  m_cons->relaxPressions(this, nombrePhases, vecPhasesO2);
  m_melangeO2->energieTotaleVersEnergieInterne(m_vecGrandeursPhysAdd); //On reconstruit l'energie interne a partir de l energie totale
  m_cons->correctionEnergie(this, nombrePhases, vecPhasesO2);
}

//***********************************************************************

void CelluleO2::alloueEtCopiePhase(const int &numeroPhase, Phase *phase)
{
  phase->alloueEtCopiePhase(&m_vecPhases[numeroPhase]);
  phase->alloueEtCopiePhase(&m_vecPhasesO2[numeroPhase]);
}

//***********************************************************************

void CelluleO2::calculsEtendus(const int &nombrePhases, Prim type)
{
  switch (type) {
  case vecPhases: //Idem cellule ordre 1
    for (int k = 0; k < nombrePhases; k++) {
      m_vecPhases[k]->calculsEtendusPhases(m_melange->getVitesse());
    }
    this->preparePhysAdd();
    m_melange->calculGrandeursMelange(m_vecPhases, nombrePhases);
    m_melange->energieInterneVersEnergieTotale(m_vecGrandeursPhysAdd); //Mis a part car depend de la physique.
    break;
  case vecPhasesO2: //Utile seulement pour le parallele ordre 2
    for (int k = 0; k < nombrePhases; k++) {
      m_vecPhasesO2[k]->calculsEtendusPhases(m_melangeO2->getVitesse());
    }
    //Les gradients ne sont pas recalculer ici sous peine de perdre l'energie capillaire.
    m_melangeO2->calculGrandeursMelange(m_vecPhasesO2, nombrePhases);
    m_melangeO2->energieInterneVersEnergieTotale(m_vecGrandeursPhysAdd); //Mis a part car depend de la physique.
    break;
  default: break;
  }
}

//***********************************************************************

void CelluleO2::calculsEtendusPourCommunications(const int &nombrePhases, Prim type)
{
  switch (type) {
  case vecPhases: //Idem cellule ordre 1
    for (int k = 0; k < nombrePhases; k++) {
      m_vecPhases[k]->calculsEtendusPhases(m_melange->getVitesse());
    }
    m_melange->calculGrandeursMelange(m_vecPhases, nombrePhases);
    break;
  case vecPhasesO2: //Utile seulement pour le parallele ordre 2
    for (int k = 0; k < nombrePhases; k++) {
      m_vecPhasesO2[k]->calculsEtendusPhases(m_melangeO2->getVitesse());
    }
    //Les gradients ne sont pas recalculer ici sous peine de perdre l'energie capillaire.
    m_melangeO2->calculGrandeursMelange(m_vecPhasesO2, nombrePhases);
    break;
  default: break;
  }
}

//***********************************************************************

void CelluleO2::projection(const Coord &normale, const Coord &tangente, const Coord &binormale, const int &nombrePhases, Prim type)
{
  switch (type) {
    case vecPhases:
      for (int k = 0; k < nombrePhases; k++) {
        m_vecPhases[k]->projection(normale, tangente, binormale);
      }
      m_melange->projection(normale, tangente, binormale);
      break;
    case vecPhasesO2:
      for (int k = 0; k < nombrePhases; k++) {
        m_vecPhasesO2[k]->projection(normale, tangente, binormale);
      }
      m_melangeO2->projection(normale, tangente, binormale);
      break;
    default: break;
  }
}

//***********************************************************************

void CelluleO2::projectionRepereAbsolu(const Coord &normale, const Coord &tangente, const Coord &binormale, const int &nombrePhases, Prim type)
{
  switch (type) {
    case vecPhases:
      for (int k = 0; k < nombrePhases; k++) {
        m_vecPhases[k]->projectionRepereAbsolu(normale, tangente, binormale);
      }
      m_melange->projectionRepereAbsolu(normale, tangente, binormale);
      break;
    case vecPhasesO2:
      for (int k = 0; k < nombrePhases; k++) {
        m_vecPhasesO2[k]->projectionRepereAbsolu(normale, tangente, binormale);
      }
      m_melangeO2->projectionRepereAbsolu(normale, tangente, binormale);
      break;
    default: break;
  }
}

//***********************************************************************

void CelluleO2::copieDansCellule(Cellule &cellSource, Prim type) const
{
  cellSource = static_cast<Cellule>(*this);
}

//***********************************************************************

Phase* CelluleO2::getPhase(const int &numeroPhase, Prim type) const
{
  switch (type){
    case vecPhases: return m_vecPhases[numeroPhase]; break;
    case vecPhasesO2: return m_vecPhasesO2[numeroPhase]; break;
    default: return 0; break;
  }
}

//***********************************************************************

Phase** CelluleO2::getPhases(Prim type) const
{
  switch (type) {
    case vecPhases: return m_vecPhases; break;
    case vecPhasesO2: return m_vecPhasesO2; break;
    default: return 0; break;
  }
}

//***********************************************************************

Melange* CelluleO2::getMelange(Prim type) const
{
  switch (type) {
    case vecPhases: return m_melange; break;
    case vecPhasesO2: return m_melangeO2; break;
  default: return 0; break;
  }
}

//***********************************************************************

Transport& CelluleO2::getTransport(const int &numTransport, Prim type) const
{
	switch (type) {
	case vecPhases: return m_vecTransports[numTransport]; break;
	case vecPhasesO2: return m_vecTransportsO2[numTransport]; break;
  default: return m_vecTransports[numTransport]; break; //FP//TODO// trouver un moyen plus intelligent de faire les renvoi par defaut sur objets.
	}
}

//***********************************************************************

Transport* CelluleO2::getTransports(Prim type) const
{
  switch (type) {
    case vecPhases: return m_vecTransports; break;
    case vecPhasesO2: return m_vecTransportsO2; break;
    default: return 0; break;
  }
}

//***********************************************************************

void CelluleO2::setTransport(double valeur, int &numTransport, Prim type)
{
  switch (type) {
  case vecPhases: m_vecTransports[numTransport].setValeur(valeur); break;
  case vecPhasesO2: m_vecTransportsO2[numTransport].setValeur(valeur); break;
  default: break;
  }
}

//****************************************************************************
//***************************** Methode AMR **********************************
//****************************************************************************

void CelluleO2::creerCelluleEnfant(const int &num, const int &lvl)
{
  m_cellulesEnfants.push_back(new CelluleO2(lvl + 1));
}

//****************************************************************************
//********************** Methode Ordre 2 Parallele ***************************
//****************************************************************************

void CelluleO2::rempliTamponPentes(double *tampon, int &compteur) const
{
	//Mise a zero de la pente locale a transmettre
	//--------------------------------------------
	int bordRef(0);
	double coeff(0.), sommeCoeff(0.);
	for (int k = 0; k < m_nombrePhases; k++) {
		pentesPhasesLocal1[k]->miseAZero();
	}
	pentesMelangeLocal1->miseAZero();
	for (int k = 0; k < m_nombreTransports; k++) {
		pentesTransportLocal1[k] = 0.;
	}

	//Boucle sur les bords pour la détermination du bord de reference (liaison avec CelluleO2Ghost)
	//---------------------------------------------------------------------------------------------
	for (unsigned int b = 0; b < m_bords.size(); b++) {
		if (m_bords[b]->getCellGauche()->estCelluleO2Ghost() || m_bords[b]->getCellDroite()->estCelluleO2Ghost()) {
			bordRef = b;
			break;
		}
	}

	//Boucle sur les bords pour la détermination de la pente a transmettre
	//--------------------------------------------------------------------
	for (unsigned int b = 0; b < m_bords.size(); b++) {
		if (b == !bordRef) {
			//Produit scalaire des normales avec celle de reference
			coeff = m_bords[b]->getFace()->getNormale().scalaire(m_bords[bordRef]->getFace()->getNormale());
			if (coeff > 1.e-6) {
				for (int k = 0; k < m_nombrePhases; k++) { pentesPhasesLocal1[k]->multiplieEtAjoute(*m_bords[b]->getPentesPhase(k), coeff); }
				pentesMelangeLocal1->multiplieEtAjoute(*m_bords[b]->getPentesMelange(), coeff);
				for (int k = 0; k < m_nombreTransports; k++) { pentesTransportLocal1[k] += coeff*(m_bords[b]->getPentesTransport(k)->getValeur()); }
				sommeCoeff += coeff;
			}
		}
	}
	
	//Normalisation de la pente
	//-------------------------
	if (sommeCoeff > 1.e-5) {
		for (int k = 0; k < m_nombrePhases; k++) { pentesPhasesLocal1[k]->diviser(sommeCoeff); }
		pentesMelangeLocal1->diviser(sommeCoeff);
		for (int k = 0; k < m_nombreTransports; k++) { pentesTransportLocal1[k] /= sommeCoeff; }
	}

	//Rempli tampon pour envoi
	//------------------------
	for (int k = 0; k < m_nombrePhases; k++) {
		pentesPhasesLocal1[k]->rempliTamponPentes(tampon, compteur);
	}
	pentesMelangeLocal1->rempliTamponPentes(tampon, compteur);
	for (int k = 0; k < m_nombreTransports; k++) {
		tampon[++compteur] = pentesTransportLocal1[k];
	}
}

//***********************************************************************

void CelluleO2::rempliTamponPentesAMRjeSuisCpuGauche(double *tampon, int &compteur, const int &lvl) const
{
	if (m_lvl == lvl) {
		//Mise a zero de la pente locale a transmettre
		//--------------------------------------------
		double epsilon(1.e-8), sommeCoeff(0.);
		for (int k = 0; k < m_nombrePhases; k++) {
			pentesPhasesLocal1[k]->miseAZero();
		}
		pentesMelangeLocal1->miseAZero();
		for (int k = 0; k < m_nombreTransports; k++) {
			pentesTransportLocal1[k] = 0.;
		}

		//Boucle sur les bords pour la détermination de la pente a transmettre
		//--------------------------------------------------------------------
		for (unsigned int b = 0; b < m_bords.size(); b++) {
      if (!m_bords[b]->getSplit()) {
        if (fabs(m_bords[b]->getFace()->getNormale().getX()) > 0.99) {
          if ((m_bords[b]->getFace()->getPos().getX() + epsilon) < m_element->getPosition().getX()) {
            for (int k = 0; k < m_nombrePhases; k++) { pentesPhasesLocal1[k]->multiplieEtAjoute(*m_bords[b]->getPentesPhase(k), 1.); }
            pentesMelangeLocal1->multiplieEtAjoute(*m_bords[b]->getPentesMelange(), 1.);
            for (int k = 0; k < m_nombreTransports; k++) { pentesTransportLocal1[k] += (m_bords[b]->getPentesTransport(k)->getValeur()); }
            sommeCoeff += 1.;
          }
        }
      }
		}

		//Normalisation de la pente
		//-------------------------
		if (sommeCoeff > 1.e-5) {
			for (int k = 0; k < m_nombrePhases; k++) { pentesPhasesLocal1[k]->diviser(sommeCoeff); }
			pentesMelangeLocal1->diviser(sommeCoeff);
			for (int k = 0; k < m_nombreTransports; k++) { pentesTransportLocal1[k] /= sommeCoeff; }
		}

		//Rempli tampon pour envoi
		//------------------------
		for (int k = 0; k < m_nombrePhases; k++) {
			pentesPhasesLocal1[k]->rempliTamponPentes(tampon, compteur);
		}
		pentesMelangeLocal1->rempliTamponPentes(tampon, compteur);
		for (int k = 0; k < m_nombreTransports; k++) {
			tampon[++compteur] = pentesTransportLocal1[k];
		}
	}

	else {
		for (unsigned int i = 0; i < m_cellulesEnfants.size(); i++) {
			//Je suis CPU gauche, j'envoie donc la moyenne des pentes cote gauche des elements limites de droite
			if ((i % 2) == 1) { m_cellulesEnfants[i]->rempliTamponPentesAMRjeSuisCpuGauche(tampon, compteur, lvl); }
		}
	}
}

//***********************************************************************

void CelluleO2::rempliTamponPentesAMRjeSuisCpuDroite(double *tampon, int &compteur, const int &lvl) const
{
	if (m_lvl == lvl) {
		//Mise a zero de la pente locale a transmettre
		//--------------------------------------------
		double epsilon(1.e-8), sommeCoeff(0.);
		for (int k = 0; k < m_nombrePhases; k++) {
			pentesPhasesLocal1[k]->miseAZero();
		}
		pentesMelangeLocal1->miseAZero();
		for (int k = 0; k < m_nombreTransports; k++) {
			pentesTransportLocal1[k] = 0.;
		}

		//Boucle sur les bords pour la détermination de la pente a transmettre
		//--------------------------------------------------------------------
		for (unsigned int b = 0; b < m_bords.size(); b++) {
      if (!m_bords[b]->getSplit()) {
        if (fabs(m_bords[b]->getFace()->getNormale().getX()) > 0.99) {
          if ((m_bords[b]->getFace()->getPos().getX() - epsilon) > m_element->getPosition().getX()) {
            for (int k = 0; k < m_nombrePhases; k++) { pentesPhasesLocal1[k]->multiplieEtAjoute(*m_bords[b]->getPentesPhase(k), 1.); }
            pentesMelangeLocal1->multiplieEtAjoute(*m_bords[b]->getPentesMelange(), 1.);
            for (int k = 0; k < m_nombreTransports; k++) { pentesTransportLocal1[k] += (m_bords[b]->getPentesTransport(k)->getValeur()); }
            sommeCoeff += 1.;
          }
        }
      }
		}

		//Normalisation de la pente
		//-------------------------
		if (sommeCoeff > 1.e-5) {
			for (int k = 0; k < m_nombrePhases; k++) { pentesPhasesLocal1[k]->diviser(sommeCoeff); }
			pentesMelangeLocal1->diviser(sommeCoeff);
			for (int k = 0; k < m_nombreTransports; k++) { pentesTransportLocal1[k] /= sommeCoeff; }
		}

		//Rempli tampon pour envoi
		//------------------------
		for (int k = 0; k < m_nombrePhases; k++) {
			pentesPhasesLocal1[k]->rempliTamponPentes(tampon, compteur);
		}
		pentesMelangeLocal1->rempliTamponPentes(tampon, compteur);
		for (int k = 0; k < m_nombreTransports; k++) {
			tampon[++compteur] = pentesTransportLocal1[k];
		}
	}

	else {
		for (unsigned int i = 0; i < m_cellulesEnfants.size(); i++) {
			//Je suis CPU droite, j'envoie donc la moyenne des pentes cote droite des elements limites de gauche
			if ((i % 2) == 0) { m_cellulesEnfants[i]->rempliTamponPentesAMRjeSuisCpuDroite(tampon, compteur, lvl); }
		}
	}
}

//***********************************************************************
