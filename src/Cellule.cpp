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

#include "Cellule.h"

using namespace std;

//***********************************************************************

Cellule::Cellule() : m_vecPhases(0), m_melange(0), m_cons(0), m_vecTransports(0), m_consTransports(0), m_cellulesEnfants(0)
{
  m_lvl = 0;
  m_xi = 0.;
	m_split = false;
}

//***********************************************************************

Cellule::Cellule(int lvl) : m_vecPhases(0), m_melange(0), m_cons(0), m_vecTransports(0), m_consTransports(0), m_cellulesEnfants(0)
{
  m_lvl = lvl;
  m_xi = 0.;
	m_split = false;
}

//***********************************************************************

Cellule::~Cellule()
{
  for (int k = 0; k < m_nombrePhases; k++) {
    delete m_vecPhases[k];
  }
  delete[] m_vecPhases;
  delete m_melange;
  delete[] m_vecTransports;
  delete m_cons;
  delete[] m_consTransports;
  for (unsigned int gpa = 0; gpa < m_vecGrandeursPhysAdd.size(); gpa++) {
    delete m_vecGrandeursPhysAdd[gpa];
  }
  for (unsigned int i = 0; i < m_bordsEnfantsInternes.size(); i++) {
    m_bordsEnfantsInternes[i]->finaliseFace();
    delete m_bordsEnfantsInternes[i];
  }
  m_bordsEnfantsInternes.clear();
  for (unsigned int i = 0; i < m_cellulesEnfants.size(); i++) {
    delete m_cellulesEnfants[i];
  }
  m_cellulesEnfants.clear();
}

//***********************************************************************

void Cellule::ajouteBord(BordDeMaille *bord)
{
  m_bords.push_back(bord);
}

//***********************************************************************

void Cellule::supprimeBord(BordDeMaille *bord)
{
  for (unsigned int b = 0; b < m_bords.size(); b++) {
    if (m_bords[b] == bord) { m_bords.erase(m_bords.begin() + b); }
  }
}

//***********************************************************************

int Cellule::getBordsSize() const
{
  return m_bords.size();
}

//***********************************************************************

BordDeMaille* Cellule::getBord(int &b)
{
  return m_bords[b];
}

//***********************************************************************

void Cellule::alloue(const int &nombrePhases, const int &nombreTransports, const std::vector<PhysAdd*> &physAdd, Modele *modele)
{
  m_nombrePhases = nombrePhases;
  m_nombreTransports = nombreTransports;
  m_vecPhases = new Phase*[nombrePhases];
  for (int k = 0; k < nombrePhases; k++) {
    modele->allouePhase(&m_vecPhases[k]);
  }
  modele->alloueMelange(&m_melange);
  modele->alloueCons(&m_cons,nombrePhases);
  if (nombreTransports > 0) {
    m_vecTransports = new Transport[nombreTransports];
    m_consTransports = new Transport[nombreTransports];
  }
  for (unsigned int k = 0; k < physAdd.size(); k++) {
    physAdd[k]->ajouteGrandeurPhysAdd(this);
  }
}

//***********************************************************************

void Cellule::alloueEosYk(const int &nombrePhases, Modele *modele)
{
	modele->alloueEos(*this, nombrePhases);
	m_melange->alloueYk(nombrePhases);
}

//***********************************************************************

void Cellule::rempli(vector<DomaineGeometrique*> &domaines)
{
  Coord coordonnees;
  coordonnees = m_element->getPosition();
  for (unsigned int geom = 0; geom < domaines.size(); geom++) {
    if (domaines[geom]->appartient(coordonnees)) {
      domaines[geom]->rempli(this, m_nombrePhases, m_nombreTransports);
    }
  }
}

//***********************************************************************

void Cellule::alloueEtCopiePhase(const int &numeroPhase, Phase *phase)
{
  phase->alloueEtCopiePhase(&m_vecPhases[numeroPhase]);
}

//***********************************************************************

void Cellule::copiePhase(const int &numeroPhase, Phase *phase)
{
  m_vecPhases[numeroPhase]->copiePhase(*phase);
}

//***********************************************************************

void Cellule::copieMelange(Melange *melange)
{
  m_melange->copieMelange(*melange);
}

//***********************************************************************

void Cellule::miseAZeroCons(const int &nombrePhases, const int &nombreTransports)
{
  m_cons->miseAZero(nombrePhases);
  for (int k = 0; k < nombreTransports; k++) {
    m_consTransports[k].setValeur(0.);
  }
}

//***********************************************************************

void Cellule::miseAZeroConsGenerale(const int &nombrePhases, const int &nombreTransports)
{
  if (!m_split) {
    m_cons->miseAZero(nombrePhases);
    for (int k = 0; k < nombreTransports; k++) {
      m_consTransports[k].setValeur(0.);
    }
  }
  else {
    for (unsigned int i = 0; i < m_cellulesEnfants.size(); i++) {
      m_cellulesEnfants[i]->miseAZeroConsGenerale(nombrePhases, nombreTransports);
    }
  }
}

//***********************************************************************

void Cellule::miseAZeroFluxTemp(const int &nombrePhases)
{
  m_cons->miseAZeroFluxTemp(nombrePhases);
}

//***********************************************************************

void Cellule::evolutionTemporelle(const double &dt, const int &nombrePhases, const int &nombreTransports)
{
  m_cons->miseEnTampon(*this, nombrePhases);   //On determine Un grace au vecteur primitif des phases et de melange que l on stocke dans fluxTempXXX
  m_cons->multiplie(dt, nombrePhases);         //On multiplie m_cons par dt
  m_cons->ajoutFlux(1., nombrePhases);         //On y ajoute fluxTempXXX (Un) -> on obtient Un+1 dans m_cons

  //Idem pour transport (juste pas besoin de construire Un)
  for (int k = 0; k < nombreTransports; k++) {
    m_consTransports[k].multiplie(dt);
    m_vecTransports[k].ajoute(m_consTransports[k].getValeur());
  }
}

//***********************************************************************

void Cellule::construitPrim(const int &nombrePhases)
{
  m_cons->construitPrim(m_vecPhases, m_melange, nombrePhases);
}

//***********************************************************************

void Cellule::construitCons(const int &nombrePhases)
{
  m_cons->construitCons(m_vecPhases, nombrePhases, m_melange);
}

//***********************************************************************

void Cellule::relaxPressions(const int &nombrePhases)
{
  m_cons->relaxPressions(this, nombrePhases);
}

//***********************************************************************

void Cellule::relaxPTMu(const int &nombrePhases)
{
  m_cons->relaxPTMu(this, nombrePhases);
}

//***********************************************************************

void Cellule::correctionEnergie(const int &nombrePhases)
{
  m_melange->energieTotaleVersEnergieInterne(m_vecGrandeursPhysAdd); //On reconstruit l'energie interne a partir de l energie totale
  m_cons->correctionEnergie(this, nombrePhases);                     //On recalcule la pression
}

//***********************************************************************

void Cellule::ecritPhasesMelange(const int &nombrePhases, const int &nombreTransports, ofstream &fluxFichier) const
{
  for (int k = 0; k < nombrePhases; k++) { m_vecPhases[k]->ecritPhase(fluxFichier); }
  m_melange->ecritMelange(fluxFichier);
  for (int k = 0; k < nombreTransports; k++) { fluxFichier << m_vecTransports[k].getValeur() << " "; }
}

//***********************************************************************

void Cellule::calculsEtendus(const int &nombrePhases, Prim type)
{
  for (int k = 0; k < nombrePhases; k++) {
    m_vecPhases[k]->calculsEtendusPhases(m_melange->getVitesse());
  }
  this->preparePhysAdd();
  m_melange->calculGrandeursMelange(m_vecPhases, nombrePhases);
  m_melange->energieInterneVersEnergieTotale(m_vecGrandeursPhysAdd); //Mis a part car depend de la physique.
}

//***********************************************************************

void Cellule::calculsEtendusPourRiemann(const int &nombrePhases)
{
  for (int k = 0; k < nombrePhases; k++) {
    m_vecPhases[k]->calculsEtendusPhases(m_melange->getVitesse());
  }
  m_melange->calculGrandeursMelange(m_vecPhases, nombrePhases);
  m_melange->calculEnergieTotale();
}

//***********************************************************************

void Cellule::calculsEtendusPourCommunications(const int &nombrePhases, Prim type)
{
  for (int k = 0; k < nombrePhases; k++) {
    m_vecPhases[k]->calculsEtendusPhases(m_melange->getVitesse());
  }
  m_melange->calculGrandeursMelange(m_vecPhases, nombrePhases);
}

//***********************************************************************

void Cellule::projection(const Coord &normale, const Coord &tangente, const Coord &binormale, const int &nombrePhases, Prim type)
{
  for (int k = 0; k < nombrePhases; k++) {
    m_vecPhases[k]->projection(normale, tangente, binormale);
  }
  m_melange->projection(normale, tangente, binormale);
}

//***********************************************************************

void Cellule::projectionRepereAbsolu(const Coord &normale, const Coord &tangente, const Coord &binormale, const int &nombrePhases, Prim type)
{
  for (int k = 0; k < nombrePhases; k++) {
    m_vecPhases[k]->projectionRepereAbsolu(normale, tangente, binormale);
  }
  m_melange->projectionRepereAbsolu(normale, tangente, binormale);
}

//***********************************************************************

void Cellule::copieVec(Phase **vecPhases, Melange *melange, Transport *vecTransports)
{
  for (int k = 0; k < m_nombrePhases; k++) {
    m_vecPhases[k]->copiePhase(*vecPhases[k]);
  }
  m_melange->copieMelange(*melange);
  //On stocke l'energie capillaire dans l'energie totale ici, utile ensuite pour le calculsEtendusPourRiemann
  m_melange->calculEnergieCapillaire();
  for (int k = 0; k < m_nombreTransports; k++) {
    m_vecTransports[k] = vecTransports[k];
  }
}

//***********************************************************************

//void Cellule::ecritureCoupe1Dde2D(std::ofstream &fluxFichier, std::string variableConstanteCoupe, const double &valeurCoupe, const double &dL)
//{
//  if (m_cellulesEnfants.size() == 0) {
//    bool imprX(false), imprY(false), imprZ(false);
//    double dLsur2, position, epsilon;
//    dLsur2 = dL / pow(2., (double)m_lvl) / 2.;
//    epsilon = 1.e-3*dLsur2;
//    if (variableConstanteCoupe == "x") {
//      imprY = true;
//      position = m_element->getPosition().getX();
//    }
//    else if (variableConstanteCoupe == "y") {
//      imprX = true;
//      position = m_element->getPosition().getY();
//    }
//
//    //imprX = true;
//    //imprY = false;
//    //if (fabs(m_element->getPosition().getX() - m_element->getPosition().getY()) < 1.e-6) {
//    if (fabs(position - valeurCoupe - epsilon) < dLsur2) {
//      m_element->ecritPos(fluxFichier, imprX, imprY, imprZ);
//      this->ecritPhasesMelange(m_nombrePhases, m_nombreTransports, fluxFichier);
//      fluxFichier << m_lvl << " " << m_xi << " ";
//      fluxFichier << endl;
//    }
//  }
//  else {
//    for (unsigned int i = 0; i < m_cellulesEnfants.size(); i++) {
//      m_cellulesEnfants[i]->ecritureCoupe1Dde2D(fluxFichier, variableConstanteCoupe, valeurCoupe, dL);
//    }
//  }
//}
//
////***********************************************************************
//
//void Cellule::ecritureCoupe1Dde3D(std::ofstream &fluxFichier, std::string variableConstanteCoupe1, std::string variableConstanteCoupe2,
//  const double &valeurCoupe1, const double &valeurCoupe2, const double &dL1, const double &dL2)
//{
//  if (m_cellulesEnfants.size() == 0) {
//    bool imprX(true), imprY(true), imprZ(true);
//    double dL1sur2, dL2sur2, position1, position2, epsilon1, epsilon2;
//    dL1sur2 = dL1 / pow(2., (double)m_lvl) / 2.;
//    dL2sur2 = dL2 / pow(2., (double)m_lvl) / 2.;
//    epsilon1 = 1.e-3*dL1sur2;
//    epsilon2 = 1.e-3*dL2sur2;
//
//    if (variableConstanteCoupe1 == "x") {
//      imprX = false;
//      position1 = m_element->getPosition().getX();
//    }
//    else if (variableConstanteCoupe1 == "y") {
//      imprY = false;
//      position1 = m_element->getPosition().getY();
//    }
//    else {
//      imprZ = false;
//      position1 = m_element->getPosition().getZ();
//    }
//
//    if (variableConstanteCoupe2 == "x") {
//      imprX = false;
//      position2 = m_element->getPosition().getX();
//    }
//    else if (variableConstanteCoupe2 == "y") {
//      imprY = false;
//      position2 = m_element->getPosition().getY();
//    }
//    else {
//      imprZ = false;
//      position2 = m_element->getPosition().getZ();
//    }
//
//    if ((fabs(position1 - valeurCoupe1 - epsilon1) < dL1sur2) && (fabs(position2 - valeurCoupe2 - epsilon2) < dL2sur2)) {
//      m_element->ecritPos(fluxFichier, imprX, imprY, imprZ);
//      this->ecritPhasesMelange(m_nombrePhases, m_nombreTransports, fluxFichier);
//      fluxFichier << m_lvl << " " << m_xi << " ";
//      fluxFichier << endl;
//    }
//  }
//  else {
//    for (unsigned int i = 0; i < m_cellulesEnfants.size(); i++) {
//      m_cellulesEnfants[i]->ecritureCoupe1Dde3D(fluxFichier, variableConstanteCoupe1, variableConstanteCoupe2, valeurCoupe1, valeurCoupe2, dL1, dL2);
//    }
//  }
//}

//****************************************************************************
//***********************Physiques additionnelles*****************************
//****************************************************************************

void Cellule::preparePhysAdd()
{
  for (unsigned int gpa = 0; gpa < m_vecGrandeursPhysAdd.size(); gpa++) {
    m_vecGrandeursPhysAdd[gpa]->calculGrandeurs(this);
  }
}

//***********************************************************************

double Cellule::selectionneScalaire(string nomVariable, int num) const
{
  //Selection scalaire concerne
  if (nomVariable == "TR") {
    return m_vecTransports[num].getValeur();
  }
  else if (nomVariable == "P") {
    if (m_nombrePhases > 1) {
      return m_melange->getPression();
    }
    else {
      return m_vecPhases[num]->getPression();
    }
  }
  else if (nomVariable == "RHO") {
    if (m_nombrePhases > 1) {
      return m_melange->getDensite();
    }
    else {
      return m_vecPhases[num]->getDensite();
    }
  }
  else if (nomVariable == "ALPHA") {
    if (m_nombrePhases > 1) {
      return m_vecPhases[num]->getAlpha();
    }
    else { return 1.; }
  }
  else if (nomVariable == "u") {
    if (m_nombrePhases > 1) {
      return m_melange->getVitesse().getX();
    }
    else {
      return m_vecPhases[num]->getU();
    }
  }
  else if (nomVariable == "v") {
    if (m_nombrePhases > 1) {
      return m_melange->getVitesse().getY();
    }
    else {
      return m_vecPhases[num]->getV();
    }
  }
  else if (nomVariable == "w") {
    if (m_nombrePhases > 1) {
      return m_melange->getVitesse().getZ();
    }
    else {
      return m_vecPhases[num]->getW();
    }
  }
  else if (nomVariable == "T") {
    return m_vecPhases[num]->getTemperature();
  }
  //FP//TODO// faire GPA et Phases
  else { Erreurs::messageErreur("nomVariable inconnu dans selectionneScalaire (lie a GrandeursPhysAdd)"); return 0; }
}

//***********************************************************************

void Cellule::setScalaire(string nomVariable, const double &valeur, int num, int indice)
{
  //Selection vecteur concerne
  if (nomVariable == "TR") { //transport
    m_vecTransports[num].setValeur(valeur);
  }
  //FP//TODO// faire GPA et Phases
  else { Erreurs::messageErreur("nomVariable inconnu dans setScalaire (lie a GrandeursPhysAdd)"); }
}

//***********************************************************************

Coord Cellule::selectionneVecteur(string nomVecteur, int num, int indice) const
{
  //Selection vecteur concerne
  if (nomVecteur == "GPA") { //Physique additionelle
    return m_vecGrandeursPhysAdd[num]->getGrad(indice);
  }
  //FP//TODO// faire vecteur pour les Phases
  else { Erreurs::messageErreur("nomVecteur inconnu dans selectionneVecteur (lie a GrandeursPhysAdd)"); return 0; }
}

//***********************************************************************

void Cellule::setVecteur(std::string nomVecteur, const Coord &valeur, int num, int indice)
{
  //Selection vecteur concerne
  if (nomVecteur == "GPA") { //Physique additionelle
    m_vecGrandeursPhysAdd[num]->setGrad(valeur, indice);
  }
  //FP//TODO// faire vecteur pour les Phases
  else { Erreurs::messageErreur("nomVecteur inconnu dans setVecteur (lie a GrandeursPhysAdd)"); }
}

//***********************************************************************

//KS// A mettre a jour pour prendre en compte le cas ou il y a des murs

Coord Cellule::calculGradient(string nomVariable, int num)
{
  double sommeDistanceX = 0.;
  double sommeDistanceY = 0.;
  double sommeDistanceZ = 0.;
  Coord grad = 0.;             /*!< vecteur du gradient de la variable en question sur la cellule*/
  double gradBord(0.);         /*!< gradient de la variable en question sur un bord de maille dans la direction de la face*/
  Coord gradBordProjete = 0.;  /*!< vecteur du gradient de la variable en question sur un bord de maille dans le repere absolu*/
  Cellule *cellGauche(0);
  Cellule *cellDroite(0);
  for (unsigned int b = 0; b < m_bords.size(); b++) {
    if (!m_bords[b]->getSplit()) {
      if (m_bords[b]->quiSuisJe() == 0) //Bord de type BordDeMaille/O2
      {
        // Recuperation des valeurs de la variable en question a gauche et a droite pour chaque bord de maille
        // et calcul du gradient sur la normale a la face
        cellGauche = m_bords[b]->getCellGauche();
        cellDroite = m_bords[b]->getCellDroite();
        double cg = cellGauche->selectionneScalaire(nomVariable, num);
        double cd = cellDroite->selectionneScalaire(nomVariable, num);

        double distance(cellGauche->distance(cellDroite));
        gradBord = (cd - cg) / distance;

        // Projection du gradient du bord de maille dans le repere absolu
        Coord normale = m_bords[b]->getFace()->getNormale();
        gradBordProjete.setX(normale.getX()*gradBord);
        gradBordProjete.setY(normale.getY()*gradBord);
        gradBordProjete.setZ(normale.getZ()*gradBord);

        // Sommation de tous les bords de maille projetes dans le repere absolue ponderes par les distances dans chaque direction
        // pour obtenir le gradient de la cellule que l'on normalise apres avec la somme des distances.
        double distanceX(cellGauche->distanceX(cellDroite));
        double distanceY(cellGauche->distanceY(cellDroite));
        double distanceZ(cellGauche->distanceZ(cellDroite));
        distanceX = fabs(distanceX);
        distanceY = fabs(distanceY);
        distanceZ = fabs(distanceZ);

        gradBordProjete.setXYZ(gradBordProjete.getX()*distanceX, gradBordProjete.getY()*distanceY, gradBordProjete.getZ()*distanceZ);

        sommeDistanceX += distanceX;
        sommeDistanceY += distanceY;
        sommeDistanceZ += distanceZ;

        grad += gradBordProjete;
      }
      else if (m_bords[b]->quiSuisJe() == 2) { //Bord de type Mur, prise en compte seulement pour symetrie
        double distanceX(this->distanceX(m_bords[b]));
        double distanceY(this->distanceY(m_bords[b]));
        double distanceZ(this->distanceZ(m_bords[b]));
        distanceX = fabs(distanceX)*2.;
        distanceY = fabs(distanceY)*2.;
        distanceZ = fabs(distanceZ)*2.;
        sommeDistanceX += distanceX;
        sommeDistanceY += distanceY;
        sommeDistanceZ += distanceZ;
      }
    }
  }

  // Verification du multi-D pour eviter division par zero
  if (sommeDistanceX <= 1.e-12) { sommeDistanceX = 1.; }
  if (sommeDistanceY <= 1.e-12) { sommeDistanceY = 1.; }
  if (sommeDistanceZ <= 1.e-12) { sommeDistanceZ = 1.; }

  // Gradient finale de la variable en question sur la cellule (on normalise)
  grad.setXYZ(grad.getX() / sommeDistanceX, grad.getY() / sommeDistanceY, grad.getZ() / sommeDistanceZ);

  return grad;
}

//***********************************************************************

GrandeursPhysAdd* Cellule::getGPA(int &numGPA) const
{
  return m_vecGrandeursPhysAdd[numGPA];
}

//***********************************************************************

Coord Cellule::getGradTk(int &numPhase, int &numPhysAdd) const
{
  return m_vecGrandeursPhysAdd[numPhysAdd]->getGradTk(numPhase);
}

//***********************************************************************

void Cellule::setGradTk(int &numPhase, int &numPhysAdd, double *tampon, int &compteur)
{
  Coord grad;
  grad.setX(tampon[++compteur]);
  grad.setY(tampon[++compteur]);
  grad.setZ(tampon[++compteur]);
  m_vecGrandeursPhysAdd[numPhysAdd]->setGradTk(numPhase, grad);
}

//***********************************************************************

void Cellule::ajoutNonConsPhysAdd(const int &nombrePhases, PhysAdd &physAdd)
{
  physAdd.ajoutNonConsPhysAdd(this, nombrePhases);
}

//***********************************************************************

void Cellule::reinitialiseFonctionCouleur()
{
	m_vecTransports[0].setValeur(m_vecPhases[1]->getAlpha()); //KS// A generaliser ...
}

//****************************************************************************
//*****************************Accesseurs*************************************
//****************************************************************************

Phase* Cellule::getPhase(const int &numeroPhase, Prim type) const
{
  return m_vecPhases[numeroPhase];
}

//***********************************************************************

Phase** Cellule::getPhases(Prim type) const
{
  return m_vecPhases;
}

//***********************************************************************

Melange* Cellule::getMelange(Prim type) const
{
  return m_melange;
}

//***********************************************************************

Flux* Cellule::getCons() const
{
  return m_cons;
}

//***********************************************************************

void Cellule::setCons(Flux *cons)
{
  m_cons->setCons(cons, m_nombrePhases);
}

//***********************************************************************

Coord Cellule::getPosition() const
{
  return m_element->getPosition();
}

//***********************************************************************

void Cellule::setElement(Element *element, const int &numCellule)
{
  m_element = element;
  m_element->setCelluleAssociee(numCellule);
}

//***********************************************************************

Element* Cellule::getElement()
{
  return m_element;
}

//***********************************************************************

void Cellule::setTransport(double valeur, int &numTransport, Prim type)
{
  m_vecTransports[numTransport].setValeur(valeur);
}

//***********************************************************************

Transport& Cellule::getTransport(const int &numTransport, Prim type) const
{
	return m_vecTransports[numTransport];
}

//***********************************************************************

Transport* Cellule::getTransports(Prim type) const
{
	return m_vecTransports;
}

//***********************************************************************

Transport* Cellule::getConsTransport(const int &numTransport) const
{
  return &m_consTransports[numTransport];
}

//***********************************************************************

void Cellule::setConsTransport(double valeur, const int &numTransport)
{
  m_consTransports[numTransport].setValeur(valeur);
}

//***********************************************************************

int Cellule::getNombrePhases() const
{
  return m_nombrePhases;
}

//***********************************************************************

int Cellule::getNombreTransports() const
{
  return m_nombreTransports;
}

//***********************************************************************

double Cellule::getGradient()
{
  string nomVariable = "RHO";
  Coord gradient = 0.;
  int var = 0; //Utile seulement pour monophase
  gradient = this->calculGradient(nomVariable, var);
  return gradient.norme();
}

//***********************************************************************

vector<GrandeursPhysAdd*>& Cellule::getVecGrandeursPhysAdd()
{
  return m_vecGrandeursPhysAdd;
}

//***********************************************************************

void Cellule::afficheInfos() const
{
  m_element->afficheInfos();
}

//****************************************************************************
//******************************Distances*************************************
//****************************************************************************

double Cellule::distance(Cellule *c)
{
  return m_element->distance(c->getElement());
}

//***********************************************************************

double Cellule::distanceX(Cellule *c)
{
	return m_element->distanceX(c->getElement());
}

//***********************************************************************

double Cellule::distanceY(Cellule *c)
{
	return m_element->distanceY(c->getElement());
}

//***********************************************************************

double Cellule::distanceZ(Cellule *c)
{
	return m_element->distanceZ(c->getElement());
}

//***********************************************************************

double Cellule::distance(BordDeMaille *b)
{
  return m_element->distance(b->getFace());
}

//***********************************************************************

double Cellule::distanceX(BordDeMaille *b)
{
  return m_element->distanceX(b->getFace());
}

//***********************************************************************

double Cellule::distanceY(BordDeMaille *b)
{
  return m_element->distanceY(b->getFace());
}

//***********************************************************************

double Cellule::distanceZ(BordDeMaille *b)
{
  return m_element->distanceZ(b->getFace());
}

//***********************************************************************

bool Cellule::traverseObjet(const ObjetGeometrique &objet) const
{
  return m_element->traverseObjet(objet);
}

//****************************************************************************
//***************************** Methode AMR **********************************
//****************************************************************************

void Cellule::miseAZeroXi()
{
  m_xi = 0.;
}

//***********************************************************************

void Cellule::miseAZeroConsXi()
{
  m_consXi = 0.;
}

//***********************************************************************

void Cellule::evolutionTemporelleXi(const double &dtDiff)
{
  m_xi += dtDiff*m_consXi;
}

//***********************************************************************

void Cellule::choixRaffine(const double &xiSplit, const int &nbMaillesY, const int &nbMaillesZ,
  const double &dX, const double &dY, const double &dZ, const std::vector<PhysAdd*> &physAdd, Modele *modele, int &nbMaillesTotalAMR)
{
  if (!m_split) {
    if (m_xi >= xiSplit) {
      if (!this->lvlVoisinTropPetit()) {
        this->raffineCelluleEtBords(nbMaillesY, nbMaillesZ, dX, dY, dZ, physAdd, modele);
        nbMaillesTotalAMR += m_cellulesEnfants.size() - 1;
      }
    }
  }
}

//***********************************************************************

void Cellule::choixDeraffine(const double &xiJoin, int &nbMaillesTotalAMR)
{
  if (m_split) {
    bool deraffineGlobal(false);
    if (m_xi < xiJoin) {
      deraffineGlobal = true;
      for (unsigned int i = 0; i < m_cellulesEnfants.size(); i++) {
        //Si un seul des mes enfants a lui meme des enfants, alors je ne peux pas deraffiner
        if (m_cellulesEnfants[i]->getNombreCellulesEnfants() > 0) { deraffineGlobal = false; }
      }
      //Si un de mes voisins a un niveau trop grand, alors en deraffinant je creerais un ecart de niveau de cellules superieur a 1, donc je ne deraffine pas
      if (deraffineGlobal) { if (this->lvlVoisinTropGrand()) { deraffineGlobal = false; }; }
    }
    //Si tous mes enfants peuvent etre deraffiner, je deraffine
    if (deraffineGlobal) {
      nbMaillesTotalAMR -= m_cellulesEnfants.size() - 1;
      this->deraffineCelluleEtBords();
    }
  }
}

//***********************************************************************

void Cellule::raffineCelluleEtBords(const int &nbMaillesY, const int &nbMaillesZ, const double &dX, const double &dY, const double &dZ,
  const vector<PhysAdd*> &physAdd, Modele *modele)
{
	m_split = true;

  //--------------------------------------
  //Initialisations (enfants et dimension)
  //--------------------------------------

  double dimX(1.), dimY(0.), dimZ(0.);
  int nombreCellulesEnfants(2);
  int dim(1);
  if (nbMaillesZ != 1) {
    nombreCellulesEnfants = 8;
    dimY = 1.;
    dimZ = 1.;
    dim = 3;
  }
  else if (nbMaillesY != 1) {
    nombreCellulesEnfants = 4;
    dimY = 1.;
    dim = 2;
  }
  BordDeMaille* bordRef(0);
  for (unsigned int b = 0; b < m_bords.size(); b++) {
    if (m_bords[b]->quiSuisJe() == 0) { bordRef = m_bords[b]; break; } //Bord de type BordDeMaille/O2
  }
  double surfaceEnfant(bordRef->getFace()->getSurface() * pow(2., (double)((dim - 1)*(bordRef->getLvl() - m_lvl - 1))));
  int allouePenteLocal = 1;
  
  //------------------------
  //Raffinement des cellules
  //------------------------

  //Initialisation des donnees de maillage pour les cellules enfants
  //----------------------------------------------------------------
  double posXCelluleParent, posYCelluleParent, posZCelluleParent;
  double dXParent, dYParent, dZParent;
  double posXEnfant, posYEnfant, posZEnfant;
  posXCelluleParent = m_element->getPosition().getX();
  posYCelluleParent = m_element->getPosition().getY();
  posZCelluleParent = m_element->getPosition().getZ();
  dXParent = dX / pow(2., (double)m_lvl);
  dYParent = dY / pow(2., (double)m_lvl);
  dZParent = dZ / pow(2., (double)m_lvl);
  double volumeCelluleParent, lCFLCelluleParent;
  volumeCelluleParent = m_element->getVolume();
  lCFLCelluleParent = m_element->getLCFL();

  for (int i = 0; i < nombreCellulesEnfants; i++) {

    //Creation des cellules enfants
    //-----------------------------
    this->creerCelluleEnfant(i, m_lvl);
    m_element->creerElementEnfant();
    m_cellulesEnfants[i]->setElement(m_element->getElementEnfant(i), i);
    m_cellulesEnfants[i]->getElement()->setVolume(volumeCelluleParent / (double)nombreCellulesEnfants);
    m_cellulesEnfants[i]->getElement()->setLCFL(lCFLCelluleParent / 2.);
    posXEnfant = posXCelluleParent + dimX*dXParent*(double)(-0.25 + 0.5 * (i % 2));
    posYEnfant = posYCelluleParent + dimY*dYParent*(double)(-0.25 + 0.5 * ((i / 2) % 2));
    posZEnfant = posZCelluleParent + dimZ*dZParent*(double)(-0.25 + 0.5 * ((i / 4) % 2));
    m_cellulesEnfants[i]->getElement()->setPos(posXEnfant, posYEnfant, posZEnfant);

    //Initialisation des tableaux principaux selon modele et nombre de phases
    //et initialisation Physique : Remplissage donnees physiques, pour les cellules enfants
    //-------------------------------------------------------------------------------------
    m_cellulesEnfants[i]->alloue(m_nombrePhases, m_nombreTransports, physAdd, modele);
    for (int k = 0; k < m_nombrePhases; k++) {
      m_cellulesEnfants[i]->copiePhase(k, m_vecPhases[k]);
    }
    m_cellulesEnfants[i]->copieMelange(m_melange);
    m_cellulesEnfants[i]->getCons()->miseAZero(m_nombrePhases);
    for (int k = 0; k < m_nombreTransports; k++) { m_cellulesEnfants[i]->setTransport(m_vecTransports[k].getValeur(), k); }
    for (int k = 0; k < m_nombreTransports; k++) { m_cellulesEnfants[i]->setConsTransport(0., k); }
    m_cellulesEnfants[i]->setXi(m_xi);
  }

  //------------------------------
  //Raffinement des bords internes
  //------------------------------

  if (nbMaillesZ == 1) {
    if (nbMaillesY == 1) {

      //Cas 1D
      //------

      // |-------------|-------------|
      // 1      0      0      1      2

      //Bord enfant interne numero 0 (face selon X)
      //-------------------------------------------
      bordRef->creerBordEnfantInterne(m_lvl, &m_bordsEnfantsInternes);
      m_bordsEnfantsInternes[0]->creerFaceEnfant(bordRef);
      m_bordsEnfantsInternes[0]->getFace()->setSurface(surfaceEnfant);
      m_bordsEnfantsInternes[0]->getFace()->setNormale(1., 0., 0.);
      m_bordsEnfantsInternes[0]->getFace()->setTangente(0., 1., 0.);
      m_bordsEnfantsInternes[0]->getFace()->setBinormale(0., 0., 1.);
      m_bordsEnfantsInternes[0]->getFace()->setPos(posXCelluleParent, posYCelluleParent, posZCelluleParent);
      m_bordsEnfantsInternes[0]->initialiseGauche(m_cellulesEnfants[0]);
      m_bordsEnfantsInternes[0]->initialiseDroite(m_cellulesEnfants[1]);
      m_cellulesEnfants[0]->ajouteBord(m_bordsEnfantsInternes[0]);
      m_cellulesEnfants[1]->ajouteBord(m_bordsEnfantsInternes[0]);

      //Attribution modele et pentes
      //----------------------------
      m_bordsEnfantsInternes[0]->associeModele(modele);
      m_bordsEnfantsInternes[0]->allouePentes(m_nombrePhases, m_nombreTransports, allouePenteLocal);
    }
    else {

      //Cas 2D
      //------

      // |--------------|--------------|
      // |              |              |
      // |              |              |
      // |      2      1|      3       |
      // |              |              |
      // |              |              |
      // |--------------|--------------|
      // |      2       |      3       |
      // |              |              |
      // |      0      0|      1       |
      // |              |              |
      // |              |              |
      // |--------------|--------------|

      for (int i = 0; i < 4; i++) {
        bordRef->creerBordEnfantInterne(m_lvl, &m_bordsEnfantsInternes);
        m_bordsEnfantsInternes[i]->creerFaceEnfant(bordRef);
        m_bordsEnfantsInternes[i]->getFace()->setSurface(surfaceEnfant);
        if (i < 2) {
          //Bords interne 0 et 1 (face selon X)
          //-----------------------------------
          m_bordsEnfantsInternes[i]->getFace()->setNormale(1., 0., 0.);
          m_bordsEnfantsInternes[i]->getFace()->setTangente(0., 1., 0.);
          m_bordsEnfantsInternes[i]->getFace()->setBinormale(0., 0., 1.);
          m_bordsEnfantsInternes[i]->getFace()->setPos(posXCelluleParent, posYCelluleParent + dYParent*(-0.25 + 0.5 * (double)i), posZCelluleParent);
          m_bordsEnfantsInternes[i]->initialiseGauche(m_cellulesEnfants[2 * i]);
          m_bordsEnfantsInternes[i]->initialiseDroite(m_cellulesEnfants[1 + 2 * i]);
          m_cellulesEnfants[2 * i]->ajouteBord(m_bordsEnfantsInternes[i]);
          m_cellulesEnfants[1 + 2 * i]->ajouteBord(m_bordsEnfantsInternes[i]);
        }
        else {
          //Bords interne 2 et 3 (face selon Y)
          //-----------------------------------
          m_bordsEnfantsInternes[i]->getFace()->setNormale(0., 1., 0.);
          m_bordsEnfantsInternes[i]->getFace()->setTangente(-1., 0., 0.);
          m_bordsEnfantsInternes[i]->getFace()->setBinormale(0., 0., 1.);
          m_bordsEnfantsInternes[i]->getFace()->setPos(posXCelluleParent + dXParent*(-0.25 + 0.5 * (double)(i % 2)), posYCelluleParent, posZCelluleParent);
          m_bordsEnfantsInternes[i]->initialiseGauche(m_cellulesEnfants[i % 2]);
          m_bordsEnfantsInternes[i]->initialiseDroite(m_cellulesEnfants[2 + i % 2]);
          m_cellulesEnfants[i % 2]->ajouteBord(m_bordsEnfantsInternes[i]);
          m_cellulesEnfants[2 + i % 2]->ajouteBord(m_bordsEnfantsInternes[i]);
        }
        //Attribution modele et pentes
        //----------------------------
        m_bordsEnfantsInternes[i]->associeModele(modele);
        m_bordsEnfantsInternes[i]->allouePentes(m_nombrePhases, m_nombreTransports, allouePenteLocal);
      }
    }
  }
  else {

    //Cas 3D
    //------

    //3 fois le plan 2D selon X, Y puis Z avec la numerotation suivante :
    //- Numero au centre correspond ici au numero de la face (normale vers nous),
    //- Numero dans le coin bas gauche pour cellules devant le plan (cellule droite),
    //- Numero dans le coin haut droite pour cellules derriere le plan (cellule gauche).

    //             Selon X :                          Selon Y:                           Selon Z :
    //
    // |--------------|--------------|    |--------------|--------------|    |--------------|--------------|
    // |             6|             2|    |             4|             0|    |             2|             3|
    // |              |              |    |              |              |    |              |              |
    // |      2       |      3       |    |      2       |      3       |    |      2       |      3       |
    // |              |              |    |              |              |    |              |              |
    // |7             |3             |    |6             |2             |    |6             |7             |
    // |--------------|--------------|    |--------------|--------------|    |--------------|--------------|
    // |             4|             0|    |             5|             1|    |             0|             1|
    // |              |              |    |              |              |    |              |              |
    // |      0       |      1       |    |      0       |      1       |    |      0       |      1       |
    // |              |              |    |              |              |    |              |              |
    // |5             |1             |    |7             |3             |    |4             |5             |
    // |--------------|--------------|    |--------------|--------------|    |--------------|--------------|
    //
    //                y                              z---|                                  y
    //                |                                  |                                  |
    //            z---|                                  x                                  |---x

    //Face selon X
    for (int i = 0; i < 4; i++) {
      bordRef->creerBordEnfantInterne(m_lvl, &m_bordsEnfantsInternes);
      m_bordsEnfantsInternes[i]->creerFaceEnfant(bordRef);
      m_bordsEnfantsInternes[i]->getFace()->setSurface(surfaceEnfant);
      m_bordsEnfantsInternes[i]->getFace()->setNormale(1., 0., 0.);
      m_bordsEnfantsInternes[i]->getFace()->setTangente(0., 1., 0.);
      m_bordsEnfantsInternes[i]->getFace()->setBinormale(0., 0., 1.);
      if (i == 0) {
        m_bordsEnfantsInternes[i]->getFace()->setPos(posXCelluleParent, posYCelluleParent - 0.25*dYParent, posZCelluleParent + 0.25*dZParent);
        m_bordsEnfantsInternes[i]->initialiseGauche(m_cellulesEnfants[4]);
        m_bordsEnfantsInternes[i]->initialiseDroite(m_cellulesEnfants[5]);
        m_cellulesEnfants[4]->ajouteBord(m_bordsEnfantsInternes[i]);
        m_cellulesEnfants[5]->ajouteBord(m_bordsEnfantsInternes[i]);
      }
      else if (i == 1) {
        m_bordsEnfantsInternes[i]->getFace()->setPos(posXCelluleParent, posYCelluleParent - 0.25*dYParent, posZCelluleParent - 0.25*dZParent);
        m_bordsEnfantsInternes[i]->initialiseGauche(m_cellulesEnfants[0]);
        m_bordsEnfantsInternes[i]->initialiseDroite(m_cellulesEnfants[1]);
        m_cellulesEnfants[0]->ajouteBord(m_bordsEnfantsInternes[i]);
        m_cellulesEnfants[1]->ajouteBord(m_bordsEnfantsInternes[i]);
      }
      else if (i == 2) {
        m_bordsEnfantsInternes[i]->getFace()->setPos(posXCelluleParent, posYCelluleParent + 0.25*dYParent, posZCelluleParent + 0.25*dZParent);
        m_bordsEnfantsInternes[i]->initialiseGauche(m_cellulesEnfants[6]);
        m_bordsEnfantsInternes[i]->initialiseDroite(m_cellulesEnfants[7]);
        m_cellulesEnfants[6]->ajouteBord(m_bordsEnfantsInternes[i]);
        m_cellulesEnfants[7]->ajouteBord(m_bordsEnfantsInternes[i]);
      }
      else {
        m_bordsEnfantsInternes[i]->getFace()->setPos(posXCelluleParent, posYCelluleParent + 0.25*dYParent, posZCelluleParent - 0.25*dZParent);
        m_bordsEnfantsInternes[i]->initialiseGauche(m_cellulesEnfants[2]);
        m_bordsEnfantsInternes[i]->initialiseDroite(m_cellulesEnfants[3]);
        m_cellulesEnfants[2]->ajouteBord(m_bordsEnfantsInternes[i]);
        m_cellulesEnfants[3]->ajouteBord(m_bordsEnfantsInternes[i]);
      }
      //Attribution modele et pentes
      m_bordsEnfantsInternes[i]->associeModele(modele);
      m_bordsEnfantsInternes[i]->allouePentes(m_nombrePhases, m_nombreTransports, allouePenteLocal);
    }

    //Face selon Y
    for (int i = 4; i < 8; i++) {
      bordRef->creerBordEnfantInterne(m_lvl, &m_bordsEnfantsInternes);
      m_bordsEnfantsInternes[i]->creerFaceEnfant(bordRef);
      m_bordsEnfantsInternes[i]->getFace()->setSurface(surfaceEnfant);
      m_bordsEnfantsInternes[i]->getFace()->setNormale(0., 1., 0.);
      m_bordsEnfantsInternes[i]->getFace()->setTangente(-1., 0., 0.);
      m_bordsEnfantsInternes[i]->getFace()->setBinormale(0., 0., 1.);
      if (i == 4) {
        m_bordsEnfantsInternes[i]->getFace()->setPos(posXCelluleParent + 0.25*dXParent, posYCelluleParent, posZCelluleParent + 0.25*dZParent);
        m_bordsEnfantsInternes[i]->initialiseGauche(m_cellulesEnfants[5]);
        m_bordsEnfantsInternes[i]->initialiseDroite(m_cellulesEnfants[7]);
        m_cellulesEnfants[5]->ajouteBord(m_bordsEnfantsInternes[i]);
        m_cellulesEnfants[7]->ajouteBord(m_bordsEnfantsInternes[i]);
      }
      else if (i == 5) {
        m_bordsEnfantsInternes[i]->getFace()->setPos(posXCelluleParent + 0.25*dXParent, posYCelluleParent, posZCelluleParent - 0.25*dZParent);
        m_bordsEnfantsInternes[i]->initialiseGauche(m_cellulesEnfants[1]);
        m_bordsEnfantsInternes[i]->initialiseDroite(m_cellulesEnfants[3]);
        m_cellulesEnfants[1]->ajouteBord(m_bordsEnfantsInternes[i]);
        m_cellulesEnfants[3]->ajouteBord(m_bordsEnfantsInternes[i]);
      }
      else if (i == 6) {
        m_bordsEnfantsInternes[i]->getFace()->setPos(posXCelluleParent - 0.25*dXParent, posYCelluleParent, posZCelluleParent + 0.25*dZParent);
        m_bordsEnfantsInternes[i]->initialiseGauche(m_cellulesEnfants[4]);
        m_bordsEnfantsInternes[i]->initialiseDroite(m_cellulesEnfants[6]);
        m_cellulesEnfants[4]->ajouteBord(m_bordsEnfantsInternes[i]);
        m_cellulesEnfants[6]->ajouteBord(m_bordsEnfantsInternes[i]);
      }
      else {
        m_bordsEnfantsInternes[i]->getFace()->setPos(posXCelluleParent - 0.25*dXParent, posYCelluleParent, posZCelluleParent - 0.25*dZParent);
        m_bordsEnfantsInternes[i]->initialiseGauche(m_cellulesEnfants[0]);
        m_bordsEnfantsInternes[i]->initialiseDroite(m_cellulesEnfants[2]);
        m_cellulesEnfants[0]->ajouteBord(m_bordsEnfantsInternes[i]);
        m_cellulesEnfants[2]->ajouteBord(m_bordsEnfantsInternes[i]);
      }
      //Attribution modele et pentes
      m_bordsEnfantsInternes[i]->associeModele(modele);
      m_bordsEnfantsInternes[i]->allouePentes(m_nombrePhases, m_nombreTransports, allouePenteLocal);
    }

    //Face selon Z
    for (int i = 8; i < 12; i++) {
      bordRef->creerBordEnfantInterne(m_lvl, &m_bordsEnfantsInternes);
      m_bordsEnfantsInternes[i]->creerFaceEnfant(bordRef);
      m_bordsEnfantsInternes[i]->getFace()->setSurface(surfaceEnfant);
      m_bordsEnfantsInternes[i]->getFace()->setNormale(0., 0., 1.);
      m_bordsEnfantsInternes[i]->getFace()->setTangente(1., 0., 0.);
      m_bordsEnfantsInternes[i]->getFace()->setBinormale(0., 1., 0.);
      if (i == 8) {
        m_bordsEnfantsInternes[i]->getFace()->setPos(posXCelluleParent - 0.25*dXParent, posYCelluleParent - 0.25*dYParent, posZCelluleParent);
        m_bordsEnfantsInternes[i]->initialiseGauche(m_cellulesEnfants[0]);
        m_bordsEnfantsInternes[i]->initialiseDroite(m_cellulesEnfants[4]);
        m_cellulesEnfants[0]->ajouteBord(m_bordsEnfantsInternes[i]);
        m_cellulesEnfants[4]->ajouteBord(m_bordsEnfantsInternes[i]);
      }
      else if (i == 9) {
        m_bordsEnfantsInternes[i]->getFace()->setPos(posXCelluleParent + 0.25*dXParent, posYCelluleParent - 0.25*dYParent, posZCelluleParent);
        m_bordsEnfantsInternes[i]->initialiseGauche(m_cellulesEnfants[1]);
        m_bordsEnfantsInternes[i]->initialiseDroite(m_cellulesEnfants[5]);
        m_cellulesEnfants[1]->ajouteBord(m_bordsEnfantsInternes[i]);
        m_cellulesEnfants[5]->ajouteBord(m_bordsEnfantsInternes[i]);
      }
      else if (i == 10) {
        m_bordsEnfantsInternes[i]->getFace()->setPos(posXCelluleParent - 0.25*dXParent, posYCelluleParent + 0.25*dYParent, posZCelluleParent);
        m_bordsEnfantsInternes[i]->initialiseGauche(m_cellulesEnfants[2]);
        m_bordsEnfantsInternes[i]->initialiseDroite(m_cellulesEnfants[6]);
        m_cellulesEnfants[2]->ajouteBord(m_bordsEnfantsInternes[i]);
        m_cellulesEnfants[6]->ajouteBord(m_bordsEnfantsInternes[i]);
      }
      else {
        m_bordsEnfantsInternes[i]->getFace()->setPos(posXCelluleParent + 0.25*dXParent, posYCelluleParent + 0.25*dYParent, posZCelluleParent);
        m_bordsEnfantsInternes[i]->initialiseGauche(m_cellulesEnfants[3]);
        m_bordsEnfantsInternes[i]->initialiseDroite(m_cellulesEnfants[7]);
        m_cellulesEnfants[3]->ajouteBord(m_bordsEnfantsInternes[i]);
        m_cellulesEnfants[7]->ajouteBord(m_bordsEnfantsInternes[i]);
      }
      //Attribution modele et pentes
      m_bordsEnfantsInternes[i]->associeModele(modele);
      m_bordsEnfantsInternes[i]->allouePentes(m_nombrePhases, m_nombreTransports, allouePenteLocal);
    }
  }

  //------------------------------
  //Raffinement des bords externes
  //------------------------------

  for (unsigned int b = 0; b < m_bords.size(); b++) {
    if (!m_bords[b]->getSplit()) { m_bords[b]->raffineBordExterne(nbMaillesY, nbMaillesZ, dXParent, dYParent, dZParent, this, surfaceEnfant); }
  }
}

//***********************************************************************

void Cellule::creerCelluleEnfant(const int &num, const int &lvl)
{
  m_cellulesEnfants.push_back(new Cellule(lvl + 1));
}

//***********************************************************************

void Cellule::deraffineCelluleEtBords()
{
  //--------------------------------
  //Mise a jour de la cellule parent
  //--------------------------------

  this->moyenneEnfantsDansParent();

  //--------------------------------------
  //Suppression des bords enfants internes
  //--------------------------------------

  for (unsigned int i = 0; i < m_bordsEnfantsInternes.size(); i++) {
    m_bordsEnfantsInternes[i]->finaliseFace();
    delete m_bordsEnfantsInternes[i];
  }
  m_bordsEnfantsInternes.clear();

  //--------------------------------------
  //Suppression des bords enfants externes
  //--------------------------------------

  for (unsigned int b = 0; b < m_bords.size(); b++) {
    m_bords[b]->deraffineBordExterne(this);
  }

  //--------------------------------
  //Suppression des cellules enfants
  //--------------------------------

  for (unsigned int i = 0; i < m_cellulesEnfants.size(); i++) {
    delete m_cellulesEnfants[i];
  }
  m_cellulesEnfants.clear();
  m_element->finaliseElementsEnfants();

	m_split = false;
}

//***********************************************************************

void Cellule::moyenneEnfantsDansParent()
{
  int nombreCellulesEnfants(m_cellulesEnfants.size());
  if (nombreCellulesEnfants > 0) {
    //Reconstruction thermodynamique de la cellule parent
    //On commence par moyenner les variables conservatives
    m_cons->miseAZero(m_nombrePhases);
    for (int i = 0; i < nombreCellulesEnfants; i++) {
      m_cons->miseEnTampon(*m_cellulesEnfants[i], m_nombrePhases);              //On construit les variables conservatives des cellules enfants dans fluxTempXXX
      m_cons->ajoutFlux(1., m_nombrePhases);                                    //On les ajoute dans m_cons
    }
    m_cons->multiplie(1. / (double)nombreCellulesEnfants, m_nombrePhases);      //On divise m_cons par le nombre d enfants pour obtenir la moyenne
    //Puis ensuite on reconstruit les variables primitives en couplant avec une relaxation
    m_cons->construitPrim(m_vecPhases, m_melange, m_nombrePhases);
    m_cons->relaxPressions(this, m_nombrePhases);

    //Transport
    for (int k = 0; k < m_nombreTransports; k++) {
      double transport(0.);
      for (int i = 0; i < nombreCellulesEnfants; i++) {
        transport += m_cellulesEnfants[i]->getTransport(k).getValeur();
      }
      transport = transport / (double)nombreCellulesEnfants;
      m_vecTransports[k].setValeur(transport);
    }

    //Mise a zero de m_cons pour la suite
    m_cons->miseAZero(m_nombrePhases);
    for (int k = 0; k < m_nombreTransports; k++) {
      m_consTransports[k].setValeur(0.);
    }
  }
}

//***********************************************************************

bool Cellule::lvlVoisinTropGrand()
{
  bool critere(false);
  for (unsigned int b = 0; b < m_bords.size(); b++) {
    if (m_bords[b]->getLvl() == m_lvl) {
      for (int bEnfant = 0; bEnfant < m_bords[b]->getNombreBordsEnfants(); bEnfant++) {
        if (m_bords[b]->getBordEnfant(bEnfant)->getSplit()) { critere = true; return critere; }
      }
    }
    else {
      if (m_bords[b]->getSplit()) { critere = true; return critere; }
    }
  }
  return critere;
}

//***********************************************************************

bool Cellule::lvlVoisinTropPetit()
{
  bool critere(false);
  for (unsigned int b = 0; b < m_bords.size(); b++) {
    if (!m_bords[b]->getSplit()) {
      if (m_bords[b]->quiSuisJe() == 0) //Bord de type BordDeMaille/O2
      {
        // Recuperation des niveaux AMR a gauche et a droite pour chaque bord de maille
        int lvlg = m_bords[b]->getCellGauche()->getLvl();
        int lvld = m_bords[b]->getCellDroite()->getLvl();

        // Regarde si le voisin en question a un niveau AMR trop petit
        if (lvlg < m_lvl) { critere = true; return critere; }
        if (lvld < m_lvl) { critere = true; return critere; }
      }
      else
      {
        // Recuperation des niveaux AMR a gauche seulement pour chaque condition aux limites
        int lvlg = m_bords[b]->getCellGauche()->getLvl();

        // Regarde si le voisin en question a un niveau AMR trop petit
        if (lvlg < m_lvl) { critere = true; return critere; }
      }
    }
  }
  return critere;
}

//***********************************************************************

void Cellule::constructionTableauxCellulesLvlEtBordsInternesLvl(vector<Cellule *> *cellulesLvl, vector<BordDeMaille *> *bordsLvl)
{
  for (unsigned int i = 0; i < m_cellulesEnfants.size(); i++) {
    cellulesLvl[m_lvl + 1].push_back(m_cellulesEnfants[i]);
  }
  for (unsigned int i = 0; i < m_bordsEnfantsInternes.size(); i++) {
    bordsLvl[m_lvl + 1].push_back(m_bordsEnfantsInternes[i]);
  }
}

//***********************************************************************

void Cellule::ecritureGnuplotAMR(std::ofstream &fluxFichier, const int &dim, ObjetGeometrique *objet)
{
  bool ecrit(true);
  int dimension(dim);
  Coord position = m_element->getPosition();
  //Si ecriture issue d'une coupe
  if (objet != 0) {
    ecrit = m_element->traverseObjet(*objet);
    position = objet->projettePoint(position);
    dimension = objet->getType() + 1;
  }
  if (ecrit) {
    if (!m_split) {
      //ATTENTION : l ordre d ecriture des donnees doit rester le meme pour la conformite avec la generation des scripts gnuplot
      fluxFichier << position.getX() << " ";
      if (dimension == 2) fluxFichier << position.getY() << " ";
      if (dimension == 3)fluxFichier << position.getZ() << " ";
      this->ecritPhasesMelange(m_nombrePhases, m_nombreTransports, fluxFichier);
      fluxFichier << m_lvl << " " << m_xi << " ";
      fluxFichier << endl;
    }
    else {
      for (unsigned int i = 0; i < m_cellulesEnfants.size(); i++) {
        m_cellulesEnfants[i]->ecritureGnuplotAMR(fluxFichier, dim, objet);
      }
    }
  } //Fin ecrit
}

//***********************************************************************

int Cellule::getLvl()
{
  return m_lvl;
}

//***********************************************************************

bool Cellule::getSplit()
{
	return m_split;
}

//***********************************************************************

double Cellule::getXi()
{
  return m_xi;
}

//***********************************************************************

void Cellule::setXi(double valeur)
{
  m_xi = valeur;
}

//***********************************************************************

void Cellule::ajoutFluxXi(double valeur)
{
  m_consXi += valeur;
}

//***********************************************************************

void Cellule::retireFluxXi(double valeur)
{
  m_consXi -= valeur;
}

//***********************************************************************

int Cellule::getNombreCellulesEnfants()
{
  return m_cellulesEnfants.size();
}

//***********************************************************************

Cellule* Cellule::getCelluleEnfant(const int &num)
{
  return m_cellulesEnfants[num];
}

//****************************************************************************
//************************** Parallele non-AMR *******************************
//****************************************************************************

void Cellule::rempliTamponPrimitives(double *tampon, int &compteur, Prim type) const
{
	for (int k = 0; k < m_nombrePhases; k++) {
		this->getPhase(k, type)->rempliTampon(tampon, compteur);
	}
	this->getMelange(type)->rempliTampon(tampon, compteur);
	for (int k = 0; k < m_nombreTransports; k++) {
		tampon[++compteur] = this->getTransport(k, type).getValeur();
	}
}

//***********************************************************************

void Cellule::recupereTamponPrimitives(double *tampon, int &compteur, Eos **eos, Prim type)
{
	for (int k = 0; k < m_nombrePhases; k++) {
		this->getPhase(k, type)->recupereTampon(tampon, compteur, eos);
	}
	this->getMelange(type)->recupereTampon(tampon, compteur);
	for (int k = 0; k < m_nombreTransports; k++) {
		this->setTransport(tampon[++compteur], k, type);
	}
	this->calculsEtendusPourCommunications(m_nombrePhases, type);
}

//***********************************************************************

void Cellule::rempliTamponScalaire(double *tampon, int &compteur, string nomVariable) const
{
	tampon[++compteur] = this->selectionneScalaire(nomVariable);
}

//***********************************************************************

void Cellule::recupereTamponScalaire(double *tampon, int &compteur, string nomVariable)
{
	this->setScalaire(nomVariable, tampon[++compteur]);
}

//***********************************************************************

void Cellule::rempliTamponVecteur(double *tampon, int &compteur, const int &dim, std::string nomVecteur, int num, int indice) const
{
	tampon[++compteur] = this->selectionneVecteur(nomVecteur, num, indice).getX();
	if (dim > 1) tampon[++compteur] = this->selectionneVecteur(nomVecteur, num, indice).getY();
	if (dim > 2) tampon[++compteur] = this->selectionneVecteur(nomVecteur, num, indice).getZ();
}

//***********************************************************************

void Cellule::recupereTamponVecteur(double *tampon, int &compteur, const int &dim, std::string nomVecteur, int num, int indice)
{
	Coord temp;
	temp.setX(tampon[++compteur]);
	if (dim > 1) temp.setY(tampon[++compteur]);
	if (dim > 2) temp.setZ(tampon[++compteur]);
	this->setVecteur(nomVecteur, temp, num, indice);
}

//***********************************************************************

void Cellule::rempliTamponTransports(double *tampon, int &compteur) const
{
  for (int k = 0; k < m_nombreTransports; k++) {
    tampon[++compteur] = this->getTransport(k).getValeur();
  }
}

//***********************************************************************

void Cellule::recupereTamponTransports(double *tampon, int &compteur)
{
  for (int k = 0; k < m_nombreTransports; k++) {
    this->setTransport(tampon[++compteur], k);
  }
}

//****************************************************************************
//**************************** AMR Parallele *********************************
//****************************************************************************

void Cellule::rempliTamponPrimitivesAMRjeSuisCpuGauche(double *tampon, int &compteur, const int &lvl, Prim type) const
{
	if (m_lvl == lvl) {
		for (int k = 0; k < m_nombrePhases; k++) {
			this->getPhase(k, type)->rempliTampon(tampon, compteur);
		}
		this->getMelange(type)->rempliTampon(tampon, compteur);
		for (int k = 0; k < m_nombreTransports; k++) {
			tampon[++compteur] = this->getTransport(k, type).getValeur();
		}
	}
	else {
		for (unsigned int i = 0; i < m_cellulesEnfants.size(); i++) {
			//Je suis CPU gauche, j'envoie donc les elements limites de droite
			if ((i % 2) == 1) { m_cellulesEnfants[i]->rempliTamponPrimitivesAMRjeSuisCpuGauche(tampon, compteur, lvl, type); }
		}
	}
}

//***********************************************************************

void Cellule::rempliTamponPrimitivesAMRjeSuisCpuDroite(double *tampon, int &compteur, const int &lvl, Prim type) const
{
	if (m_lvl == lvl) {
		for (int k = 0; k < m_nombrePhases; k++) {
			this->getPhase(k, type)->rempliTampon(tampon, compteur);
		}
		this->getMelange(type)->rempliTampon(tampon, compteur);
		for (int k = 0; k < m_nombreTransports; k++) {
			tampon[++compteur] = this->getTransport(k, type).getValeur();
		}
	}
	else {
		for (unsigned int i = 0; i < m_cellulesEnfants.size(); i++) {
			//Je suis CPU droite, j'envoie donc les elements limites de gauche
			if ((i % 2) == 0) { m_cellulesEnfants[i]->rempliTamponPrimitivesAMRjeSuisCpuDroite(tampon, compteur, lvl, type); }
		}
	}
}

//***********************************************************************

void Cellule::recupereTamponPrimitivesAMR(double *tampon, int &compteur, const int &lvl, Eos **eos, Prim type)
{
	if (m_lvl == lvl) {
		for (int k = 0; k < m_nombrePhases; k++) {
			this->getPhase(k, type)->recupereTampon(tampon, compteur, eos);
		}
		this->getMelange(type)->recupereTampon(tampon, compteur);
		for (int k = 0; k < m_nombreTransports; k++) {
			this->setTransport(tampon[++compteur], k, type);
		}
		this->calculsEtendusPourCommunications(m_nombrePhases, type);
	}
	else {
		for (unsigned int i = 0; i < m_cellulesEnfants.size(); i++) {
			m_cellulesEnfants[i]->recupereTamponPrimitivesAMR(tampon, compteur, lvl, eos, type);
		}
	}
}

//***********************************************************************

void Cellule::rempliTamponScalaireAMRjeSuisCpuGauche(double *tampon, int &compteur, const int &lvl, string nomVariable) const
{
	if (m_lvl == lvl) {
		tampon[++compteur] = this->selectionneScalaire(nomVariable);
	}
	else {
		for (unsigned int i = 0; i < m_cellulesEnfants.size(); i++) {
			//Je suis CPU gauche, j'envoie donc les elements limites de droite
			if ((i % 2) == 1) { m_cellulesEnfants[i]->rempliTamponScalaireAMRjeSuisCpuGauche(tampon, compteur, lvl, nomVariable); }
		}
	}
}

//***********************************************************************

void Cellule::rempliTamponScalaireAMRjeSuisCpuDroite(double *tampon, int &compteur, const int &lvl, string nomVariable) const
{
	if (m_lvl == lvl) {
		tampon[++compteur] = this->selectionneScalaire(nomVariable);
	}
	else {
		for (unsigned int i = 0; i < m_cellulesEnfants.size(); i++) {
			//Je suis CPU droite, j'envoie donc les elements limites de gauche
			if ((i % 2) == 0) { m_cellulesEnfants[i]->rempliTamponScalaireAMRjeSuisCpuDroite(tampon, compteur, lvl, nomVariable); }
		}
	}
}

//***********************************************************************

void Cellule::recupereTamponScalaireAMR(double *tampon, int &compteur, const int &lvl, string nomVariable)
{
	if (m_lvl == lvl) {
		this->setScalaire(nomVariable, tampon[++compteur]);
	}
	else {
		for (unsigned int i = 0; i < m_cellulesEnfants.size(); i++) {
			m_cellulesEnfants[i]->recupereTamponScalaireAMR(tampon, compteur, lvl, nomVariable);
		}
	}
}

//***********************************************************************

void Cellule::rempliTamponVecteurAMRjeSuisCpuGauche(double *tampon, int &compteur, const int &lvl, const int &dim, std::string nomVecteur, int num, int indice) const
{
	if (m_lvl == lvl) {
		tampon[++compteur] = this->selectionneVecteur(nomVecteur, num, indice).getX();
		if (dim > 1) tampon[++compteur] = this->selectionneVecteur(nomVecteur, num, indice).getY();
		if (dim > 2) tampon[++compteur] = this->selectionneVecteur(nomVecteur, num, indice).getZ();
	}
	else {
		for (unsigned int i = 0; i < m_cellulesEnfants.size(); i++) {
			//Je suis CPU gauche, j'envoie donc les elements limites de droite
			if ((i % 2) == 1) { m_cellulesEnfants[i]->rempliTamponVecteurAMRjeSuisCpuGauche(tampon, compteur, lvl, dim, nomVecteur, num, indice); }
		}
	}
}

//***********************************************************************

void Cellule::rempliTamponVecteurAMRjeSuisCpuDroite(double *tampon, int &compteur, const int &lvl, const int &dim, std::string nomVecteur, int num, int indice) const
{
	if (m_lvl == lvl) {
		tampon[++compteur] = this->selectionneVecteur(nomVecteur, num, indice).getX();
		if (dim > 1) tampon[++compteur] = this->selectionneVecteur(nomVecteur, num, indice).getY();
		if (dim > 2) tampon[++compteur] = this->selectionneVecteur(nomVecteur, num, indice).getZ();
	}
	else {
		for (unsigned int i = 0; i < m_cellulesEnfants.size(); i++) {
			//Je suis CPU droite, j'envoie donc les elements limites de gauche
			if ((i % 2) == 0) { m_cellulesEnfants[i]->rempliTamponVecteurAMRjeSuisCpuDroite(tampon, compteur, lvl, dim, nomVecteur, num, indice); }
		}
	}
}

//***********************************************************************

void Cellule::recupereTamponVecteurAMR(double *tampon, int &compteur, const int &lvl, const int &dim, std::string nomVecteur, int num, int indice)
{
	if (m_lvl == lvl) {
		Coord temp;
		temp.setX(tampon[++compteur]);
		if (dim > 1) temp.setY(tampon[++compteur]);
		if (dim > 2) temp.setZ(tampon[++compteur]);
		this->setVecteur(nomVecteur, temp, num, indice);
	}
	else {
		for (unsigned int i = 0; i < m_cellulesEnfants.size(); i++) {
			m_cellulesEnfants[i]->recupereTamponVecteurAMR(tampon, compteur, lvl, dim, nomVecteur, num, indice);
		}
	}
}

//***********************************************************************

void Cellule::rempliTamponTransportsAMRjeSuisCpuGauche(double *tampon, int &compteur, const int &lvl) const
{
  if (m_lvl == lvl) {
    for (int k = 0; k < m_nombreTransports; k++) {
      tampon[++compteur] = this->getTransport(k).getValeur();
    }
  }
  else {
    for (unsigned int i = 0; i < m_cellulesEnfants.size(); i++) {
      //Je suis CPU gauche, j'envoie donc les elements limites de droite
      if ((i % 2) == 1) { m_cellulesEnfants[i]->rempliTamponTransportsAMRjeSuisCpuGauche(tampon, compteur, lvl); }
    }
  }
}

//***********************************************************************

void Cellule::rempliTamponTransportsAMRjeSuisCpuDroite(double *tampon, int &compteur, const int &lvl) const
{
  if (m_lvl == lvl) {
    for (int k = 0; k < m_nombreTransports; k++) {
      tampon[++compteur] = this->getTransport(k).getValeur();
    }
  }
  else {
    for (unsigned int i = 0; i < m_cellulesEnfants.size(); i++) {
      //Je suis CPU droite, j'envoie donc les elements limites de gauche
      if ((i % 2) == 0) { m_cellulesEnfants[i]->rempliTamponTransportsAMRjeSuisCpuDroite(tampon, compteur, lvl); }
    }
  }
}

//***********************************************************************

void Cellule::recupereTamponTransportsAMR(double *tampon, int &compteur, const int &lvl)
{
  if (m_lvl == lvl) {
    for (int k = 0; k < m_nombreTransports; k++) {
      this->setTransport(tampon[++compteur], k);
    }
  }
  else {
    for (unsigned int i = 0; i < m_cellulesEnfants.size(); i++) {
      m_cellulesEnfants[i]->recupereTamponTransportsAMR(tampon, compteur, lvl);
    }
  }
}

//***********************************************************************

void Cellule::choixRaffineDeraffineGhost(const int &nbMaillesY, const int &nbMaillesZ, const double &dX, const double &dY, const double &dZ,
	const vector<PhysAdd*> &physAdd, Modele *modele, vector<Cellule *> *cellulesLvlGhost)
{
  if (m_split) {
    if (m_cellulesEnfants.size() == 0) { this->raffineCelluleEtBordsGhost(nbMaillesY, nbMaillesZ, dX, dY, dZ, physAdd, modele); }
  }
  else {
    if (m_cellulesEnfants.size() > 0) { this->deraffineCelluleEtBordsGhost(); }
  }
  for (unsigned int i = 0; i < m_cellulesEnfants.size(); i++) {
    cellulesLvlGhost[m_lvl + 1].push_back(m_cellulesEnfants[i]);
  }
}

//***********************************************************************

void Cellule::raffineCelluleEtBordsGhost(const int &nbMaillesY, const int &nbMaillesZ, const double &dX, const double &dY, const double &dZ,
	const vector<PhysAdd*> &physAdd, Modele *modele)
{
	//--------------------------------------
	//Initialisations (enfants et dimension)
	//--------------------------------------

	//A noter que le nombre de cellules enfants pour les cellules fantomes est different de celui des cellules internes
	double dimX(1.), dimY(0.), dimZ(0.);
	int nombreCellulesEnfants(1);
	int dim(1);
	if (nbMaillesZ != 1) {
		nombreCellulesEnfants = 4;
		dimY = 1.;
		dimZ = 1.;
		dim = 3;
	}
	else if (nbMaillesY != 1) {
		nombreCellulesEnfants = 2;
		dimY = 1.;
		dim = 2;
	}
	BordDeMaille* bordRef(0);
	for (unsigned int b = 0; b < m_bords.size(); b++) {
		if (m_bords[b]->quiSuisJe() == 0) { bordRef = m_bords[b]; break; } //Bord de type BordDeMaille/O2
	}
	double surfaceEnfant(bordRef->getFace()->getSurface() * pow(2., (double)((dim - 1)*(bordRef->getLvl() - m_lvl - 1))));
	int allouePenteLocal = 1;

	//------------------------
	//Raffinement des cellules
	//------------------------

	//Initialisation des donnees de maillage pour les cellules enfants
	//----------------------------------------------------------------
	double posXCelluleParent, posYCelluleParent, posZCelluleParent;
	double dXParent, dYParent, dZParent;
	double posXEnfant, posYEnfant, posZEnfant;
	posXCelluleParent = m_element->getPosition().getX();
	posYCelluleParent = m_element->getPosition().getY();
	posZCelluleParent = m_element->getPosition().getZ();
	dXParent = dX / pow(2., (double)m_lvl);
	dYParent = dY / pow(2., (double)m_lvl);
	dZParent = dZ / pow(2., (double)m_lvl);
	double volumeCelluleParent, lCFLCelluleParent;
	volumeCelluleParent = m_element->getVolume();
	lCFLCelluleParent = m_element->getLCFL();

	for (int i = 0; i < nombreCellulesEnfants; i++) {

		//Creation des cellules enfants
		//-----------------------------
		this->creerCelluleEnfant(i, m_lvl);
		m_element->creerElementEnfant();
		m_cellulesEnfants[i]->setElement(m_element->getElementEnfant(i), i);
		m_cellulesEnfants[i]->getElement()->setVolume(volumeCelluleParent / (double)nombreCellulesEnfants / 2);
		m_cellulesEnfants[i]->getElement()->setLCFL(lCFLCelluleParent / 2.);
		if (bordRef->getFace()->getPos().getX() < m_element->getPosition().getX()) {
			//La cellule fantome est a droite des cellules internes, donc sa position est colle a gauche
			posXEnfant = posXCelluleParent - 0.25*dimX*dXParent;
		}
		else {
			//La cellule fantome est a gauche des cellules internes, donc sa position est colle a droite
			posXEnfant = posXCelluleParent + 0.25*dimX*dXParent;
		}
		posYEnfant = posYCelluleParent + dimY*dYParent*(double)(-0.25 + 0.5 * (i % 2));
		posZEnfant = posZCelluleParent + dimZ*dZParent*(double)(-0.25 + 0.5 * ((i / 2) % 2));
		m_cellulesEnfants[i]->getElement()->setPos(posXEnfant, posYEnfant, posZEnfant);

		//Initialisation des tableaux principaux selon modele et nombre de phases
		//et initialisation Physique : Remplissage donnees physiques, pour les cellules enfants
		//-------------------------------------------------------------------------------------
		m_cellulesEnfants[i]->alloue(m_nombrePhases, m_nombreTransports, physAdd, modele);
		for (int k = 0; k < m_nombrePhases; k++) {
			m_cellulesEnfants[i]->copiePhase(k, m_vecPhases[k]);
		}
		m_cellulesEnfants[i]->copieMelange(m_melange);
		m_cellulesEnfants[i]->getCons()->miseAZero(m_nombrePhases);
		for (int k = 0; k < m_nombreTransports; k++) { m_cellulesEnfants[i]->setTransport(m_vecTransports[k].getValeur(), k); }
		for (int k = 0; k < m_nombreTransports; k++) { m_cellulesEnfants[i]->setConsTransport(0., k); }
		m_cellulesEnfants[i]->setXi(m_xi);
	}

	//------------------------------
	//Raffinement des bords externes
	//------------------------------

	for (unsigned int b = 0; b < m_bords.size(); b++) {
		if (!m_bords[b]->getSplit()) { m_bords[b]->raffineBordExterneGhost(nbMaillesY, nbMaillesZ, dXParent, dYParent, dZParent, this, surfaceEnfant); }
	}
}

//***********************************************************************

void Cellule::deraffineCelluleEtBordsGhost()
{
	//--------------------------------------
	//Suppression des bords enfants externes
	//--------------------------------------

	for (unsigned int b = 0; b < m_bords.size(); b++) {
		m_bords[b]->deraffineBordExterne(this);
	}

	//--------------------------------
	//Suppression des cellules enfants
	//--------------------------------

	for (unsigned int i = 0; i < m_cellulesEnfants.size(); i++) {
		delete m_cellulesEnfants[i];
	}
	m_cellulesEnfants.clear();
	m_element->finaliseElementsEnfants();
}

//***********************************************************************

void Cellule::rempliTamponXiJeSuisCpuGauche(double *tampon, int &compteur, const int &lvl) const
{
	if (m_lvl == lvl) {
		tampon[++compteur] = m_xi;
	}
	else {
		for (unsigned int i = 0; i < m_cellulesEnfants.size(); i++) {
			//Je suis CPU gauche, j'envoie donc les elements limites de droite
			if ((i % 2) == 1) { m_cellulesEnfants[i]->rempliTamponXiJeSuisCpuGauche(tampon, compteur, lvl); }
		}
	}
}

//***********************************************************************

void Cellule::rempliTamponXiJeSuisCpuDroite(double *tampon, int &compteur, const int &lvl) const
{
	if (m_lvl == lvl) {
		tampon[++compteur] = m_xi;
	}
	else {
		for (unsigned int i = 0; i < m_cellulesEnfants.size(); i++) {
			//Je suis CPU droite, j'envoie donc les elements limites de gauche
			if ((i % 2) == 0) { m_cellulesEnfants[i]->rempliTamponXiJeSuisCpuDroite(tampon, compteur, lvl); }
		}
	}
}

//***********************************************************************

void Cellule::recupereTamponXi(double *tampon, int &compteur, const int &lvl)
{
	if (m_lvl == lvl) {
		m_xi = tampon[++compteur];
	}
	else {
		for (unsigned int i = 0; i < m_cellulesEnfants.size(); i++) {
			m_cellulesEnfants[i]->recupereTamponXi(tampon, compteur, lvl);
		}
	}
}

//***********************************************************************

void Cellule::rempliTamponSplitJeSuisCpuGauche(bool *tampon, int &compteur, const int &lvl) const
{
	if (m_lvl == lvl) {
		tampon[++compteur] = m_split;
	}
	else {
		for (unsigned int i = 0; i < m_cellulesEnfants.size(); i++) {
			//Je suis CPU gauche, j'envoie donc les elements limites de droite
			if ((i % 2) == 1) { m_cellulesEnfants[i]->rempliTamponSplitJeSuisCpuGauche(tampon, compteur, lvl); }
		}
	}
}

//***********************************************************************

void Cellule::rempliTamponSplitJeSuisCpuDroite(bool *tampon, int &compteur, const int &lvl) const
{
	if (m_lvl == lvl) {
		tampon[++compteur] = m_split;
	}
	else {
		for (unsigned int i = 0; i < m_cellulesEnfants.size(); i++) {
			//Je suis CPU droite, j'envoie donc les elements limites de gauche
			if ((i % 2) == 0) { m_cellulesEnfants[i]->rempliTamponSplitJeSuisCpuDroite(tampon, compteur, lvl); }
		}
	}
}

//***********************************************************************

void Cellule::recupereTamponSplit(bool *tampon, int &compteur, const int &lvl)
{
	if (m_lvl == lvl) {
		m_split = tampon[++compteur];
	}
	else {
		for (unsigned int i = 0; i < m_cellulesEnfants.size(); i++) {
			m_cellulesEnfants[i]->recupereTamponSplit(tampon, compteur, lvl);
		}
	}
}

//***********************************************************************

void Cellule::rempliNombreElementsAEnvoyerAVoisinJeSuisCpuGauche(int &NombreElementsAEnvoyerAVoisin, const int &lvl)
{
	if (m_lvl == lvl) {
		NombreElementsAEnvoyerAVoisin++;
	}
	else {
		for (unsigned int i = 0; i < m_cellulesEnfants.size(); i++) {
			//Je suis CPU gauche, j'envoie donc les elements limites de droite
			if ((i % 2) == 1) { m_cellulesEnfants[i]->rempliNombreElementsAEnvoyerAVoisinJeSuisCpuGauche(NombreElementsAEnvoyerAVoisin, lvl); }
		}
	}
}

//***********************************************************************

void Cellule::rempliNombreElementsAEnvoyerAVoisinJeSuisCpuDroite(int &NombreElementsAEnvoyerAVoisin, const int &lvl)
{
	if (m_lvl == lvl) {
		NombreElementsAEnvoyerAVoisin++;
	}
	else {
		for (unsigned int i = 0; i < m_cellulesEnfants.size(); i++) {
			//Je suis CPU droite, j'envoie donc les elements limites de gauche
			if ((i % 2) == 0) { m_cellulesEnfants[i]->rempliNombreElementsAEnvoyerAVoisinJeSuisCpuDroite(NombreElementsAEnvoyerAVoisin, lvl); }
		}
	}
}

//***************************************************************************