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

#include "Entrees.h"
#include "EnteteEntreesSorties.h"
#include <iostream>

using namespace std;
using namespace tinyxml2;

//***********************************************************************

Entrees::Entrees(Run *run) : m_run(run)
{
	//Attribution des numeros de version pour la lecture
	m_vMain = 5;
	m_vMaillage = 5;
	m_vCI = 4;
	m_vModele = 4;

	m_nomMain = "mainV" + IO::toString(m_vMain) + ".xml";
	m_nomMaillage = "maillageV" + IO::toString(m_vMaillage) + ".xml";
	m_nomCI = "conditionsInitialesV" + IO::toString(m_vCI) + ".xml";
	m_nomModele = "modeleV" + IO::toString(m_vModele) + ".xml";
}

//***********************************************************************

Entrees::~Entrees(){}

//***********************************************************************

void Entrees::lectureEntreesXML(vector<DomaineGeometrique*> &domaines, vector<CondLim*> &condLim)
{
  try{
    //1) Parametres generaux du calcul
    entreeMain(m_run->m_casTest);
    //2) Donnees de Maillage
    entreeMaillage(m_run->m_casTest);
    //3) Donnees modeles et fluides
    entreeModele(m_run->m_casTest);
    //4) Lecture des conditions initiales
    entreeConditionsInitiales(m_run->m_casTest, domaines, condLim);
  }
  catch (ErreurXML &e){
    cerr << e.infoErreur() << endl;
    throw;
  }
}

//***********************************************************************

void Entrees::entreeMain(string casTest)
{
  try{
    //1) Parsing du fichier XML par la bibliotheque tinyxml2
    //------------------------------------------------------
	stringstream nomFichier(casTest + m_nomMain);
    XMLDocument xmlMain;
    XMLError erreur(xmlMain.LoadFile(nomFichier.str().c_str())); //Le fichier est parse ici
    if (erreur != XML_SUCCESS) throw ErreurXML(nomFichier.str(),__FILE__, __LINE__);
    
    //2) Recuperation des donnees principales du calcul
    //-------------------------------------------------
    //Recuperation racine du document XML
    XMLNode *paramCalculs = xmlMain.FirstChildElement("paramCalculs");
    if (paramCalculs == NULL) throw ErreurXMLRacine("paramCalculs", nomFichier.str(), __FILE__, __LINE__);

    XMLElement *element, *sousElement;

    //Recuperation nom du run
    element = paramCalculs->FirstChildElement("run");
    if (element == NULL) throw ErreurXMLElement("run", nomFichier.str(), __FILE__, __LINE__);
    XMLNode* xmlNode2 = element->FirstChild();
    if (xmlNode2 == NULL) throw ErreurXMLElement("run", nomFichier.str(), __FILE__, __LINE__);
    XMLText* xmlText = xmlNode2->ToText();
    if (xmlText == NULL) throw ErreurXMLElement("run", nomFichier.str(), __FILE__, __LINE__);

    //Lecture des informations de sorties/ecritures
    element = paramCalculs->FirstChildElement("modeSortie");
    if (element == NULL) throw ErreurXMLElement("modeSortie", nomFichier.str(), __FILE__, __LINE__);
    //Lecture format sortie
    string format(element->Attribute("format"));
    if (format == "") throw ErreurXMLAttribut("format", nomFichier.str(), __FILE__, __LINE__);
    Outils::majuscule(format);
    if (format == "XML") { m_run->m_sortie = new SortiesXML(casTest, xmlText->Value(), element, nomFichier.str(), this); }
    else if (format == "GNU") { m_run->m_sortie = new SortiesGNU(casTest, xmlText->Value(), element, nomFichier.str(), this); }
    else { throw ErreurXMLDev(nomFichier.str(), __FILE__, __LINE__); }

    //Lecture des coupes 1D
    element = paramCalculs->FirstChildElement("coupe1D");
    while (element != NULL)
    {
      m_run->m_coupes.push_back(new SortiesCoupeGNU(casTest, xmlText->Value(), element, nomFichier.str(), DROITE, this));
      element = element->NextSiblingElement("coupe1D");
    }
    //Lecture des coupes 2D
    element = paramCalculs->FirstChildElement("coupe2D");
    while (element != NULL)
    {
      m_run->m_coupes.push_back(new SortiesCoupeGNU(casTest, xmlText->Value(), element, nomFichier.str(), PLAN, this));
      element = element->NextSiblingElement("coupe2D");
    }
    
    //Recuperation Iteration / temps Physique
    element = paramCalculs->FirstChildElement("modeControleTemporel");
    if (element == NULL) throw ErreurXMLElement("modeControleTemporel", nomFichier.str(), __FILE__, __LINE__);
    erreur = element->QueryBoolAttribute("iterations", &m_run->m_controleIterations);
    if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut("iterations", nomFichier.str(), __FILE__, __LINE__);
    if (m_run->m_controleIterations)
    {
      //Recuperation Iterations / Frequence
      sousElement = element->FirstChildElement("iterations");
      if (sousElement == NULL) throw ErreurXMLElement("iterations", nomFichier.str(), __FILE__, __LINE__);
      erreur = sousElement->QueryIntAttribute("nombre", &m_run->m_nbIte);
      if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut("nombre", nomFichier.str(), __FILE__, __LINE__);
      erreur = sousElement->QueryIntAttribute("freqImpr", &m_run->m_freq);
      if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut("freqImpr", nomFichier.str(), __FILE__, __LINE__);
    }
    else
    {
      //Recuperation Temps / Frequence
      sousElement = element->FirstChildElement("temporel");
      if (sousElement == NULL) throw ErreurXMLElement("temporel", nomFichier.str(), __FILE__, __LINE__);
      erreur = sousElement->QueryFloatAttribute("tempsTotal", &m_run->m_tempsFinal);
      if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut("tempsTotal", nomFichier.str(), __FILE__, __LINE__);
      erreur = sousElement->QueryFloatAttribute("freqTemp", &m_run->m_freqTemps);
      if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut("freqTemp", nomFichier.str(), __FILE__, __LINE__);
    }

    //Recuperation CFL
    element = paramCalculs->FirstChildElement("controleCalcul");
    if (element == NULL) throw ErreurXMLElement("controleCalcul", nomFichier.str(), __FILE__, __LINE__);
    erreur = element->QueryDoubleAttribute("CFL", &m_run->m_cfl);
    if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut("CFL", nomFichier.str(), __FILE__, __LINE__);

    //Lecture Ordre2
    element = paramCalculs->FirstChildElement("ordre2");
    if(element == NULL) {
      m_run->m_ordre = "ORDRE1";
      m_run->m_limiteurGlobal = new Limiteur;
    }
    else{
      m_run->m_ordre = "ORDRE2";
      XMLNode *contenu;
      //Recuperation limiteur Global
      sousElement = element->FirstChildElement("limiteurGlobal");
      if (sousElement == NULL) throw ErreurXMLElement("limiteurGlobal", nomFichier.str(), __FILE__, __LINE__);
      contenu = sousElement->FirstChild();
      if (contenu == NULL) throw ErreurXMLElement("limiteurGlobal", nomFichier.str(), __FILE__, __LINE__);
      string limiteurGlobal = contenu->ToText()->Value();
      Outils::majuscule(limiteurGlobal);
      if (limiteurGlobal == "MINMOD") { m_run->m_limiteurGlobal = new LimiteurMinmod; }
      else if (limiteurGlobal == "VANLEER") { m_run->m_limiteurGlobal = new LimiteurVanLeer; }
      else if (limiteurGlobal == "VANALBADA") { m_run->m_limiteurGlobal = new LimiteurVanAlbada; }
      else if (limiteurGlobal == "SUPERBEE") { m_run->m_limiteurGlobal = new LimiteurSuperBee; }
      else if (limiteurGlobal == "MC") { m_run->m_limiteurGlobal = new LimiteurMC; }
      else { throw ErreurXMLDev(nomFichier.str(), __FILE__, __LINE__); }
      //Recuperation limiteur Interface
      string limiteurInterface = limiteurGlobal;
      sousElement = element->FirstChildElement("limiteurInterface");
      if (sousElement != NULL) {
        contenu = sousElement->FirstChild();
        if (contenu == NULL) throw ErreurXMLElement("limiteurInterface", nomFichier.str(), __FILE__, __LINE__);
        limiteurInterface = contenu->ToText()->Value();
      }
      Outils::majuscule(limiteurInterface);
      if (limiteurInterface == "MINMOD") { m_run->m_limiteurInterface = new LimiteurMinmod; }
      else if (limiteurInterface == "VANLEER") { m_run->m_limiteurInterface = new LimiteurVanLeer; }
      else if (limiteurInterface == "VANALBADA") { m_run->m_limiteurInterface = new LimiteurVanAlbada; }
      else if (limiteurInterface == "SUPERBEE") { m_run->m_limiteurInterface = new LimiteurSuperBee; }
      else if (limiteurInterface == "MC") { m_run->m_limiteurInterface = new LimiteurMC; }
      else { throw ErreurXMLDev(nomFichier.str(), __FILE__, __LINE__); }

    }

    //Reprise de Calcul depuis fichier resultat
    element = paramCalculs->FirstChildElement("repriseFichier");
    if (element != NULL) { 
      erreur = element->QueryIntAttribute("numero", &m_run->m_repriseFichier);
      if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut("numero", nomFichier.str(), __FILE__, __LINE__);
    }

  }
  catch (ErreurXML &){ throw; } // Renvoi au niveau suivant
}

//***********************************************************************

void Entrees::entreeMaillage(string casTest)
{
  try{
    //Methode AMR: initialisation des variables
    m_run->m_lvlMax = 0;
		double critereVar(1.e10);
		bool varRho(false), varP(false), varU(false), varAlpha(false);
		double xiSplit(1.), xiJoin(1.);

    //1) Parsing du fichier XML par la bibliotheque tinyxml2
    //------------------------------------------------------
    stringstream nomFichier(casTest + m_nomMaillage);
    XMLDocument xmlMaillage;
    XMLError erreur(xmlMaillage.LoadFile(nomFichier.str().c_str())); //Le fichier est parse ici
    if (erreur != XML_SUCCESS) throw ErreurXML(nomFichier.str(), __FILE__, __LINE__);

    //2) Recuperation des donnees principales du calcul
    //-------------------------------------------------
    //Recuperation racine du document XML
    XMLNode *maillage = xmlMaillage.FirstChildElement("maillage");
    if (maillage == NULL) throw ErreurXMLRacine("maillage", nomFichier.str(), __FILE__, __LINE__);

    XMLElement *element;
    //Recuperation du type de maillage
    element = maillage->FirstChildElement("type");
    if (element == NULL) throw ErreurXMLElement("type", nomFichier.str(), __FILE__, __LINE__);
    string structureMaillage(element->Attribute("structure"));
    if (structureMaillage == "") throw ErreurXMLAttribut("structure", nomFichier.str(), __FILE__, __LINE__);
    Outils::majuscule(structureMaillage);
    if (structureMaillage == "NONSTRUCTURE")
    {
      //----------------MAILLAGE NON STRUCTURE ---------------------
      XMLElement *maillageNS;
      maillageNS = maillage->FirstChildElement("maillageNonStructure");
      if (maillageNS == NULL) throw ErreurXMLElement("maillageNonStructure", nomFichier.str(), __FILE__, __LINE__);
      //Lecture nom du fichier contenant les informations de maillage (.msh)
      element = maillageNS->FirstChildElement("fichier");
      if (element == NULL) throw ErreurXMLElement("fichier", nomFichier.str(), __FILE__, __LINE__);
      string fichierMaillage(element->Attribute("nom"));
      if (fichierMaillage == "") throw ErreurXMLAttribut("nom", nomFichier.str(), __FILE__, __LINE__);
      m_run->m_maillage = new MaillageNonStruct(fichierMaillage);
      //Recuperation pretraitement parallele
      element = maillageNS->FirstChildElement("modeParallele");
      if (element != NULL) {
        erreur = element->QueryBoolAttribute("pretraitementGMSH", &m_run->m_pretraitementParallele);
        if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut("pretraitementGMSH", nomFichier.str(), __FILE__, __LINE__);
      }
      //Methode AMR non possible avec maillage non structure
      element = maillageNS->FirstChildElement("AMR");
      if (element != NULL) { throw ErreurXMLAttribut("Methode AMR non possible avec maillage non structure", nomFichier.str(), __FILE__, __LINE__); }
    }
    else if (structureMaillage == "CARTESIEN")
    {
      //----------------MAILLAGE CARTESIEN ---------------------
      XMLElement *maillageCartesien;
      maillageCartesien = maillage->FirstChildElement("maillageCartesien");
      if (maillageCartesien == NULL) throw ErreurXMLElement("maillageCartesien", nomFichier.str(), __FILE__, __LINE__);
      //Recuperation des dimensions
      double lX, lY, lZ;
      element = maillageCartesien->FirstChildElement("dimensions");
      if (element == NULL) throw ErreurXMLElement("dimensions", nomFichier.str(), __FILE__, __LINE__);
      erreur = element->QueryDoubleAttribute("x", &lX);
      if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut("x", nomFichier.str(), __FILE__, __LINE__);
      erreur = element->QueryDoubleAttribute("y", &lY);
      if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut("y", nomFichier.str(), __FILE__, __LINE__);
      erreur = element->QueryDoubleAttribute("z", &lZ);
      if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut("z", nomFichier.str(), __FILE__, __LINE__);
      //Recuperation des nombres de mailles
      int nbX, nbY, nbZ;
      element = maillageCartesien->FirstChildElement("nombreMailles");
      if (element == NULL) throw ErreurXMLElement("nombreMailles", nomFichier.str(), __FILE__, __LINE__);
      erreur = element->QueryIntAttribute("x", &nbX);
      if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut("x", nomFichier.str(), __FILE__, __LINE__);
      erreur = element->QueryIntAttribute("y", &nbY);
      if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut("y", nomFichier.str(), __FILE__, __LINE__);
      erreur = element->QueryIntAttribute("z", &nbZ);
      if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut("z", nomFichier.str(), __FILE__, __LINE__);
      //Recuperation des variables pour methode AMR
      element = maillageCartesien->FirstChildElement("AMR");
      if (element != NULL) {
        //Interdiction parallele
        erreur = element->QueryIntAttribute("lvlMax", &m_run->m_lvlMax);
        if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut("lvlMax", nomFichier.str(), __FILE__, __LINE__);
        if (m_run->m_repriseFichier && m_run->m_lvlMax > 0) throw ErreurXMLAttribut("Reprise fichier non prevue pour maillage considere (AMR)", nomFichier.str(), __FILE__, __LINE__);
        erreur = element->QueryDoubleAttribute("critereVar", &critereVar);
        if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut("critereVar", nomFichier.str(), __FILE__, __LINE__);
        erreur = element->QueryBoolAttribute("varRho", &varRho);
        if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut("varRho", nomFichier.str(), __FILE__, __LINE__);
        erreur = element->QueryBoolAttribute("varP", &varP);
        if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut("varP", nomFichier.str(), __FILE__, __LINE__);
        erreur = element->QueryBoolAttribute("varU", &varU);
        if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut("varU", nomFichier.str(), __FILE__, __LINE__);
        erreur = element->QueryBoolAttribute("varAlpha", &varAlpha);
        if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut("varAlpha", nomFichier.str(), __FILE__, __LINE__);
        erreur = element->QueryDoubleAttribute("xiSplit", &xiSplit);
        if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut("xiSplit", nomFichier.str(), __FILE__, __LINE__);
        erreur = element->QueryDoubleAttribute("xiJoin", &xiJoin);
        if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut("xiJoin", nomFichier.str(), __FILE__, __LINE__);
        m_run->m_maillage = new MaillageCartesienAMR(lX, nbX, lY, nbY, lZ, nbZ, m_run->m_lvlMax, critereVar, varRho, varP, varU, varAlpha, xiSplit, xiJoin);
				if (Ncpu > 1) {
					if (nbX >= nbY && nbX >= nbZ) {} //Decoupage selon X
					else { throw ErreurXML(nomFichier.str(), __FILE__, __LINE__); } // Erreurs::messageErreur("decoupage en Y ou Z non gere pour maillage AMR"); }
				}
      }
      else {
        m_run->m_maillage = new MaillageCartesien(lX, nbX, lY, nbY, lZ, nbZ);
      }

    }
    else{ throw ErreurXMLDev(nomFichier.str(), __FILE__, __LINE__); }
  }
  catch (ErreurXML &){ throw; } // Renvoi au niveau suivant
}

//***********************************************************************

void Entrees::entreeModele(string casTest)
{
  try{
    //1) Parsing du fichier XML par la bibliotheque tinyxml2
    //------------------------------------------------------
    stringstream nomFichier(casTest + m_nomModele);
    XMLDocument xmlModele;
    XMLError erreur(xmlModele.LoadFile(nomFichier.str().c_str())); //Le fichier est parse ici
    if (erreur != XML_SUCCESS) throw ErreurXML(nomFichier.str(), __FILE__, __LINE__);

    //2) Recuperation des donnees du Modele
    //-------------------------------------
    //Recuperation racine du document XML
    XMLNode *xmlNode = xmlModele.FirstChildElement("modele");
    if (xmlNode == NULL) throw ErreurXMLRacine("modele", nomFichier.str(), __FILE__, __LINE__);
    //Recuperation nom du modele resolu
    XMLElement *element(xmlNode->FirstChildElement("modeleResolu"));
    if (element == NULL) throw ErreurXMLElement("modeleResolu", nomFichier.str(), __FILE__, __LINE__);
    string modele(element->Attribute("nom"));
    if (modele == "") throw ErreurXMLAttribut("nom", nomFichier.str(), __FILE__, __LINE__);
    m_run->m_nombreTransports = 0;
    erreur = element->QueryIntAttribute("nombreTransports", &m_run->m_nombreTransports);
    //if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut("nombreTransports", nomFichier.str(), __FILE__, __LINE__);
    Outils::majuscule(modele);
    //Switch selon modele
    if (modele == "EULER"){ m_run->m_modele = new ModEuler(m_run->m_nombreTransports); m_run->m_nombrePhases = 1; }
    else if (modele == "EULERHOMOGENE") { 
      int liquide, vapeur;
      erreur = element->QueryIntAttribute("liquide", &liquide);
      erreur = element->QueryIntAttribute("vapeur", &vapeur);
      m_run->m_modele = new ModEulerHomogene(m_run->m_nombreTransports, liquide, vapeur); m_run->m_nombrePhases = 2; 
    }
    else if (modele == "KAPILA")
    { 
      erreur = element->QueryIntAttribute("nombrePhases", &m_run->m_nombrePhases);
      m_run->m_modele = new ModKapila(m_run->m_nombreTransports, m_run->m_nombrePhases);
      if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut("nombrePhases", nomFichier.str(), __FILE__, __LINE__);
    }
    else { throw ErreurXMLDev(nomFichier.str(), __FILE__, __LINE__); }
    
    //Thermodynamique
    vector<string> nomEOS;
    int EOSTrouvee(0);
    XMLElement *sousElement(xmlNode->FirstChildElement("equationEtat"));
    while (sousElement != NULL)
    {
      //Lecture Eos
      nomEOS.push_back(sousElement->Attribute("nom"));
      if (nomEOS[EOSTrouvee] == "") throw ErreurXMLAttribut("nom", nomFichier.str(), __FILE__, __LINE__);
      EOSTrouvee++;
      sousElement = sousElement->NextSiblingElement("equationEtat");
    }
    m_run->m_nombreEos = nomEOS.size();
    if (m_run->m_nombreEos == 0) throw ErreurXMLEOS(nomFichier.str(), __FILE__, __LINE__);
    m_run->m_eos = new Eos*[m_run->m_nombrePhases];
    int numeroEos(0);
    for (int i = 0; i < m_run->m_nombreEos; i++){ m_run->m_eos[i] = entreeEOS(nomEOS[i], numeroEos); }

    //Transports
    int transportsTrouvee(0);
    XMLElement *elementTransport(xmlNode->FirstChildElement("transport"));
    while (elementTransport != NULL)
    {
      //Lecture nom transport
      m_run->m_nomGTR.push_back(elementTransport->Attribute("nom"));
      if (m_run->m_nomGTR[transportsTrouvee] == "") throw ErreurXMLAttribut("nom", nomFichier.str(), __FILE__, __LINE__);
      transportsTrouvee++;
      elementTransport = elementTransport->NextSiblingElement("transport");
    }
    //Verifications
    if(transportsTrouvee < m_run->m_nombreTransports) throw ErreurXML(nomFichier.str(), __FILE__, __LINE__);

    //Lecture des Physiques Additionnelles
    int nombreGPA(0);
    int physiqueAddTrouvee(0);
    element = xmlNode->FirstChildElement("physiqueAdd");
    while (element != NULL) {
      physiqueAddTrouvee++;
      //Lecture physique additionnelle
      string typePhysAdd(element->Attribute("type"));
      if (typePhysAdd == "") throw ErreurXMLAttribut("type", nomFichier.str(), __FILE__, __LINE__);
      Outils::majuscule(typePhysAdd);
      //switch sur le type de physique additionelle
      if (typePhysAdd == "CAPILLARITE") { 
        //switch modele d ecoulement
        if (modele == "KAPILA") { m_run->m_physAdd.push_back(new PAKCapillarite(element, nombreGPA, m_run->m_nomGTR, nomFichier.str())); }
        else { throw ErreurXMLDev(nomFichier.str(), __FILE__, __LINE__); }
      }
      else if (typePhysAdd == "VISCOSITE") {
        //switch modele d ecoulement
        if (modele == "KAPILA") { m_run->m_physAdd.push_back(new PAKViscosite(nombreGPA, m_run->m_eos, m_run->m_nombrePhases, nomFichier.str())); }
        else { throw ErreurXMLDev(nomFichier.str(), __FILE__, __LINE__); }
      }
      else if (typePhysAdd == "CONDUCTIVITE") {
        //switch modele d ecoulement
        if (modele == "KAPILA") { m_run->m_physAdd.push_back(new PAKConductivite(nombreGPA, m_run->m_eos, m_run->m_nombrePhases, nomFichier.str())); }
        else { throw ErreurXMLDev(nomFichier.str(), __FILE__, __LINE__); }
      }
      else { throw ErreurXMLDev(nomFichier.str(), __FILE__, __LINE__); }
      element = element->NextSiblingElement("physiqueAdd");
    }
    m_run->m_nombrePhysAdd = physiqueAddTrouvee;

    //Lecture des termes Sources
    int sourceTrouve(0);
    element = xmlNode->FirstChildElement("termeSource");
    while (element != NULL) {
      sourceTrouve++;
      //Lecture source
      string typeSource(element->Attribute("type"));
      if (typeSource == "") throw ErreurXMLAttribut("type", nomFichier.str(), __FILE__, __LINE__);
      Outils::majuscule(typeSource);
      //switch sur le type de source
      if (typeSource == "GRAVITE") { m_run->m_sources.push_back(new SourceGravite(element, nomFichier.str())); }
      else if (typeSource == "AXI") { m_run->m_sources.push_back(new SourceAxi(element, nomFichier.str())); }
      else { throw ErreurXMLDev(nomFichier.str(), __FILE__, __LINE__); }
      element = element->NextSiblingElement("termeSource");
    }
    m_run->m_nombreSources = sourceTrouve;

    //Lecture des relaxations
    element = xmlNode->FirstChildElement("relaxation");
    while (element != NULL) {
      //Lecture source
      string typeRelax(element->Attribute("type"));
      if (typeRelax == "") throw ErreurXMLAttribut("type", nomFichier.str(), __FILE__, __LINE__);
      Outils::majuscule(typeRelax);
      //switch sur le type de source
      if (typeRelax == "PTMU") { m_run->m_evaporation = 1; }
      else { throw ErreurXMLDev(nomFichier.str(), __FILE__, __LINE__); }
      element = element->NextSiblingElement("relaxation");
    }

  }
  catch (ErreurXML &){ throw; } // Renvoi au niveau suivant
}

//***********************************************************************

Eos* Entrees::entreeEOS(string EOS, int &numeroEOS)
{
  try{
    //1) Parsing du fichier XML par la bibliotheque tinyxml2
    //------------------------------------------------------
    stringstream nomFichier("./libMateriaux/" + EOS);
    XMLDocument xmlEOS;
    XMLError erreur(xmlEOS.LoadFile(nomFichier.str().c_str())); //Le fichier est parse ici
    if (erreur != XML_SUCCESS) throw ErreurXML(nomFichier.str(), __FILE__, __LINE__);

    //2) Recuperation des donnees de l'EOS
    //------------------------------------
    //Recuperation racine du document XML
    XMLNode *xmlNode = xmlEOS.FirstChildElement("parametresEOS");
    if (xmlNode == NULL) throw ErreurXMLRacine("parametresEOS", nomFichier.str(), __FILE__, __LINE__);
    //Recuperation type d'EOS
    XMLElement *element;
    element = xmlNode->FirstChildElement("EOS");
    if (element == NULL) throw ErreurXMLElement("EOS", nomFichier.str(), __FILE__, __LINE__);
    string typeEOS(element->Attribute("type"));
    if (typeEOS == "") throw ErreurXMLAttribut("type", nomFichier.str(), __FILE__, __LINE__);
    Outils::majuscule(typeEOS);
    //Switch selon EOS
    Eos *eos;
    vector<string> NomsParametresEos;
    if      (typeEOS == "GP"){ eos = new EosGP(NomsParametresEos, numeroEOS); }
    else if (typeEOS == "SG"){ eos = new EosSG(NomsParametresEos, numeroEOS); }
    else{ throw ErreurXMLEOSInconnue(typeEOS, nomFichier.str(), __FILE__, __LINE__); } //Cas ou la loi etat est inconnue

    //Recuperation des parametres de l'EOS
    element = xmlNode->FirstChildElement("parametres");
    if (element == NULL) throw ErreurXMLElement("parametres", nomFichier.str(), __FILE__, __LINE__);
    vector<double> parametresEos(NomsParametresEos.size());
    for (unsigned int p = 0; p < NomsParametresEos.size(); p++)
    {
      erreur = element->QueryDoubleAttribute(NomsParametresEos[p].c_str(), &parametresEos[p]);
      if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut(NomsParametresEos[p].c_str(), nomFichier.str(), __FILE__, __LINE__);
    }
    eos->attributParametresEos(EOS.c_str(), parametresEos);
    //Lecture des parametres physiques (viscosite, conductivite, etc.)
    eos->lectureParamPhysiques(xmlNode, nomFichier.str());

    return eos;
  }
  catch (ErreurXML &){ throw; } // Renvoi au niveau suivant
}

//***********************************************************************

void Entrees::entreeConditionsInitiales(string casTest, vector<DomaineGeometrique*> &domaines, vector<CondLim*> &condLim)
{
  try{
    //1) Parsing du fichier XML par la bibliotheque tinyxml2
    //------------------------------------------------------
    stringstream nomFichier(casTest + m_nomCI);
    XMLDocument xmlModele;
    XMLError erreur(xmlModele.LoadFile(nomFichier.str().c_str())); //Le fichier est parse ici
    if (erreur != XML_SUCCESS) throw ErreurXML(nomFichier.str(), __FILE__, __LINE__);

    //------------------------------ CONDITIONS INITIALES -------------------------------

    //2) Recuperation racine CI du document XML
    //-----------------------------------------
    XMLNode *xmlNode = xmlModele.FirstChildElement("CI");
    if (xmlNode == NULL) throw ErreurXMLRacine("CI", nomFichier.str(), __FILE__, __LINE__);
    XMLElement *elementDomaine(xmlNode->FirstChildElement("domainesPhysiques"));
    if(elementDomaine == NULL) throw ErreurXMLElement("domainesPhysiques", nomFichier.str(), __FILE__, __LINE__);

    //3) Recuperation des domaines definis
    //------------------------------------
    string nomDomaine, etatDomaine;
    int domaineTrouve(0);
    XMLElement *element(elementDomaine->FirstChildElement("domaine"));
    while (element != NULL)
    {
      //A)Lecture nom domaine
      //*********************
      nomDomaine = element->Attribute("nom");
      if (nomDomaine == "") throw ErreurXMLAttribut("nom", nomFichier.str(), __FILE__, __LINE__);
      Outils::majuscule(nomDomaine);

      //B)Lecture etat domaine
      //**********************
      etatDomaine = element->Attribute("etat");
      if (etatDomaine == "") throw ErreurXMLAttribut("etat", nomFichier.str(), __FILE__, __LINE__);
      Outils::majuscule(etatDomaine);
      //Recherche de l etat associe au domaine
      XMLElement *etat(xmlNode->FirstChildElement("etat"));
      bool trouve(false); string nomEtat;
      while (etat != NULL)
      {
        nomEtat = etat->Attribute("nom");
        Outils::majuscule(nomEtat);
        if (nomEtat == etatDomaine){ trouve = true; break; }
        etat = etat->NextSiblingElement("etat");
      }
      if (!trouve){ throw ErreurXMLEtat(etatDomaine, nomFichier.str(), __FILE__, __LINE__); }
      //On peut lire chaque materiau de l etat trouve
      vector<Phase*> etatsPhases;
      Melange *etatMelange(0);
      string typeMateriau; string nomEOS;
      XMLElement* materiau(etat->FirstChildElement("materiau"));
      int nbMateriauxEtat(0);
      while (materiau != NULL)
      {
        typeMateriau = materiau->Attribute("type");
        Outils::majuscule(typeMateriau);

        //LECTURE FLUIDE
        if (typeMateriau == "FLUIDE") { 
          nbMateriauxEtat++;
          //Recuperation de l EOS
          nomEOS = materiau->Attribute("EOS");
          int e(0);//  bool eosTrouvee(false);
          for (e = 0; e < m_run->m_nombrePhases; e++) {
            if (nomEOS == m_run->m_eos[e]->retourneNom()) { break; }
          }
          if (e == m_run->m_nombrePhases) { throw ErreurXMLEOSInconnue(nomEOS, nomFichier.str(), __FILE__, __LINE__); }
          if (m_run->m_modele->quiSuisJe() == "EULER") { etatsPhases.push_back(new PhaseEuler(materiau, m_run->m_eos[e], nomFichier.str())); }
          else if (m_run->m_modele->quiSuisJe() == "KAPILA") { etatsPhases.push_back(new PhaseKapila(materiau, m_run->m_eos[e], nomFichier.str())); }
          else { throw ErreurXMLElement("Couplage Modele-Fluide non valide", nomFichier.str(), __FILE__, __LINE__); }
        }

        //LECTURE MELANGE BINAIRE LIQUIDE VAPEUR
        else if (typeMateriau == "MELANGEBINAIRELIQVAP") {
          nbMateriauxEtat += 2;
          //FP//TODO// A simplifier car le liquide et la vapeur sont deja connue par le modele
          string nomEOSLiq = materiau->Attribute("EOSLiq");
          int eLiq(0);//  bool eosTrouvee(false);
          for (eLiq = 0; eLiq < m_run->m_nombrePhases; eLiq++) {
            if (nomEOSLiq == m_run->m_eos[eLiq]->retourneNom()) { break; }
          }
          if (eLiq == m_run->m_nombrePhases) { throw ErreurXMLEOSInconnue(nomEOSLiq, nomFichier.str(), __FILE__, __LINE__); }
          string nomEOSVap = materiau->Attribute("EOSVap");
          int eVap(0);//  bool eosTrouvee(false);
          for (eVap = 0; eVap < m_run->m_nombrePhases; eVap++) {
            if (nomEOSVap == m_run->m_eos[eVap]->retourneNom()) { break; }
          }
          if (eVap == m_run->m_nombrePhases) { throw ErreurXMLEOSInconnue(nomEOSVap, nomFichier.str(), __FILE__, __LINE__); }

          XMLElement *donneesMelange(materiau->FirstChildElement("donneesMelangeBinaireLiqVap"));
          double alphaLiq, pression;
          erreur = donneesMelange->QueryDoubleAttribute("alphaLiq", &alphaLiq);
          erreur = donneesMelange->QueryDoubleAttribute("pression", &pression);
 
          //Calcul des donnees fluides a faire a partir de alpha et p ou alpha et T ou alpha et rho, etc.
          MelEulerHomogene melange;
          double Tsat = melange.calculTsat(m_run->m_eos[eLiq], m_run->m_eos[eVap], pression);
          double alphaVap = 1. - alphaLiq;
          double densiteLiq = m_run->m_eos[eLiq]->calculDensiteSaturation(pression, Tsat, Tsat);
          double densiteVap = m_run->m_eos[eVap]->calculDensiteSaturation(pression, Tsat, Tsat);

          //vitesse
          XMLElement *vitesse(donneesMelange->FirstChildElement("vitesse"));
          //if (vitesse == NULL) throw ErreurXMLElement("vitesse", nomFichier, __FILE__, __LINE__);
          double vitesseX(0.), vitesseY(0.), vitesseZ(0.);
          erreur = vitesse->QueryDoubleAttribute("x", &vitesseX);
          //if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut("x", nomFichier, __FILE__, __LINE__);
          erreur = vitesse->QueryDoubleAttribute("y", &vitesseY);
          //if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut("y", nomFichier, __FILE__, __LINE__);
          erreur = vitesse->QueryDoubleAttribute("z", &vitesseZ);
          //if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut("z", nomFichier, __FILE__, __LINE__);

          if (m_run->m_modele->quiSuisJe() == "EULERHOMOGENE") {
            etatsPhases.push_back(new PhaseEulerHomogene(alphaLiq, densiteLiq, pression, vitesseX, vitesseY, vitesseZ, m_run->m_eos[eLiq])); //Liq
            etatsPhases.push_back(new PhaseEulerHomogene(alphaVap, densiteVap, pression, vitesseX, vitesseY, vitesseZ, m_run->m_eos[eVap])); //Vap
          }
          else { throw ErreurXMLElement("Couplage Modele-MelangeBinaireLiqVap non valide", nomFichier.str(), __FILE__, __LINE__); }
        }

        //MATERIAU INCONNU
        else { throw ErreurXMLMateriauInconnu(typeMateriau, nomFichier.str(), __FILE__, __LINE__); } //Cas ou le type de materiau n a pas ete implemente
        materiau = materiau->NextSiblingElement("materiau");
      }
      if(nbMateriauxEtat!=m_run->m_nombrePhases) throw ErreurXMLEtat(etatDomaine, nomFichier.str(), __FILE__, __LINE__);

      //Etat melange selon modele
      if (m_run->m_modele->quiSuisJe() == "EULER") { etatMelange = new MelEuler(); }
      else if (m_run->m_modele->quiSuisJe() == "KAPILA") { etatMelange = new MelKapila(etat, nomFichier.str()); }
      else  if (m_run->m_modele->quiSuisJe() == "EULERHOMOGENE") { etatMelange = new MelEulerHomogene(materiau, nomFichier.str()); }
      else { throw ErreurXMLElement("Couplage Modele-Fluide non valide", nomFichier.str(), __FILE__, __LINE__); }

      //Lecture des variables transportees pour l'etat trouve
      vector<Transport> etatsTransport(m_run->m_nombreTransports);
      string nomTransport; double valeurTransport(0.);
      XMLElement *elementTransport(etat->FirstChildElement("transport"));
      int nbTransports(0);
      while (elementTransport != NULL) {
        nomTransport = elementTransport->Attribute("nom");
        elementTransport->QueryDoubleAttribute("valeur", &valeurTransport);
        int e(0);
        for (e = 0; e < m_run->m_nombreTransports; e++) {
          if (nomTransport == m_run->m_nomGTR[e]) { break; }
        }
        if (e != m_run->m_nombreTransports) {
          etatsTransport[e].setValeur(valeurTransport);
          nbTransports++;
        }
        elementTransport = elementTransport->NextSiblingElement("transport");
      }

      //C)Lecture type de domaine
      //*************************
      string typeDomaine(element->Attribute("type"));
      Outils::majuscule(typeDomaine);
      vector<string> NomsParametresDomaine;
      if      (typeDomaine == "DOMAINECOMPLET")  { domaines.push_back(new DGDomaineComplet(nomDomaine, etatsPhases, etatMelange, etatsTransport)); }
      else if (typeDomaine == "DEMIESPACE")      { domaines.push_back(new DGDemiEspace(nomDomaine, etatsPhases, etatMelange, etatsTransport, element, nomFichier.str())); }
      else if (typeDomaine == "DISQUE")          { domaines.push_back(new DGDisque(nomDomaine, etatsPhases, etatMelange, etatsTransport, element, nomFichier.str())); }
      else if (typeDomaine == "RECTANGLE")       { domaines.push_back(new DGRectangle(nomDomaine, etatsPhases, etatMelange, etatsTransport, element, nomFichier.str())); }
      else if (typeDomaine == "PAVE")            { domaines.push_back(new DGPave(nomDomaine, etatsPhases, etatMelange, etatsTransport, element, nomFichier.str())); }
      else if (typeDomaine == "SPHERE")          { domaines.push_back(new DGSphere(nomDomaine, etatsPhases, etatMelange, etatsTransport, element, nomFichier.str())); }
      else{ throw ErreurXMLDomaineInconnu(typeDomaine, nomFichier.str(), __FILE__, __LINE__); } //Cas ou le domaine n a pas ete implemente
      domaineTrouve++;
      //Domaine suivant
      element = element->NextSiblingElement("domaine");

      for (int k = 0; k < m_run->m_nombrePhases; k++) { delete etatsPhases[k]; }
      delete etatMelange;
    }

    //------------------------------ CONDITIONS AUX LIMITES -------------------------------

    //4) Recuperation racine CL du document XML
    //-----------------------------------------
    XMLElement *elementLimites(xmlNode->FirstChildElement("conditionsLimites"));
    if (elementLimites == NULL) throw ErreurXMLElement("conditionsLimites", nomFichier.str(), __FILE__, __LINE__);

    //5) Recuperation des condLim definies
    //------------------------------------
    string nomCondLim, etatLimite, typeCondLim;
    int numCondLim;
    int condLimTrouve(0);
    element = elementLimites->FirstChildElement("condLim");
    while (element != NULL)
    {
      //A)Lecture nom condLim
      //*********************
      nomCondLim = element->Attribute("nom");
      if (nomCondLim == "") throw ErreurXMLAttribut("nom", nomFichier.str(), __FILE__, __LINE__);
      Outils::majuscule(nomCondLim);

      //B)Lecture type condLim
      //**********************
      typeCondLim = element->Attribute("type");
      if (typeCondLim == "") throw ErreurXMLAttribut("type", nomFichier.str(), __FILE__, __LINE__);
      Outils::majuscule(typeCondLim);

      //C)numero de la condLim associe a la geometrie maillage
      //******************************************************
      erreur = element->QueryIntAttribute("numero", &numCondLim);
      if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut("numero", nomFichier.str(), __FILE__, __LINE__);

      //D)Lecture etat condition Limite si necessaire
      //*********************************************
      vector<Phase*> etats;
      if (element->Attribute("etat") != NULL) {
        etatLimite = element->Attribute("etat");
        if (etatLimite == "") throw ErreurXMLAttribut("etat", nomFichier.str(), __FILE__, __LINE__);
        Outils::majuscule(etatLimite);
        //Recherche de l etat associe au domaine
        XMLElement *etat(xmlNode->FirstChildElement("etat"));
        bool trouve(false); string nomEtat;
        while (etat != NULL)
        {
          nomEtat = etat->Attribute("nom");
          Outils::majuscule(nomEtat);
          if (nomEtat == etatLimite) { trouve = true; break; }
          etat = etat->NextSiblingElement("etat");
        }
        if (!trouve) { throw ErreurXMLEtat(etatLimite, nomFichier.str(), __FILE__, __LINE__); }
        //On peut lire chaque materiau de l etat trouve
        string typeMateriau; string nomEOS;
        XMLElement* materiau(etat->FirstChildElement("materiau"));
        while (materiau != NULL)
        {
          //Recuperation de l EOS
          nomEOS = materiau->Attribute("EOS");
          //Outils::majuscule(nomEOS);
          int e(0);  bool eosTrouvee(false);
          for (e = 0; e < m_run->m_nombrePhases; e++){
            if (nomEOS == m_run->m_eos[e]->retourneNom()) { break; }
          }
          if (e == m_run->m_nombrePhases) { throw ErreurXMLEOSInconnue(nomEOS, nomFichier.str(), __FILE__, __LINE__); }
          typeMateriau = materiau->Attribute("type");
          Outils::majuscule(typeMateriau);
          if (typeMateriau == "FLUIDE") { 
            if (m_run->m_modele->quiSuisJe() == "EULER") { etats.push_back(new PhaseEuler(materiau, m_run->m_eos[e], nomFichier.str())); }
            else if (m_run->m_modele->quiSuisJe() == "KAPILA") { etats.push_back(new PhaseKapila(materiau, m_run->m_eos[e], nomFichier.str())); }
            else { throw ErreurXMLElement("Couplage Modele-Fluide non valide", nomFichier.str(), __FILE__, __LINE__); }
          }
          else { throw ErreurXMLMateriauInconnu(typeMateriau, nomFichier.str(), __FILE__, __LINE__); } //Cas ou le type de materiau n a pas ete implemente
          materiau = materiau->NextSiblingElement("materiau");
        }
      }

      //E)Lecture du type de condLim
      //****************************
      if      (typeCondLim == "ABS") { condLim.push_back(new CondLimAbs(numCondLim)); }
      else if (typeCondLim == "MUR") { if (m_run->m_ordre == "ORDRE1") { condLim.push_back(new CondLimMur(numCondLim)); } else { condLim.push_back(new CondLimMurO2(numCondLim)); } }
      else if (typeCondLim == "INJECTION") { condLim.push_back(new CondLimInj(numCondLim, element, etats, m_run->m_nombreTransports, m_run->m_nomGTR, nomFichier.str())); }
      else if (typeCondLim == "RESERVOIR") { condLim.push_back(new CondLimRes(numCondLim, element, etats, m_run->m_nombreTransports, m_run->m_nomGTR, nomFichier.str())); }
      else if (typeCondLim == "SORTIE") { condLim.push_back(new CondLimSortie(numCondLim, element, m_run->m_nombrePhases, m_run->m_nombreTransports, m_run->m_nomGTR, nomFichier.str())); }
      else { throw ErreurXMLCondLimInconnue(typeCondLim, nomFichier.str(), __FILE__, __LINE__); } //Cas ou la limite n a pas ete implemente
      condLimTrouve++;
      //Domaine suivant
      element = element->NextSiblingElement("condLim");

    }

  }
  catch (ErreurXML &){ throw; } // Renvoi au niveau suivant
}

//***********************************************************************
