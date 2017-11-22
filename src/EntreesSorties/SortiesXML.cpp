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

#include "SortiesXML.h"
#include "../Run.h"

using namespace std;
using namespace tinyxml2;

//***********************************************************************

SortiesXML::SortiesXML(){}

//***********************************************************************

SortiesXML::SortiesXML(string casTest, string run, XMLElement *element, string nomFichier, Entrees *entree) :
  Sorties(casTest, run, element, nomFichier, entree)
{}

//***********************************************************************

SortiesXML::~SortiesXML(){}

//***********************************************************************

void SortiesXML::prepareSortieSpecifique()
{
  try {
    ofstream fluxFichier;
    //Création du fichier de sortie collection
    m_fichierCollection = m_dossierSortie + creationNomFichierXML(m_nomFichierCollection.c_str());
    fluxFichier.open(m_fichierCollection.c_str());
    if (!fluxFichier) { throw ErreurECOGEN("Impossible d ouvrir le fichier " + m_fichierCollection, __FILE__, __LINE__); }
    fluxFichier << "<?xml version=\"1.0\"?>" << endl;
    fluxFichier << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"";
    if (!m_ecritBinaire) { fluxFichier << "LittleEndian\" "; }
    else { fluxFichier << m_endianMode.c_str() << "\" "; }
    fluxFichier << "compressor=\"vtkZLibDataCompressor\">";
    fluxFichier << endl << "  <Collection>" << endl;
    fluxFichier << endl << "  </Collection>" << endl;
    fluxFichier << "</VTKFile>" << endl;
    fluxFichier.close();
  }
  catch (ErreurECOGEN &) { throw; }
}

//***********************************************************************

void SortiesXML::ecritSolutionSpecifique(Maillage* maillage, vector<Cellule *> *cellulesLvl)
{
  try {
    //Ecriture des fichiers de sortie au format XML
    ecritSolutionXML(maillage, cellulesLvl);
    //Ajout du fichier Collection pour grouper les niveaux, les temps, les CPU, etc.
    if (rang == 0) { ecritCollectionXML(maillage); }
  }
  catch (ErreurECOGEN &) { throw; } // Renvoi au niveau suivant
}

//***********************************************************************

string SortiesXML::creationNomFichierXML(const char* nom, Maillage *maillage, int lvl, int proc, int numFichier, string nomVariable)
{
  stringstream num;
  string prefix;
  if (proc==-2) { prefix = "p"; }
  else { prefix = ""; }

  try {
    //FP//TODO Donnees separees a faire...
    if (m_donneesSeparees) { throw ErreurECOGEN("SortiesXML::creationNomFichierXML : donnees Separees non prevu", __FILE__, __LINE__); }
    num << nom;
    //Gestion nomVariable
    if (nomVariable != "defaut") num << "_" << nomVariable << "_";
    //Gestion binaire
    if (m_ecritBinaire) num << "B64";
    //Gestion cpu
    if (proc > -1) num << "_CPU" << proc;
    //Gestion niveau AMR
    //if (lvl != -1 && maillage->getType()==AMR) 
    if (lvl != -1) num << "_AMR" << lvl;
    //Gestion numero de fichier resultat
    if (numFichier != -1) num << "_TIME" << numFichier;
    //Gestion extension
    if(maillage==0) num << ".pvd"; //Extension pour la collection
    else {
      num << "." << prefix;
      switch (maillage->getType()) {
      case REC:
        num << "vtr"; break;
      case UNS:
        num << "vtu"; break;
      case AMR:
        num << "vtp"; break;
      default:
        throw ErreurECOGEN("SortiesXML::creationNomFichierXML : type maillage inconnu", __FILE__, __LINE__);
      }
    }
  }
  catch (ErreurECOGEN &) { throw; }

  return num.str();
}

//***********************************************************************

void SortiesXML::ecritSolutionXML(Maillage* maillage, vector<Cellule *> *cellulesLvl)
{
  ofstream fluxFichier;

  try {
    //On balaye les niveau pour AMR
    for (int lvl = 0; lvl <= maillage->getLvlMax(); lvl++) {

      //1) Ouverture / creation fichier
      //-------------------------------
      string fichier = m_dossierSortie + creationNomFichierXML(m_nomFichierResultats.c_str(), maillage, lvl, rang, m_numFichier);
      fluxFichier.open(fichier.c_str());
      if (!fluxFichier) { throw ErreurECOGEN("Impossible d ouvrir le fichier " + fichier, __FILE__, __LINE__); }
      fluxFichier << "<?xml version=\"1.0\"?>" << endl;

      //2) Ecriture du maillage
      //-----------------------
      switch (maillage->getType()) {
      case REC:
        ecritMaillageRectilinearXML(maillage, cellulesLvl, fluxFichier); break;
      case UNS:
        ecritMaillageUnstructuredXML(maillage, cellulesLvl, fluxFichier, lvl); break;
      case AMR:
        ecritMaillagePolyDataXML(maillage, cellulesLvl, fluxFichier, lvl); break;
      default:
        throw ErreurECOGEN("Sorties::ecritSolutionXML : type maillage inconnu", __FILE__, __LINE__); break;
      }

      //3) Ecriture des donnees phases fluides
      //--------------------------------------
      ecritDonneesPhysiquesXML(maillage, cellulesLvl, fluxFichier, lvl);

      //4) Finalisation fichier
      //-----------------------
      switch (maillage->getType()) {
      case REC:
        ecritFinFichierRectilinearXML(fluxFichier); break;
      case UNS:
        ecritFinFichierUnstructuredXML(fluxFichier); break;
      case AMR:
        ecritFinFichierPolyDataXML(fluxFichier); break;
      default:
        throw ErreurECOGEN("Sorties::ecritSolutionXML : type maillage inconnu", __FILE__, __LINE__); break;
      }
      fluxFichier.close();

    } //Fin lvl
  } //Fin try
  catch (ErreurECOGEN &) { throw; }
}

//***********************************************************************

void SortiesXML::ecritCollectionXML(Maillage *maillage)
{
  try {
    //1) Parsing du fichier XML par la bibliotheque tinyxml2
    //------------------------------------------------------
    XMLDocument xmlCollection;
    XMLError erreur(xmlCollection.LoadFile(m_fichierCollection.c_str())); //Le fichier est parse ici
    if (erreur != XML_SUCCESS) {
      throw ErreurXML(m_fichierCollection, __FILE__, __LINE__);
    }
    XMLNode *VTKFile = xmlCollection.FirstChildElement("VTKFile");
    if (VTKFile == NULL) throw ErreurXMLRacine("VTKFile", m_fichierCollection, __FILE__, __LINE__);
    XMLElement *collection = VTKFile->FirstChildElement("Collection");
    if (collection == NULL) throw ErreurXMLElement("Collection", m_fichierCollection, __FILE__, __LINE__);
    //2) Insertion des fichiers resultats
    //-----------------------------------
    for (int p = 0; p < Ncpu; p++) {
      for (int lvl = 0; lvl <= maillage->getLvlMax(); lvl++) {
        XMLElement *dataSet = xmlCollection.NewElement("DataSet");
        if (dataSet == NULL) throw ErreurXMLElement("DataSet", m_fichierCollection, __FILE__, __LINE__);
        dataSet->SetAttribute("timestep", m_numFichier);
        dataSet->SetAttribute("part", p);
        //dataSet->SetAttribute("group", p);
        //Nom du fichier a jouter
        string fichier = creationNomFichierXML(m_nomFichierResultats.c_str(), maillage, lvl, p, m_numFichier);
        dataSet->SetAttribute("file", fichier.c_str());
        collection->InsertEndChild(dataSet);
      }
    }
    xmlCollection.SaveFile(m_fichierCollection.c_str());
  }
  catch (ErreurXML &) { throw; } // Renvoi au niveau suivant
}

//***********************************************************************

void SortiesXML::ecritDonneesPhysiquesXML(Maillage *maillage, vector<Cellule *> *cellulesLvl, std::ofstream &fluxFichier, const int &lvl, bool parallel)
{
  vector<double> jeuDonnees;

  string prefix;
  if (parallel) { prefix = "P"; }
  else { prefix = ""; }

  string format = "ascii";
  if (m_ecritBinaire) format = "binary";

  fluxFichier << "      <" << prefix << "CellData>" << endl;

  //1) Ecriture des variables des phases
  //------------------------------------
  for (int phase = 0; phase < m_run->getNombrePhases(); phase++)
  {
    //Ecriture des variables scalaires
    for (int var = 1; var <= m_celluleRef.getPhase(phase)->getNombreScalaires(); var++) {
      fluxFichier << "        <" << prefix << "DataArray type=\"Float32\" Name=\"F" << phase << "_" << m_celluleRef.getPhase(phase)->renvoieNomScalaire(var) << "_" << m_celluleRef.getPhase(phase)->getEos()->retourneNom() << "\"";
      if (!parallel) {
        fluxFichier << " format=\"" << format << "\">" << endl;
        jeuDonnees.clear();
        maillage->recupereDonnees(cellulesLvl, jeuDonnees, var, phase, lvl);
        this->ecritJeuDonnees(jeuDonnees, fluxFichier, FLOAT);
        fluxFichier << endl;
        fluxFichier << "        </" << prefix << "DataArray>" << endl;
      }
      else { fluxFichier << "\"/>" << endl; }
    }
    //Ecriture des variables vectorielles
    for (int var = 1; var <= m_celluleRef.getPhase(phase)->getNombreVecteurs(); var++)
    {
      fluxFichier << "        <" << prefix << "DataArray type=\"Float32\" Name=\"F" << phase << "_" << m_celluleRef.getPhase(phase)->renvoieNomVecteur(var) << "_" << m_celluleRef.getPhase(phase)->getEos()->retourneNom() << "\" NumberOfComponents=\"3\"";
      if (!parallel) {
        fluxFichier << " format=\"" << format << "\">" << endl;
        jeuDonnees.clear();
        maillage->recupereDonnees(cellulesLvl, jeuDonnees, -var, phase, lvl);
        this->ecritJeuDonnees(jeuDonnees, fluxFichier, FLOAT);
        fluxFichier << endl;
        fluxFichier << "        </" << prefix << "DataArray>" << endl;
      }
      else { fluxFichier << "\"/>" << endl; }
    }
  } //Fin phase

  //2) Ecriture des donnees melange
  //-------------------------------
  if (m_run->m_nombrePhases > 1)
  {
    int melange = -1;
    //Ecriture des variables scalaires du melange
    for (int var = 1; var <= m_celluleRef.getMelange()->getNombreScalaires(); var++)
    {
      fluxFichier << "        <" << prefix << "DataArray type=\"Float32\" Name=\"" << m_celluleRef.getMelange()->renvoieNomScalaire(var) << "\"";
      if (!parallel) {
        fluxFichier << " format=\"" << format << "\">" << endl;
        jeuDonnees.clear();
        maillage->recupereDonnees(cellulesLvl, jeuDonnees, var, melange, lvl);
        this->ecritJeuDonnees(jeuDonnees, fluxFichier, FLOAT);
        fluxFichier << endl;
        fluxFichier << "        </" << prefix << "DataArray>" << endl;
      }
      else { fluxFichier << "\"/>" << endl; }
    }
    //Ecriture des variables vectorielles du melange
    for (int var = 1; var <= m_celluleRef.getMelange()->getNombreVecteurs(); var++)
    {
      fluxFichier << "        <" << prefix << "DataArray type=\"Float32\" Name=\"" << m_celluleRef.getMelange()->renvoieNomVecteur(var) << "\" NumberOfComponents=\"3\"";
      if (!parallel) {
        fluxFichier << " format=\"" << format << "\">" << endl;
        jeuDonnees.clear();
        maillage->recupereDonnees(cellulesLvl, jeuDonnees, -var, melange, lvl);
        this->ecritJeuDonnees(jeuDonnees, fluxFichier, FLOAT);
        fluxFichier << endl;
        fluxFichier << "        </" << prefix << "DataArray>" << endl;
      }
      else { fluxFichier << "\"/>" << endl; }
    }
  } //Fin melange

  //3) Ecriture des transports et autres...
  //---------------------------------------
  int transport = -2;
  for (int var = 1; var <= m_run->m_nombreTransports; var++)
  {
    fluxFichier << "        <" << prefix << "DataArray type=\"Float32\" Name=\"T" << var << "\"";
    if (!parallel) {
      fluxFichier << " format=\"" << format << "\">" << endl;
      jeuDonnees.clear();
      maillage->recupereDonnees(cellulesLvl, jeuDonnees, var, transport, lvl);
      this->ecritJeuDonnees(jeuDonnees, fluxFichier, FLOAT);
      fluxFichier << endl;
      fluxFichier << "        </" << prefix << "DataArray>" << endl;
    }
    else { fluxFichier << "\"/>" << endl; }
  }

  //4) Ecriture indicateur xi
  //-------------------------
  int xi = -3;
  fluxFichier << "        <" << prefix << "DataArray type=\"Float32\" Name=\"Xi\"";
  if (!parallel) {
    fluxFichier << " format=\"" << format << "\">" << endl;
    jeuDonnees.clear();
    maillage->recupereDonnees(cellulesLvl, jeuDonnees, 1, xi, lvl);
    this->ecritJeuDonnees(jeuDonnees, fluxFichier, FLOAT);
    fluxFichier << endl;
    fluxFichier << "        </" << prefix << "DataArray>" << endl;
  }
  else { fluxFichier << "\"/>" << endl; }

  //5) Ecriture gradient rho
  //------------------------
  int gradRho = -4;
  fluxFichier << "        <" << prefix << "DataArray type=\"Float32\" Name=\"gradRho\"";
  if (!parallel) {
    fluxFichier << " format=\"" << format << "\">" << endl;
    jeuDonnees.clear();
    maillage->recupereDonnees(cellulesLvl, jeuDonnees, 1, gradRho, lvl);
    this->ecritJeuDonnees(jeuDonnees, fluxFichier, FLOAT);
    fluxFichier << endl;
    fluxFichier << "        </" << prefix << "DataArray>" << endl;
  }
  else { fluxFichier << "\"/>" << endl; }

  //Fin
  fluxFichier << "      </" << prefix << "CellData>" << endl;
}

//***********************************************************************

void SortiesXML::ecritMaillageRectilinearXML(Maillage *maillage, vector<Cellule *> *cellulesLvl, std::ofstream &fluxFichier, bool parallel)
{
  vector<double> jeuDonnees;

  string prefix;
  if (parallel) { prefix = "P"; }
  else { prefix = ""; }

  //0) Entete
  //---------
  fluxFichier << "<VTKFile type=\"" << prefix << "RectilinearGrid\" version=\"0.1\" byte_order=\"";
  if (!m_ecritBinaire) fluxFichier << "LittleEndian\">" << endl;
  else fluxFichier << m_endianMode.c_str() << "\">" << endl;
  if (!parallel) {
    fluxFichier << "  <RectilinearGrid WholeExtent=\"" << maillage->recupereChaineExtent(rang) << "\">" << endl;
    fluxFichier << "    <Piece Extent=\"" << maillage->recupereChaineExtent(rang) << "\">" << endl;
  }
  else {
    fluxFichier << "  <PRectilinearGrid WholeExtent = \"" << maillage->recupereChaineExtent(rang, true) << "\" GhostLevel=\"0\">" << endl;
  }

  //1) Ecriture des Coordonnees des noeuds
  //--------------------------------------
  fluxFichier << "      <" << prefix << "Coordinates>" << endl;
  //Coordonnees en X
  fluxFichier << "        <" << prefix << "DataArray type=\"Float32\" ";
  if (!parallel) {
    if (!m_ecritBinaire) { fluxFichier << "format=\"ascii\">" << endl << "          "; }
    else { fluxFichier << "format=\"binary\">" << endl; }
    jeuDonnees.clear();
    maillage->recupereCoord(cellulesLvl, jeuDonnees, X);
    ecritJeuDonnees(jeuDonnees, fluxFichier, FLOAT);
    fluxFichier << endl;
    fluxFichier << "        </" << prefix << "DataArray>" << endl;
  }
  else { fluxFichier << "/>" << endl; }
  //Coordonnees en Y
  fluxFichier << "        <" << prefix << "DataArray type=\"Float32\" ";
  if (!parallel) {
    if (!m_ecritBinaire) { fluxFichier << "format=\"ascii\">" << endl << "          "; }
    else { fluxFichier << "format=\"binary\">" << endl; }
    jeuDonnees.clear();
    maillage->recupereCoord(cellulesLvl, jeuDonnees, Y);
    ecritJeuDonnees(jeuDonnees, fluxFichier, FLOAT);
    fluxFichier << endl;
    fluxFichier << "        </" << prefix << "DataArray>" << endl;
  }
  else { fluxFichier << "/>" << endl; }
  //Coordonnees en Z
  fluxFichier << "        <" << prefix << "DataArray type=\"Float32\" ";
  if (!parallel) {
    if (!m_ecritBinaire) { fluxFichier << "format=\"ascii\">" << endl << "          "; }
    else { fluxFichier << "format=\"binary\">" << endl; }
    jeuDonnees.clear();
    maillage->recupereCoord(cellulesLvl, jeuDonnees, Z);
    ecritJeuDonnees(jeuDonnees, fluxFichier, FLOAT);
    fluxFichier << endl;
    fluxFichier << "        </" << prefix << "DataArray>" << endl;
  }
  else { fluxFichier << "/>" << endl; }
  fluxFichier << "      </" << prefix << "Coordinates>" << endl;
}

//***********************************************************************

void SortiesXML::ecritMaillageUnstructuredXML(Maillage *maillage, vector<Cellule *> *cellulesLvl, std::ofstream &fluxFichier, const int &lvl, bool parallel)
{
  vector<double> jeuDonnees;

  string prefix;
  if (parallel) { prefix = "P"; }
  else { prefix = ""; }

  //0) Entete
  //---------
  fluxFichier << "<VTKFile type=\"" << prefix << "UnstructuredGrid\" version=\"0.1\" byte_order=\"";
  if (!m_ecritBinaire) fluxFichier << "LittleEndian\">" << endl;
  else fluxFichier << m_endianMode.c_str() << "\">" << endl;

  if (parallel) {
    fluxFichier << "  <PUnstructuredGrid GhostLevel=\"0\">" << endl;
  }
  else {
    fluxFichier << "  <UnstructuredGrid>" << endl;
    maillage->ecritEntetePiece(fluxFichier, cellulesLvl, lvl);
  }

  //1) Ecriture des Noeuds
  //----------------------
  fluxFichier << "      <" << prefix << "Points>" << endl;
  fluxFichier << "        <" << prefix << "DataArray type=\"Float32\" NumberOfComponents=\"3\" ";
  if (!parallel) {
    if (!m_ecritBinaire) { fluxFichier << "format=\"ascii\">" << endl << "          "; }
    else { fluxFichier << "format=\"binary\">" << endl; }
    jeuDonnees.clear();
    maillage->recupereNoeuds(jeuDonnees, lvl);
    ecritJeuDonnees(jeuDonnees, fluxFichier, FLOAT);
    fluxFichier << endl;
    fluxFichier << "        </" << prefix << "DataArray>" << endl;
  }
  else { fluxFichier << "/>" << endl; }
  fluxFichier << "      </" << prefix << "Points>" << endl;

  //2) Ecriture des Cellules
  //------------------------
  fluxFichier << "      <" << prefix << "Cells>" << endl;
  //Connectivite
  fluxFichier << "        <" << prefix << "DataArray type=\"Int32\" Name=\"connectivity\" ";
  if (!parallel) {
    if (!m_ecritBinaire) { fluxFichier << "format=\"ascii\">" << endl << "          "; }
    else { fluxFichier << "format=\"binary\">" << endl; }
    jeuDonnees.clear();
    maillage->recupereConnectivite(jeuDonnees, lvl);
    ecritJeuDonnees(jeuDonnees, fluxFichier, INT);
    fluxFichier << endl;
    fluxFichier << "        </" << prefix << "DataArray>" << endl;
  }
  else { fluxFichier << "/>" << endl; }
  //Offsets
  fluxFichier << "        <" << prefix << "DataArray type=\"Int32\" Name=\"offsets\" ";
  if (!parallel) {
    if (!m_ecritBinaire) { fluxFichier << "format=\"ascii\">" << endl << "          "; }
    else { fluxFichier << "format=\"binary\">" << endl; }
    jeuDonnees.clear();
    maillage->recupereOffsets(jeuDonnees, lvl);
    ecritJeuDonnees(jeuDonnees, fluxFichier, INT);
    fluxFichier << endl;
    fluxFichier << "        </" << prefix << "DataArray>" << endl;
  }
  else { fluxFichier << "/>" << endl; }
  //Type de cellules
  fluxFichier << "        <" << prefix << "DataArray type=\"UInt8\" Name=\"types\" ";
  if (!parallel) {
    if (!m_ecritBinaire) { fluxFichier << "format=\"ascii\">" << endl << "          "; }
    else { fluxFichier << "format=\"binary\">" << endl; }
    jeuDonnees.clear();
    maillage->recupereTypeCell(jeuDonnees, lvl);
    ecritJeuDonnees(jeuDonnees, fluxFichier, CHAR);
    fluxFichier << endl;
    fluxFichier << "        </" << prefix << "DataArray>" << endl;
  }
  else { fluxFichier << "/>" << endl; }
  fluxFichier << "      </" << prefix << "Cells>" << endl;
}

//***********************************************************************

void SortiesXML::ecritMaillagePolyDataXML(Maillage *maillage, vector<Cellule *> *cellulesLvl, std::ofstream &fluxFichier, const int &lvl, bool parallel)
{
  vector<double> jeuDonnees;

  string prefix;
  if (parallel) { prefix = "P"; }
  else { prefix = ""; }

  //0) Entete
  //---------
  fluxFichier << "<VTKFile type=\"" << prefix << "PolyData\" version=\"0.1\" byte_order=\"";
  if (!m_ecritBinaire) fluxFichier << "LittleEndian\">" << endl;
  else fluxFichier << m_endianMode.c_str() << "\">" << endl;

  if (parallel) {
    fluxFichier << "  <PPolyData GhostLevel=\"0\">" << endl;
  }
  else {
    fluxFichier << "  <PolyData>" << endl;
    maillage->ecritEntetePiece(fluxFichier, cellulesLvl, lvl);
  }

  //1) Ecriture des Noeuds
  //----------------------
  fluxFichier << "      <" << prefix << "Points>" << endl;
  fluxFichier << "        <" << prefix << "DataArray type=\"Float32\" NumberOfComponents=\"3\" ";
  if (!parallel) {
    if (!m_ecritBinaire) { fluxFichier << "format=\"ascii\">" << endl << "          "; }
    else { fluxFichier << "format=\"binary\">" << endl; }
    jeuDonnees.clear();
    maillage->recupereNoeuds(jeuDonnees, lvl);
    ecritJeuDonnees(jeuDonnees, fluxFichier, FLOAT);
    fluxFichier << endl;
    fluxFichier << "        </" << prefix << "DataArray>" << endl;
  }
  else { fluxFichier << "/>" << endl; }
  fluxFichier << "      </" << prefix << "Points>" << endl;

  //2) Ecriture des Polys
  //---------------------
  fluxFichier << "      <" << prefix << "Polys>" << endl;
  //Connectivite
  fluxFichier << "        <" << prefix << "DataArray type=\"Int32\" Name=\"connectivity\" ";
  if (!parallel) {
    if (!m_ecritBinaire) { fluxFichier << "format=\"ascii\">" << endl << "          "; }
    else { fluxFichier << "format=\"binary\">" << endl; }
    jeuDonnees.clear();
    maillage->recupereConnectivite(jeuDonnees, lvl);
    ecritJeuDonnees(jeuDonnees, fluxFichier, INT);
    fluxFichier << endl;
    fluxFichier << "        </" << prefix << "DataArray>" << endl;
  }
  else { fluxFichier << "/>" << endl; }
  //Offsets
  fluxFichier << "        <" << prefix << "DataArray type=\"Int32\" Name=\"offsets\" ";
  if (!parallel) {
    if (!m_ecritBinaire) { fluxFichier << "format=\"ascii\">" << endl << "          "; }
    else { fluxFichier << "format=\"binary\">" << endl; }
    jeuDonnees.clear();
    maillage->recupereOffsets(jeuDonnees, lvl);
    ecritJeuDonnees(jeuDonnees, fluxFichier, INT);
    fluxFichier << endl;
    fluxFichier << "        </" << prefix << "DataArray>" << endl;
  }
  else { fluxFichier << "/>" << endl; }
  fluxFichier << "      </" << prefix << "Polys>" << endl;
}

//***********************************************************************

void SortiesXML::ecritFinFichierRectilinearXML(std::ofstream &fluxFichier, bool parallel)
{
  string prefix;
  if (parallel) { prefix = "P"; }
  else { prefix = ""; }

  if (!parallel) fluxFichier << "    </Piece>" << endl;
  fluxFichier << "  </" << prefix << "RectilinearGrid>" << endl;
  fluxFichier << "</VTKFile>" << endl;
}

//***********************************************************************

void SortiesXML::ecritFinFichierUnstructuredXML(std::ofstream &fluxFichier, bool parallel)
{
  string prefix;
  if (parallel) { prefix = "P"; }
  else { prefix = ""; }

  if (!parallel) fluxFichier << "    </Piece>" << endl;
  fluxFichier << "  </" << prefix << "UnstructuredGrid>" << endl;
  fluxFichier << "</VTKFile>" << endl;
}

//***********************************************************************

void SortiesXML::ecritFinFichierPolyDataXML(std::ofstream &fluxFichier, bool parallel)
{
  string prefix;
  if (parallel) { prefix = "P"; }
  else { prefix = ""; }

  if (!parallel) fluxFichier << "    </Piece>" << endl;
  fluxFichier << "  </" << prefix << "PolyData>" << endl;
  fluxFichier << "</VTKFile>" << endl;
}

//***********************************************************************

//Old

//***********************************************************************

void SortiesXML::ecritFichierParallelXML(Maillage *maillage, vector<Cellule *> *cellulesLvl)
{
  ofstream fluxFichier;
  stringstream num;

  bool parallel(true);

  try {
    //Preparation AMR
    for (int lvl = 0; lvl <= maillage->getLvlMax(); lvl++) {

      //1) Ouverture / creation fichier
      //-------------------------------
      if (!m_donneesSeparees) {
        if (!m_ecritBinaire) { num << m_dossierSortie << "/result_" << m_numFichier; }
        else { num << m_dossierSortie << "/resultB64_" << m_numFichier; }
        num << "_lvl" << lvl;
      }
      else { throw ErreurECOGEN("Sorties::ecritFichierParallelXML : donnees Separees non prevu", __FILE__, __LINE__); }
      switch (maillage->getType()) {
      case REC:
        num << ".pvtr"; break;
      case UNS:
        num << ".pvtu"; break;
      case AMR:
        num << ".pvtp"; break;
      default:
        throw ErreurECOGEN("Sorties::ecritFichierParallelXML : type maillage inconnu", __FILE__, __LINE__); break;
      }

      fluxFichier.open((num.str()).c_str());
      if (!fluxFichier) { throw ErreurECOGEN("Impossible d ouvrir le fichier " + num.str(), __FILE__, __LINE__); }
      fluxFichier << "<?xml version=\"1.0\"?>" << endl;

      //2) Ecriture des infos maillage
      //------------------------------
      switch (maillage->getType()) {
      case REC:
        ecritMaillageRectilinearXML(maillage, cellulesLvl, fluxFichier, parallel); break;
      case UNS:
        ecritMaillageUnstructuredXML(maillage, cellulesLvl, fluxFichier, lvl, parallel); break;
      case AMR:
        ecritMaillagePolyDataXML(maillage, cellulesLvl, fluxFichier, lvl, parallel); break;
      default:
        throw ErreurECOGEN("Sorties::ecritSolutionXML : type maillage inconnu", __FILE__, __LINE__); break;
      }

      //3) Ecriture des donnees phases fluides
      //--------------------------------------
      ecritDonneesPhysiquesXML(maillage, cellulesLvl, fluxFichier, lvl, parallel);

      //4) Ecriture des nom des fichiers
      //--------------------------------
      for (int p = 0; p < Ncpu; p++)
      {
        stringstream fichTemp;
        if (m_ecritBinaire) { fichTemp << "resultB64_"; }
        else { fichTemp << "result_"; }
        //fichTemp << m_numFichier << "_lvl" << lvl << "_CPU" << rang;
        fichTemp << "CPU" << rang << "_lvl" << lvl << "_" << m_numFichier;

        switch (maillage->getType()) {
        case REC:
          fichTemp << ".vtr";
          fluxFichier << "    <Piece Extent=\"" << maillage->recupereChaineExtent(p) << "\" Source=\"" << fichTemp.str() << "\"/>" << endl;
          break;
        case UNS:
          fichTemp << ".vtu";
          fluxFichier << "    <Piece Source=\"" << fichTemp.str() << "\"/>" << endl;
          break;
        case AMR:
          fichTemp << ".vtp";
          fluxFichier << "    <Piece Source=\"" << fichTemp.str() << "\"/>" << endl;
          break;
        default:
          throw ErreurECOGEN("Sorties::ecritSolutionXML : type maillage inconnu", __FILE__, __LINE__); break;
        }
      }

      //5) Finalisation fichier
      //-----------------------
      switch (maillage->getType()) {
      case REC:
        ecritFinFichierRectilinearXML(fluxFichier, parallel); break;
      case UNS:
        ecritFinFichierUnstructuredXML(fluxFichier, parallel); break;
      case AMR:
        ecritFinFichierPolyDataXML(fluxFichier, parallel); break;
      default:
        throw ErreurECOGEN("Sorties::ecritSolutionXML : type maillage inconnu", __FILE__, __LINE__); break;
      }
      fluxFichier.close();

    }

  } //Fin try
  catch (ErreurECOGEN &) { throw; }
}

//***********************************************************************
