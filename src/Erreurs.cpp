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

#include "Erreurs.h"
#include "Run.h"

using namespace std;

std::vector<Erreurs> erreurs;

//***********************************************************************

Erreurs::Erreurs() : m_etat(0), m_ligne(0), m_valeur(0.), m_message("non renseigne")
{}

//***********************************************************************

Erreurs::Erreurs(const string &message, const char* fichierSource, int numeroLigne) :
  m_message(message), m_fichier(fichierSource), m_ligne(numeroLigne), m_valeur(0.)
{
  m_etat = 1;
}

//***********************************************************************

Erreurs::~Erreurs(){}

//***********************************************************************

void Erreurs::messageErreur(const string &message)
{
  if (rang == 0) 
  { 
    stringstream numCPU;
    numCPU << rang;

    cerr << endl << "-------------------------------------------------------" << endl;
    cerr << "ERREUR sur CPU " + numCPU.str() + " : " << message.c_str() << endl; 
    cerr << "Corriger erreur et relancer le code" << endl;
    cerr << "-------------------------------------------------------" << endl;
  }
  //run.finalise();
  exit(0);
}

//***********************************************************************

void Erreurs::messageErreur(const std::string &message, double valeur)
{
  if (rang == 0)
  {
    stringstream numCPU;
    numCPU << rang;

    cout << endl << "-------------------------------------------------------" << endl;
    cout << "ERREUR sur CPU " + numCPU.str() + " : " << message.c_str() << " " << valeur << endl;
    cout << "Corriger erreur et relancer le code" << endl;
    cout << "-------------------------------------------------------" << endl;
  }
  //run.finalise();
  exit(0);
}

//***********************************************************************

void Erreurs::setErreur(const std::string &message, const char* fichierSource, int numeroLigne)
{
  m_message = message;
  m_etat = 1;
  m_fichier = fichierSource;
  m_ligne = numeroLigne;
}

//***********************************************************************

void Erreurs::setErreur(const std::string &message, const double valeur)
{
  m_message = message;
  m_etat = 1;
  m_valeur = valeur;
}

//***********************************************************************

int Erreurs::getEtat()
{
  return m_etat;
}

//***********************************************************************

void Erreurs::afficheErreur()
{
  cout << endl << "-------------------------------------------------------" << endl;
  cout << "        ERREUR NECESSITANT ARRET DU PROGRAMME" << endl;
  cout << " - numero du CPU : " << rang << endl;
  cout << " - fichier :  " << m_fichier.c_str() << " ligne : " << m_ligne << endl;
  cout << " - message : " << m_message.c_str() << " " << m_valeur << endl;
  cout << "=> Corriger erreur et relancer le run" << endl;
  cout << "-------------------------------------------------------" << endl;
}

//***********************************************************************

void Erreurs::ecritErreurFichier()
{
  ofstream fluxFichier;
  stringstream num;

  num << "erreur_CPU" << rang;
  fluxFichier.open((num.str()).c_str());

  fluxFichier << "-------------------------------------------------------" << endl;
  fluxFichier << "        ERREUR NECESSITANT ARRET DU PROGRAMME" << endl;
  fluxFichier << " - numero du CPU : " << rang << endl;
  fluxFichier << " - fichier :  " << m_fichier.c_str() << " ligne : " << m_ligne << endl;
  fluxFichier << " - message : " << m_message.c_str() << " " << m_valeur << endl;
  fluxFichier << "=> Corriger erreur et relancer le run" << endl;
  fluxFichier << "-------------------------------------------------------" << endl;

  fluxFichier.close();
}

//***********************************************************************

void Erreurs::arretCodeApresErreur(vector<Erreurs> &erreurs)
{
  //Affichage ecran des erreurs
  for (unsigned int e = 0; e < erreurs.size(); e++) {
    erreurs[e].afficheErreur();
  }
  //Arret propre
  //run.finalise();
  exit(0);
}

//***********************************************************************
