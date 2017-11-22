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

#include "CondLimSortie.h"

using namespace std;
using namespace tinyxml2;

//****************************************************************************

CondLimSortie::CondLimSortie(){}

//****************************************************************************

CondLimSortie::CondLimSortie(int numPhysique, XMLElement *element, int &nombrePhases, int &nombreTransports, std::vector<std::string> nomTransports, string nomFichier) :
  CondLim(numPhysique)
{
  //Lecture de la pression en sortie
  XMLElement *sousElement(element->FirstChildElement("donneesSortie"));
  if (sousElement == NULL) throw ErreurXMLElement("donneesSortie", nomFichier, __FILE__, __LINE__);
  //Recuperation des attributs
  //--------------------------
  XMLError erreur;
  erreur = sousElement->QueryDoubleAttribute("p0", &m_p0);
  if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut("p0", nomFichier, __FILE__, __LINE__);

  //Lecture des transports
  int couleurTrouvee(0);
  m_valeurTransport = new double[nombreTransports];
  XMLElement *elementTransport(sousElement->FirstChildElement("transport"));
  string nomTransport;
  while (elementTransport != NULL)
  {
    nomTransport = elementTransport->Attribute("nom");
    if (nomTransport == "") throw ErreurXMLAttribut("nom", nomFichier, __FILE__, __LINE__);
    int e(0);
    for (e = 0; e < nombreTransports; e++) {
      if (nomTransport == nomTransports[e]) { break; }
    }
    if (e != nombreTransports) {
      erreur = elementTransport->QueryDoubleAttribute("valeur", &m_valeurTransport[e]);
      if (erreur != XML_NO_ERROR) throw ErreurXMLAttribut("valeur", nomFichier, __FILE__, __LINE__);
      couleurTrouvee++;
    }
    //Transport suivant
    elementTransport = elementTransport->NextSiblingElement("transport");
  }
  if (nombreTransports > couleurTrouvee) throw ErreurXMLAttribut("Pas assez d equations de tansport dans CL inj", nomFichier, __FILE__, __LINE__);
  m_nombreTransports = nombreTransports;

  //Allocation pour stocker les debits
  m_debits = new double[nombrePhases];
  m_nombrePhases = nombrePhases;
}

//****************************************************************************

CondLimSortie::CondLimSortie(double p0) : m_p0(p0){}

//****************************************************************************

CondLimSortie::CondLimSortie(const CondLimSortie& Source, const int lvl) : CondLim(Source)
{
  m_nombrePhases = Source.m_nombrePhases;
  m_nombreTransports = Source.m_nombreTransports;

  m_p0 = Source.m_p0;

  m_valeurTransport = new double[m_nombreTransports];
  for (int k = 0; k < m_nombreTransports; k++) {
    m_valeurTransport[k] = Source.m_valeurTransport[k];
  }
  
  m_debits = new double[m_nombrePhases];
  for (int k = 0; k < m_nombrePhases; k++) {
    m_debits[k] = 0.;
  }

  m_lvl = lvl;
}

//****************************************************************************

CondLimSortie::~CondLimSortie()
{
  delete[] m_valeurTransport;
  delete[] m_debits;
}

//****************************************************************************

void CondLimSortie::creeLimite(BordDeMaille **face)
{
  *face = new CondLimSortie(*(this));
}

//****************************************************************************

void CondLimSortie::resolRiemannLimite(Cellule &cellGauche, const int & nombrePhases, const double & dxGauche, double & dtMax)
{
  m_mod->resolRiemannSortie(cellGauche, nombrePhases, dxGauche, dtMax, m_p0, m_debits);
  for (int k = 0; k < nombrePhases; k++) {
    m_debits[k] *= this->getFace()->getSurface();
    //if (1) m_debits[k] *= 3.14*2.*this->getFace()->getPos().getY();
  }
}

//****************************************************************************

void CondLimSortie::resolRiemannTransportLimite(Cellule &cellGauche, const int & nombreTransports) const
{
	m_mod->resolRiemannTransportSortie(cellGauche, nombreTransports, m_valeurTransport);
}

//****************************************************************************

void CondLimSortie::afficheInfos()
{
  cout << m_numPhysique << endl;
  cout << m_p0 << endl;
}

//****************************************************************************

//double CondLimSortie::getDebit(int numPhase) const 
//{
//  return m_debits[numPhase];
//}

//****************************************************************************
//******************************Methode AMR***********************************
//****************************************************************************

void CondLimSortie::creerBordEnfant()
{
  m_bordsEnfants.push_back(new CondLimSortie(*this, m_lvl + 1));
}

//****************************************************************************