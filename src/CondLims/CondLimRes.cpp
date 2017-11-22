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

#include "CondLimRes.h"

using namespace std;
using namespace tinyxml2;

//****************************************************************************

CondLimRes::CondLimRes(){}

//****************************************************************************

CondLimRes::CondLimRes(int numPhysique, XMLElement *element, vector<Phase*> vecPhases, int &nombreTransports, std::vector<std::string> nomTransports, string nomFichier) :
  CondLim(numPhysique)
{
  m_nombrePhase = vecPhases.size();
  m_ak0 = new double[m_nombrePhase];
  m_rhok0 = new double[m_nombrePhase];
  m_pk0 = new double[m_nombrePhase];
  for (int k = 0; k < m_nombrePhase; k++)
  {
    m_ak0[k] = vecPhases[k]->getAlpha();
    m_rhok0[k] = vecPhases[k]->getDensite();
    m_pk0[k] = vecPhases[k]->getPression();
  }

  //Lecture des transports
  if (nombreTransports) {
    XMLElement *sousElement(element->FirstChildElement("donneesReservoir"));
    if (sousElement == NULL) throw ErreurXMLElement("donneesReservoir", nomFichier, __FILE__, __LINE__);
    XMLError erreur;

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
    if (nombreTransports > couleurTrouvee) throw ErreurXMLAttribut("phase", nomFichier, __FILE__, __LINE__);
  }
  m_nombreTransports = nombreTransports;
}

//****************************************************************************

CondLimRes::CondLimRes(double *ak0, double *rhok0, double *pk0, const int nombrePhases)
{
  m_ak0 = new double[nombrePhases];
  m_rhok0 = new double[nombrePhases];
  m_pk0 = new double[nombrePhases];

  for (int k = 0; k < nombrePhases; k++)
  {
    m_ak0[k] = ak0[k];
    m_rhok0[k] = rhok0[k];
    m_pk0[k] = pk0[k];
  }
}

//****************************************************************************

CondLimRes::CondLimRes(const CondLimRes &Source, const int lvl) : CondLim(Source)
{
  m_nombrePhase = Source.m_nombrePhase;
  m_nombreTransports = Source.m_nombreTransports;
  m_ak0 = new double[m_nombrePhase];
  m_rhok0 = new double[m_nombrePhase];
  m_pk0 = new double[m_nombrePhase];

  for (int k = 0; k < m_nombrePhase; k++)
  {
    m_ak0[k] = Source.m_ak0[k];
    m_rhok0[k] = Source.m_rhok0[k];
    m_pk0[k] = Source.m_pk0[k];
  }

  m_valeurTransport = new double[Source.m_nombreTransports];
  for (int k = 0; k < Source.m_nombreTransports; k++) {
    m_valeurTransport[k] = Source.m_valeurTransport[k];
  }

  m_lvl = lvl;
}

//****************************************************************************

CondLimRes::~CondLimRes()
{
  delete[] m_ak0;
  delete[] m_rhok0;
  delete[] m_pk0;
  delete[] m_valeurTransport;
}

//****************************************************************************

void CondLimRes::creeLimite(BordDeMaille **face)
{
  *face = new CondLimRes(*(this));
}

//****************************************************************************

void CondLimRes::resolRiemannLimite(Cellule &cellGauche, const int & nombrePhases, const double & dxGauche, double & dtMax)
{
  m_mod->resolRiemannRes(cellGauche, nombrePhases, dxGauche, dtMax, m_ak0, m_rhok0, m_pk0);
}

//****************************************************************************

void CondLimRes::resolRiemannTransportLimite(Cellule &cellGauche, const int & nombreTransports) const
{
	m_mod->resolRiemannTransportRes(cellGauche, nombreTransports, m_valeurTransport);
}

//****************************************************************************

void CondLimRes::afficheInfos()
{
  cout << m_numPhysique << endl;
  cout << m_rhok0[0] << endl;
}

//****************************************************************************
//******************************Methode AMR***********************************
//****************************************************************************

void CondLimRes::creerBordEnfant()
{
  m_bordsEnfants.push_back(new CondLimRes(*this, m_lvl + 1));
}

//****************************************************************************
