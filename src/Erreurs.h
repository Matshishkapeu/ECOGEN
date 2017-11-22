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

#ifndef ERREURS_H
#define ERREURS_H

//Definitions de classes d'Erreurs tout type
//Erreur de base
//Exceptions sur lecture XML

#include <iostream>
#include <sstream>
#include <cassert>
#include <vector>

// Macro permettant de faire un assert avec un message d'erreur.
# ifndef NDEBUG
#   define assertM(condition, message) \
  do {\
    if (!(condition)) {  \
    std::cerr << "---------------------------------------------------------" \
      << std::endl << "Erreur assertion non verifie" << std::endl << \
      "  fichier : " << __FILE__ << std::endl << \
      "  ligne : " << __LINE__ << std::endl << \
      "  assertion : `" #condition "` a echouee" << std::endl \
      << "  " << message << std::endl \
      << "---------------------------------------------------------" << std::endl; \
      std::exit(EXIT_FAILURE); \
    } \
  } while (false)
#else
#   define ASSERT(condition, message) do { } while (false)
#endif

class Erreurs
{
public:
  Erreurs();
  Erreurs(const std::string &message, const char* fichierSource = "non renseigne", int numeroLigne = -1);
  virtual ~Erreurs();

  static void messageErreur(const std::string &message);
  static void messageErreur(const std::string &message, double valeur);

  void setErreur(const std::string &message, const char* fichierSource = "", int numeroLigne = -1);
  void setErreur(const std::string &message, const double valeur);
  void afficheErreur();
  void ecritErreurFichier();
  static void arretCodeApresErreur(std::vector<Erreurs> &erreurs);

  //Accesseur
  int getEtat();

private:
  std::string m_message;
  int m_etat;
  std::string m_fichier;
  int m_ligne;
  double m_valeur; //!< permet de faire remonter une information en plus
};

extern std::vector<Erreurs> erreurs;

//Gestion des exceptions sur erreur code ECOGEN
//---------------------------------------------
class ErreurECOGEN : public std::exception
{
public:
  //***************
  ErreurECOGEN(std::string infoErreur = "", const char* fichierSource = "", int numeroLigne = -1):
    std::exception(), m_infoErreur(infoErreur), m_numeroLigne(numeroLigne), m_fichierSource(fichierSource) {}
  virtual ~ErreurECOGEN() throw() {}
  //***************
  virtual const char *what(void) const throw()
  {
    return "Exception ECOGEN : corriger et relancer le run";
  }
  //***************
  std::string infoErreur(void) const throw()
  {
    std::stringstream message;
    message << "--------------------------------------------------" << std::endl;
    message << this->what() << std::endl;
    message << "****************************************" << std::endl;
    if (m_fichierSource != "")
    {
      message << " infos sur exception code source :" << std::endl;
      message << "  fichier : '" << m_fichierSource << "'" << std::endl;
      if (m_numeroLigne != -1) message << "  ligne : " << m_numeroLigne << std::endl;
    }
    message << this->infosAdditionelles() << std::endl;
    message << "--------------------------------------------------" << std::endl;
    return message.str();
  }
  //***************
  virtual std::string infosAdditionelles(void) const throw()
  {
    return m_infoErreur;
  }
  //***************

private:
  int m_numeroLigne;
  std::string m_fichierSource;
  std::string m_infoErreur;
};


//Gestion des exceptions sur fichiers entrees XML
//---------------------------------------------------------------

class ErreurXML : public ErreurECOGEN
{
public:
  //***************
  ErreurXML(std::string fichierXML = "", const char* fichierSource = "", int numeroLigne = -1) :
    ErreurECOGEN(), m_fichierXML(fichierXML), m_fichierSource(fichierSource), m_numeroLigne(numeroLigne){}
  virtual ~ErreurXML() throw(){}
  //***************
  virtual const char *what(void) const throw()
  {
    return "Exception sur lecture fichier XML : fichier introuvable ou structure incorrecte";
  }
  //***************
  std::string infoErreur(void) const throw()
  {
    std::stringstream message;
    message << "--------------------------------------------------" << std::endl;
    message << this->what() << std::endl;
    message << "****************************************" << std::endl;
    if (m_fichierXML != "") { message << " fichier XML concerne : '" << m_fichierXML << "'" << std::endl; }
    if (m_fichierSource != "")
    {
      message << " infos sur exception code source :" << std::endl;
      message << "  fichier : '" << m_fichierSource << "'" << std::endl;
      if (m_numeroLigne != -1) message << "  ligne : " << m_numeroLigne << std::endl;
    }
    message << this->infosAdditionelles();
    message << "--------------------------------------------------" << std::endl;
    return message.str();
  }
  //***************
  virtual std::string infosAdditionelles(void) const throw()
  {
    return "";
  }
  //***************
  int numeroLigne() const throw () { return m_numeroLigne; }
  //***************
  std::string nomFichier() const throw () { return m_fichierSource; }
  //***************
private:
  int m_numeroLigne;
  std::string m_fichierSource;
  std::string m_fichierXML;
};

//---------------------------------------------------------------

class ErreurXMLRacine : public ErreurXML
{
public:
  //***************
  ErreurXMLRacine(std::string racine = "", std::string fichierXML = "", const char* fichierSource = "", int numeroLigne = -1) :
    ErreurXML(fichierXML, fichierSource, numeroLigne), m_racine(racine){}
  virtual ~ErreurXMLRacine() throw(){}
  //***************
  virtual const char *what(void) const throw()
  {
    return "Exception sur fichier XML : racine introuvable";
  }
  //***************
  virtual std::string infosAdditionelles(void) const throw()
  {
    std::stringstream message;
    message << " nom de la racine recherchee : '" << m_racine << "'" << std::endl;
    return message.str();
  }
  //***************
  std::string racine() const throw () { return m_racine; }
  //***************
private:
  std::string m_racine;
};

//---------------------------------------------------------------

class ErreurXMLElement : public ErreurXML
{
public:
  //***************
  ErreurXMLElement(std::string element = "", std::string fichierXML = "", const char* fichierSource = "", int numeroLigne = -1) :
    ErreurXML(fichierXML, fichierSource, numeroLigne), m_element(element){}
  virtual ~ErreurXMLElement() throw(){}
  //***************
  virtual const char *what(void) const throw()
  {
    return "Exception sur fichier XML : element introuvable ";
  }
  //***************
  virtual std::string infosAdditionelles(void) const throw()
  {
    std::stringstream message;
    message << " nom de l element recherche : '" << m_element << "'" << std::endl;
    return message.str();
  }
  //***************
  std::string element() const throw () { return m_element; }
  //***************
private:
  std::string m_element;
};

//---------------------------------------------------------------

class ErreurXMLAttribut : public ErreurXML
{
public:
  //***************
  ErreurXMLAttribut(std::string attribut = "", std::string fichierXML = "", const char* fichierSource = "", int numeroLigne = -1) :
    ErreurXML(fichierXML, fichierSource, numeroLigne), m_attribut(attribut){}
  virtual ~ErreurXMLAttribut() throw(){}
  //***************
  virtual const char *what(void) const throw()
  {
    return "Exception sur fichier XML : attribut introuvable ";
  }
  //***************
  virtual std::string infosAdditionelles(void) const throw()
  {
    std::stringstream message;
    message << " nom de l attribut recherche : '" << m_attribut << "'" << std::endl;
    return message.str();
  }
  //***************
  std::string attribut() const throw () { return m_attribut; }
  //***************
private:
  std::string m_attribut;
};

//---------------------------------------------------------------

class ErreurXMLDev : public ErreurXML
{
public:
  //***************
  ErreurXMLDev(std::string fichierXML = "", const char* fichierSource = "", int numeroLigne = -1) :
    ErreurXML(fichierXML, fichierSource, numeroLigne){}
  virtual ~ErreurXMLDev() throw(){}
  //***************
  virtual const char *what(void) const throw()
  {
    return "Exception sur fichier XML : morceaux de code en cours de developpement ";
  }
  //***************
private:
};

//---------------------------------------------------------------

class ErreurXMLLimite : public ErreurXML
{
public:
  //***************
  ErreurXMLLimite(std::string fichierXML = "", const char* fichierSource = "", int numeroLigne = -1) :
    ErreurXML(fichierXML, fichierSource, numeroLigne){}
  virtual ~ErreurXMLLimite() throw(){}
  //***************
  virtual const char *what(void) const throw()
  {
    return "Exception sur fichier XML : erreur sur conditions aux limites ";
  }
  //***************
private:
};

//---------------------------------------------------------------

class ErreurXMLTermeSource : public ErreurXML
{
public:
  //***************
  ErreurXMLTermeSource(std::string fichierXML = "", const char* fichierSource = "", int numeroLigne = -1) :
    ErreurXML(fichierXML, fichierSource, numeroLigne){}
  virtual ~ErreurXMLTermeSource() throw(){}
  //***************
  virtual const char *what(void) const throw()
  {
    return "Exception sur fichier XML : erreur sur choix des termes sources ";
  }
  //***************
private:
};

//---------------------------------------------------------------

class ErreurXMLEOS : public ErreurXML
{
public:
  //***************
  ErreurXMLEOS(std::string fichierXML = "", const char* fichierSource = "", int numeroLigne = -1) :
    ErreurXML(fichierXML, fichierSource, numeroLigne){}
  virtual ~ErreurXMLEOS() throw(){}
  //***************
  virtual const char *what(void) const throw()
  {
    return "Exception sur fichier XML : erreur sur equation etat ";
  }
  //***************
private:
};

//---------------------------------------------------------------

class ErreurXMLEOSInconnue : public ErreurXML
{
public:
  //***************
  ErreurXMLEOSInconnue(std::string typeEOS = "", std::string fichierXML = "", const char* fichierSource = "", int numeroLigne = -1) :
    ErreurXML(fichierXML, fichierSource, numeroLigne), m_typeEOS(typeEOS){}
  virtual ~ErreurXMLEOSInconnue() throw(){}
  //***************
  virtual const char *what(void) const throw()
  {
    return "Exception sur fichier XML : erreur sur equation etat ";
  }
  //***************
  virtual std::string infosAdditionelles(void) const throw()
  {
    std::stringstream message;
    message << " type EOS : '" << m_typeEOS << "' inconnu" << std::endl;
    return message.str();
  }
  //***************
private:
  std::string m_typeEOS;
};

//---------------------------------------------------------------

class ErreurXMLDomaineInconnu : public ErreurXML
{
public:
  //***************
  ErreurXMLDomaineInconnu(std::string typeDomaine = "", std::string fichierXML = "", const char* fichierSource = "", int numeroLigne = -1) :
    ErreurXML(fichierXML, fichierSource, numeroLigne), m_typeDomaine(typeDomaine){}
  virtual ~ErreurXMLDomaineInconnu() throw(){}
  //***************
  virtual const char *what(void) const throw()
  {
    return "Exception sur fichier XML : erreur sur domaine CI ";
  }
  //***************
  virtual std::string infosAdditionelles(void) const throw()
  {
    std::stringstream message;
    message << " type Domaine : '" << m_typeDomaine << "' inconnu" << std::endl;
    return message.str();
  }
  //***************
private:
  std::string m_typeDomaine;
};

//---------------------------------------------------------------

class ErreurXMLCondLimInconnue : public ErreurXML
{
public:
  //***************
  ErreurXMLCondLimInconnue(std::string typeCondLim = "", std::string fichierXML = "", const char* fichierSource = "", int numeroLigne = -1) :
    ErreurXML(fichierXML, fichierSource, numeroLigne), m_typeCondLim(typeCondLim) {}
  virtual ~ErreurXMLCondLimInconnue() throw() {}
  //***************
  virtual const char *what(void) const throw()
  {
    return "Exception sur fichier XML : erreur sur Condition Limite CL ";
  }
  //***************
  virtual std::string infosAdditionelles(void) const throw()
  {
    std::stringstream message;
    message << " type Limite : '" << m_typeCondLim << "' inconnu" << std::endl;
    return message.str();
  }
  //***************
private:
  std::string m_typeCondLim;
};

//---------------------------------------------------------------

class ErreurXMLEtat : public ErreurXML
{
public:
  //***************
  ErreurXMLEtat(std::string nomEtat = "", std::string fichierXML = "", const char* fichierSource = "", int numeroLigne = -1) :
    ErreurXML(fichierXML, fichierSource, numeroLigne), m_nomEtat(nomEtat){}
  virtual ~ErreurXMLEtat() throw(){}
  //***************
  virtual const char *what(void) const throw()
  {
    return "Exception sur fichier XML : erreur sur etat ";
  }
  //***************
  virtual std::string infosAdditionelles(void) const throw()
  {
    std::stringstream message;
    message << " Etat : '" << m_nomEtat << "' non trouve ou incomplet" << std::endl;
    return message.str();
  }
  //***************
private:
  std::string m_nomEtat;
};

//---------------------------------------------------------------

class ErreurXMLMateriauInconnu : public ErreurXML
{
public:
  //***************
  ErreurXMLMateriauInconnu(std::string nomMateriau = "", std::string fichierXML = "", const char* fichierSource = "", int numeroLigne = -1) :
    ErreurXML(fichierXML, fichierSource, numeroLigne), m_nomMateriau(nomMateriau){}
  virtual ~ErreurXMLMateriauInconnu() throw(){}
  //***************
  virtual const char *what(void) const throw()
  {
    return "Exception sur fichier XML : erreur sur etat CI ";
  }
  //***************
  virtual std::string infosAdditionelles(void) const throw()
  {
    std::stringstream message;
    message << " type Materiau : '" << m_nomMateriau << "' inconnu" << std::endl;
    return message.str();
  }
  //***************
private:
  std::string m_nomMateriau;
};

//---------------------------------------------------------------

#endif // ERREURS_H 