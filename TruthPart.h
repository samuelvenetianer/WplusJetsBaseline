#ifndef TruthPart_h
#define TruthPart_h


#include <vector>


class TruthPart
{  
 private :
 

  float _E;
  float _px;
  float _py;
  float _pz;
  float _m;
  float _pt;
  float _eta;
  float _phi;
  float _charge;
  int _index;
  int _status;
  int _statushepmc;
  int _pdgid;


 public :

  void Set_E ( const float & V_E) 	               {_E=V_E;};
  void Set_Px( const float & V_Px)                     {_px=V_Px;};
  void Set_Py( const float & V_Py)                     {_py=V_Py;};
  void Set_Pz( const float & V_Pz)                     {_pz=V_Pz;};
  void Set_M( const float & V_M)                       {_m=V_M;};
  void Set_Pt( const float & V_Pt)                     {_pt=V_Pt;};
  void Set_Eta( const float & V_Eta)                   {_eta=V_Eta;};
  void Set_Phi( const float & V_Phi)                   {_phi=V_Phi;};
  void Set_Charge( const float & V_Charge)             {_charge=V_Charge;};
  void Set_Index( const int & V_Index)                 {_index=V_Index;};
  void Set_Status( const int & V_Status)               {_status=V_Status;};
  void Set_StatusHepMC( const int & V_StatusHepMC)     {_statushepmc=V_StatusHepMC;};
  void Set_Pdgid( const int & V_Pdgid)                 {_pdgid=V_Pdgid;};
 

  float E() const           {return _E;};
  float Px() const          {return _px;};
  float Py() const          {return _py;};
  float Pz() const          {return _pz;};
  float M() const           {return _m;};
  float Pt() const          {return _pt;};
  float Eta() const         {return _eta;};
  float Phi() const         {return _phi;};
  float Charge() const      {return _charge ;};
  int Index() const         {return _index ;};
  int Status() const        {return _status ;};
  int StatusHepMC() const   {return _statushepmc ;};
  int Pdgid() const         {return _pdgid ;};


  TruthPart(const float &E=0,
	    const float &Px=0,
	    const float &Py=0,
	    const float &Pz=0,
	    const float &M=0,
	    const float &Pt=0,
	    const float &Eta=0,
	    const float &Phi=0,
	    const float &Charge =0,
	    const int &Index =0,
	    const int &Status =0,
	    const int &StatusHepMC =0,
	    const int &Pdgid =0
     )
  {  
    Set(E,
	Px,
	Py,
	Pz,
	M,
	Pt,
	Eta,
	Phi,
        Charge,
        Index, 
        Status, 
        StatusHepMC,
        Pdgid);
  }
  
  
  void Set(const float &E,
	   const float &Px,
	   const float &Py,
	   const float &Pz,
	   const float &M,
	   const float &Pt,
	   const float &Eta,
	   const float &Phi,
	   const float &Charge,
           const int &Index, 
           const int &Status, 
           const int &StatusHepMC,
           const int &Pdgid
	   )
  { _E=E;
    _px=Px;
    _py=Py;
    _pz=Pz;
    _m=M;
    _pt=Pt;
    _eta=Eta;
    _phi=Phi;
    _charge=Charge;
    _index=Index;
    _status=Status;
    _statushepmc=StatusHepMC;
    _pdgid=Pdgid;
  }	
};


#endif
