! -----------------------------------------------------------------------------  
! This file was automatically created by SARAH version 4.12.3 
! SARAH References: arXiv:0806.0538, 0909.2863, 1002.0840, 1207.0906, 1309.7223  
! (c) Florian Staub, 2013  
! ------------------------------------------------------------------------------  
! File created at 18:54 on 23.10.2021   
! ----------------------------------------------------------------------  
 
 
Module TreeLevelMasses_U1B 
 
Use Control 
Use Mathematics 
Use MathematicsQP 
Use Settings 
Use Model_Data_U1B 
!Use StandardModel 
 
 
Logical :: SignOfMassChanged =.False.  
Logical :: SignOfMuChanged =.False.  
Contains 
 
Subroutine TreeMasses(MAh,MAh2,MFd,MFd2,MFe,MFe2,MFeD,MFeD2,MFu,MFu2,MFv,             & 
& MFv2,MFvv,MFvv2,MFx,MFx2,Mhh,Mhh2,MHm,MHm2,Mss,Mss2,MVWm,MVWm2,MVZ,MVZ2,               & 
& MVZp,MVZp2,TW,TWp,ZDR,ZER,ZERD,ZUR,ZDL,ZEL,ZELD,ZUL,ZA,ZH,UV,ZP,ZW,ZZ,v,               & 
& vX,g1,gBY,g2,g3,g1p,gYB,lamh,lam2,lam3,lam4,lam5,lam6,lam7,Yu,Yd,Ye,YB1,               & 
& YC1,YA1,YB2,YC2,YA2,YRH,LE,YRG,LD,Fi,Fj,muh,mui,muj,m1,m2,GenerationMixing,kont)

Implicit None 
 
Real(dp),Intent(in) :: g1,gBY,g2,g3,g1p,gYB,YB1,YC1,YA1,YB2,YC2,YA2,YRH,Fi(2,2),Fj(2,2)

Complex(dp),Intent(in) :: lamh,lam2,lam3,lam4,lam5,lam6(2,2),lam7(2,2),Yu(3,3),Yd(3,3),Ye(3,3),LE(2,2),         & 
& YRG,LD(2,3),muh,mui,muj,m1(2,2),m2(2,2)

Real(dp),Intent(out) :: MAh(2),MAh2(2),MFd(3),MFd2(3),MFe(3),MFe2(3),MFeD(2),MFeD2(2),MFu(3),MFu2(3),         & 
& MFv(5),MFv2(5),MFvv,MFvv2,MFx,MFx2,Mhh(2),Mhh2(2),MHm(5),MHm2(5),Mss,Mss2,             & 
& MVWm,MVWm2,MVZ,MVZ2,MVZp,MVZp2,TW,TWp,ZA(2,2),ZH(2,2),ZP(5,5),ZZ(3,3)

Complex(dp),Intent(out) :: ZDR(3,3),ZER(3,3),ZERD(2,2),ZUR(3,3),ZDL(3,3),ZEL(3,3),ZELD(2,2),ZUL(3,3),            & 
& UV(5,5),ZW(2,2)

Real(dp),Intent(in) :: v,vX

Logical, Intent(in) :: GenerationMixing 
Integer, Intent(inout) :: kont 
Integer :: i1,i2,i3,i4,j1,j2,j3,kontSave 
Complex(dp) :: MFxC 
Complex(dp) :: MFvvC 
Iname = Iname + 1 
NameOfUnit(Iname) = 'TreeMassesU1B'
 
kont = 0 
Call CalculateVPVZVZp(g1,g2,g1p,gBY,gYB,v,vX,ZZ,MVZ,MVZp,MVZ2,MVZp2,TW,TWp,kont)

Call CalculateVWm(g2,v,ZW,MVWm,MVWm2,kont)

! ------------------------------- 
! Mass of Fx 
MFxC = (vX*YA1)/sqrt(2._dp) 
If (RotateNegativeFermionMasses) Then 
MFx = Abs(MFxC) 
MFx2 = MFx**2 
Else 
MFx = Real(MFxC,dp) 
MFx2 = MFx**2 
End if
! ------------------------------- 
! ------------------------------- 
! Mass of Fvv 
MFvvC = (vX*YC1)/sqrt(2._dp) 
If (RotateNegativeFermionMasses) Then 
MFvv = Abs(MFvvC) 
MFvv2 = MFvv**2 
Else 
MFvv = Real(MFvvC,dp) 
MFvv2 = MFvv**2 
End if
! ------------------------------- 
Mss2 = muj - (lam5*v**2)/2._dp 
  If (Mss2.ne.Mss2) Then 
      Write(*,*) 'NaN appearing in Mss2' 
      Call TerminateProgram 
    End If 
  If (Mss2.Ge.0._dp) Then 
Mss = sqrt(Mss2) 
  Else 
    If (ErrorLevel.Ge.0) Then 
      Write(10,*) 'Warning from Mss2' 
      Write(10,*) 'mass squarred is negative: ',Mss2 
    End If 
     Write(ErrCan,*) 'Warning from routine Mss2' 
     Write(ErrCan,*) 'in the calculation of the masses' 
     Write(ErrCan,*) 'occurred a negative mass squared!' 
     Write(ErrCan,*) Mss2 
     Write(*,*) 'Warning from routine Mss2' 
     Write(*,*) 'in the calculation of the masses' 
     Write(*,*) 'occurred a negative mass squared!' 
     Write(*,*) Mss2 
        Mss= 1._dp 
        Mss2= 1._dp 
   SignOfMassChanged = .True. 
 End if 


! ------------------------------- 
Call CalculateMhh(muh,lamh,mui,lam2,lam3,v,vX,ZH,Mhh,Mhh2,kont)

kontSave = kont 
Call CalculateMAh(g1,g2,g1p,gBY,gYB,muh,lamh,mui,lam2,lam3,v,vX,TW,TWp,               & 
& ZA,MAh,MAh2,kont)

kont = kontSave 
kontSave = kont 
Call CalculateMHm(g2,muh,lamh,lam3,m1,lam6,m2,lam7,Fi,v,vX,ZP,MHm,MHm2,kont)

kont = kontSave 
Call CalculateMFd(Yd,v,ZDL,ZDR,MFd,kont)

MFd2 = MFd**2 
Call CalculateMFu(Yu,v,ZUL,ZUR,MFu,kont)

MFu2 = MFu**2 
Call CalculateMFe(Ye,v,ZEL,ZER,MFe,kont)

MFe2 = MFe**2 
Call CalculateMFv(UV,MFv,kont)

MFv2 = MFv**2 
Call CalculateMFeD(YB1,YC1,YRH,YRG,v,vX,ZELD,ZERD,MFeD,kont)

MFeD2 = MFeD**2 

 
 Call SortGoldstones(MAh,MAh2,MFd,MFd2,MFe,MFe2,MFeD,MFeD2,MFu,MFu2,MFv,               & 
& MFv2,MFvv,MFvv2,MFx,MFx2,Mhh,Mhh2,MHm,MHm2,Mss,Mss2,MVWm,MVWm2,MVZ,MVZ2,               & 
& MVZp,MVZp2,TW,TWp,ZDR,ZER,ZERD,ZUR,ZDL,ZEL,ZELD,ZUL,ZA,ZH,UV,ZP,ZW,ZZ,kont)

If (SignOfMassChanged) Then  
 If (.Not.IgnoreNegativeMasses) Then 
  Write(*,*) " Stopping calculation because of negative mass squared." 
  Call TerminateProgram 
 Else 
  SignOfMassChanged= .False. 
  kont=0  
 End If 
End If 
If (SignOfMuChanged) Then 
 If (.Not.IgnoreMuSignFlip) Then 
  Write(*,*) " Stopping calculation because of negative mass squared in tadpoles." 
  Call TerminateProgram 
 Else 
  SignOfMuChanged= .False. 
  kont=0 
 End If 
End If 

 ! -------------------------------- 
! Setting Goldstone masses 
! -------------------------------- 
 
MAh(1)=MVZ*sqrt(RXiZ)
MAh2(1)=MVZ2*RXiZ
MAh(2)=MVZp*sqrt(RXiZp)
MAh2(2)=MVZp2*RXiZp
MHm(1)=MVWm*sqrt(RXiWm)
MHm2(1)=MVWm2*RXiWm
TW = ACos(Abs(ZZ(1,1)))
TWp = ACos(Abs(ZZ(3,3)))
Iname = Iname - 1 
 
End Subroutine  TreeMasses 
 
 
Subroutine RunningFermionMasses(MFeIN,MFe2IN,MFdIN,MFd2IN,MFuIN,MFu2IN,               & 
& v,vX,g1,gBY,g2,g3,g1p,gYB,lamh,lam2,lam3,lam4,lam5,lam6,lam7,Yu,Yd,Ye,YB1,             & 
& YC1,YA1,YB2,YC2,YA2,YRH,LE,YRG,LD,Fi,Fj,muh,mui,muj,m1,m2,kont)

Implicit None 
 
Integer, Intent(inout) :: kont 
Real(dp),Intent(in) :: g1,gBY,g2,g3,g1p,gYB,YB1,YC1,YA1,YB2,YC2,YA2,YRH,Fi(2,2),Fj(2,2)

Complex(dp),Intent(in) :: lamh,lam2,lam3,lam4,lam5,lam6(2,2),lam7(2,2),Yu(3,3),Yd(3,3),Ye(3,3),LE(2,2),         & 
& YRG,LD(2,3),muh,mui,muj,m1(2,2),m2(2,2)

Real(dp),Intent(in) :: v,vX

Real(dp),Intent(inout) :: MFeIN(3),MFe2IN(3),MFdIN(3),MFd2IN(3),MFuIN(3),MFu2IN(3)

Real(dp) :: MFe(3),MFe2(3),MFd(3),MFd2(3),MFu(3),MFu2(3)

Call CalculateMFe(Ye,v,ZEL,ZER,MFe,kont)

MFe2 = MFe**2 
MFeIN(1:2) = MFe(1:2) 
MFe2IN(1:2) = MFe2(1:2) 
Call CalculateMFd(Yd,v,ZDL,ZDR,MFd,kont)

MFd2 = MFd**2 
MFdIN(1:2) = MFd(1:2) 
MFd2IN(1:2) = MFd2(1:2) 
Call CalculateMFu(Yu,v,ZUL,ZUR,MFu,kont)

MFu2 = MFu**2 
MFuIN(1:2) = MFu(1:2) 
MFu2IN(1:2) = MFu2(1:2) 
End Subroutine RunningFermionMasses 

Subroutine TreeMassesEffPot(MAh,MAh2,MFd,MFd2,MFe,MFe2,MFeD,MFeD2,MFu,MFu2,           & 
& MFv,MFv2,MFvv,MFvv2,MFx,MFx2,Mhh,Mhh2,MHm,MHm2,Mss,Mss2,MVWm,MVWm2,MVZ,MVZ2,           & 
& MVZp,MVZp2,TW,TWp,ZDR,ZER,ZERD,ZUR,ZDL,ZEL,ZELD,ZUL,ZA,ZH,UV,ZP,ZW,ZZ,v,               & 
& vX,g1,gBY,g2,g3,g1p,gYB,lamh,lam2,lam3,lam4,lam5,lam6,lam7,Yu,Yd,Ye,YB1,               & 
& YC1,YA1,YB2,YC2,YA2,YRH,LE,YRG,LD,Fi,Fj,muh,mui,muj,m1,m2,GenerationMixing,kont)

Implicit None 
 
Real(dp),Intent(in) :: g1,gBY,g2,g3,g1p,gYB,YB1,YC1,YA1,YB2,YC2,YA2,YRH,Fi(2,2),Fj(2,2)

Complex(dp),Intent(in) :: lamh,lam2,lam3,lam4,lam5,lam6(2,2),lam7(2,2),Yu(3,3),Yd(3,3),Ye(3,3),LE(2,2),         & 
& YRG,LD(2,3),muh,mui,muj,m1(2,2),m2(2,2)

Real(dp),Intent(out) :: MAh(2),MAh2(2),MFd(3),MFd2(3),MFe(3),MFe2(3),MFeD(2),MFeD2(2),MFu(3),MFu2(3),         & 
& MFv(5),MFv2(5),MFvv,MFvv2,MFx,MFx2,Mhh(2),Mhh2(2),MHm(5),MHm2(5),Mss,Mss2,             & 
& MVWm,MVWm2,MVZ,MVZ2,MVZp,MVZp2,TW,TWp,ZA(2,2),ZH(2,2),ZP(5,5),ZZ(3,3)

Complex(dp),Intent(out) :: ZDR(3,3),ZER(3,3),ZERD(2,2),ZUR(3,3),ZDL(3,3),ZEL(3,3),ZELD(2,2),ZUL(3,3),            & 
& UV(5,5),ZW(2,2)

Real(dp),Intent(in) :: v,vX

Logical, Intent(in) :: GenerationMixing 
Integer, Intent(inout) :: kont 
Integer :: i1,i2,i3,i4,j1,j2,j3,kontSave 
Complex(dp) :: MFxC 
Complex(dp) :: MFvvC 
Iname = Iname + 1 
NameOfUnit(Iname) = 'TreeMassesU1B'
 
kont = 0 
Call CalculateVPVZVZpEffPot(g1,g2,g1p,gBY,gYB,v,vX,ZZ,MVZ,MVZp,MVZ2,MVZp2,            & 
& TW,TWp,kont)

Call CalculateVWmEffPot(g2,v,ZW,MVWm,MVWm2,kont)

! ------------------------------- 
! Mass of Fx 
MFxC = (vX*YA1)/sqrt(2._dp) 
If (RotateNegativeFermionMasses) Then 
MFx = Abs(MFxC) 
MFx2 = MFx**2 
Else 
MFx = Real(MFxC,dp) 
MFx2 = MFx**2 
End if
! ------------------------------- 
! ------------------------------- 
! Mass of Fvv 
MFvvC = (vX*YC1)/sqrt(2._dp) 
If (RotateNegativeFermionMasses) Then 
MFvv = Abs(MFvvC) 
MFvv2 = MFvv**2 
Else 
MFvv = Real(MFvvC,dp) 
MFvv2 = MFvv**2 
End if
! ------------------------------- 
Mss2 = muj - (lam5*v**2)/2._dp 
  If (Mss2.ne.Mss2) Then 
      Write(*,*) 'NaN appearing in Mss2' 
      Call TerminateProgram 
    End If 
  If (Mss2.Ge.0._dp) Then 
Mss = sqrt(Mss2) 
  Else 
 End if 


! ------------------------------- 
Call CalculateMhhEffPot(muh,lamh,mui,lam2,lam3,v,vX,ZH,Mhh,Mhh2,kont)

kontSave = kont 
Call CalculateMAhEffPot(g1,g2,g1p,gBY,gYB,muh,lamh,mui,lam2,lam3,v,vX,TW,             & 
& TWp,ZA,MAh,MAh2,kont)

kont = kontSave 
kontSave = kont 
Call CalculateMHmEffPot(g2,muh,lamh,lam3,m1,lam6,m2,lam7,Fi,v,vX,ZP,MHm,              & 
& MHm2,kont)

kont = kontSave 
Call CalculateMFdEffPot(Yd,v,ZDL,ZDR,MFd,kont)

MFd2 = MFd**2 
Call CalculateMFuEffPot(Yu,v,ZUL,ZUR,MFu,kont)

MFu2 = MFu**2 
Call CalculateMFeEffPot(Ye,v,ZEL,ZER,MFe,kont)

MFe2 = MFe**2 
Call CalculateMFvEffPot(UV,MFv,kont)

MFv2 = MFv**2 
Call CalculateMFeDEffPot(YB1,YC1,YRH,YRG,v,vX,ZELD,ZERD,MFeD,kont)

MFeD2 = MFeD**2 

 
 If (SignOfMassChanged) Then  
 If (.Not.IgnoreNegativeMasses) Then 
  Write(*,*) " Stopping calculation because of negative mass squared." 
  Call TerminateProgram 
 Else 
  SignOfMassChanged= .False. 
  kont=0  
 End If 
End If 
If (SignOfMuChanged) Then 
 If (.Not.IgnoreMuSignFlip) Then 
  Write(*,*) " Stopping calculation because of negative mass squared in tadpoles." 
  Call TerminateProgram 
 Else 
  SignOfMuChanged= .False. 
  kont=0 
 End If 
End If 
Iname = Iname - 1 
 
End Subroutine  TreeMassesEffPot 
 
 
Subroutine CalculateMhh(muh,lamh,mui,lam2,lam3,v,vX,ZH,Mhh,Mhh2,kont)

Real(dp), Intent(in) :: v,vX

Complex(dp), Intent(in) :: muh,lamh,mui,lam2,lam3

Integer, Intent(inout) :: kont 
Integer :: i1,i2,i3,i4, ierr 
Integer :: j1,j2,j3,j4, pos 
Real(dp), Intent(out) :: Mhh(2), Mhh2(2) 
Real(dp), Intent(out) :: ZH(2,2) 
 
Real(dp) :: mat(2,2)  

Real(dp) ::  test(2) 

Iname = Iname + 1 
NameOfUnit(Iname) = 'CalculateMhh'
 
mat(1,1) = 0._dp 
mat(1,1) = mat(1,1)+muh
mat(1,1) = mat(1,1)-3*lamh*v**2
mat(1,1) = mat(1,1)-(lam3*vX**2)/2._dp
mat(1,2) = 0._dp 
mat(1,2) = mat(1,2)-(lam3*v*vX)
mat(2,2) = 0._dp 
mat(2,2) = mat(2,2)+mui
mat(2,2) = mat(2,2)-(lam3*v**2)/2._dp
mat(2,2) = mat(2,2)-3*lam2*vX**2

 
 Do i1=2,2
  Do i2 = 1, i1-1 
  mat(i1,i2) = mat(i2,i1) 
  End do 
End do 

 
Call EigenSystem(mat,Mhh2,ZH,ierr,test) 
 
 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
  End If 
  ierr = 0 
End If 
 
If ((ierr.Ne.0.).And.(ErrorLevel.Ge.-1)) Then 
  Write(10,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
  Write(10,*) 'Diagonalization failed, ierr : ',ierr 
  kont = ierr 
  Iname = Iname - 1 
  Return 
End If 


Do i1=1,2
  If (Abs(Mhh2(i1)).Le.MaxMassNumericalZero) Mhh2(i1) = 1.E-10_dp 
  If (Mhh2(i1).ne.Mhh2(i1)) Then 
      Write(*,*) 'NaN appearing in '//NameOfUnit(Iname) 
      Call TerminateProgram 
    End If 
  If (Mhh2(i1).Ge.0._dp) Then 
  Mhh(i1)=Sqrt(Mhh2(i1) ) 
  Else 
    If (ErrorLevel.Ge.0) Then 
      Write(10,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
      Write(10,*) 'a mass squarred is negative: ',i1,Mhh2(i1) 
    End If 
  Mhh(i1) = 1._dp 
     Write(ErrCan,*) 'Warning from routine '//NameOfUnit(Iname) 
     Write(ErrCan,*) 'in the calculation of the masses' 
     Write(ErrCan,*) 'occurred a negative mass squared!' 
     Write(ErrCan,*) i1,Mhh2(i1) 
     Write(*,*) 'Warning from routine '//NameOfUnit(Iname) 
     Write(*,*) 'in the calculation of the masses' 
     Write(*,*) 'occurred a negative mass squared!' 
     Write(*,*) i1,Mhh2(i1) 
  Mhh2(i1) = 1._dp 
   SignOfMassChanged = .True. 
! kont = -104 
 End if 
End Do 
Iname = Iname - 1 
 
End Subroutine CalculateMhh 

Subroutine CalculateMAh(g1,g2,g1p,gBY,gYB,muh,lamh,mui,lam2,lam3,v,vX,TW,             & 
& TWp,ZA,MAh,MAh2,kont)

Real(dp), Intent(in) :: g1,g2,g1p,gBY,gYB,v,vX,TW,TWp

Complex(dp), Intent(in) :: muh,lamh,mui,lam2,lam3

Integer, Intent(inout) :: kont 
Integer :: i1,i2,i3,i4, ierr 
Integer :: j1,j2,j3,j4, pos 
Real(dp), Intent(out) :: MAh(2), MAh2(2) 
Real(dp), Intent(out) :: ZA(2,2) 
 
Real(dp) :: mat(2,2)  

Real(dp) ::  test(2) 

Iname = Iname + 1 
NameOfUnit(Iname) = 'CalculateMAh'
 
mat(1,1) = 0._dp 
mat(1,1) = mat(1,1)+muh
mat(1,1) = mat(1,1)-(lamh*v**2)
mat(1,1) = mat(1,1)-(lam3*vX**2)/2._dp
mat(1,1) = mat(1,1)+(g2**2*v**2*Cos(TW)**2*Cos(TWp)**2*RXiZ)/4._dp
mat(1,1) = mat(1,1)+(gBY**2*v**2*Cos(TWp)**2*RXiZp)/4._dp
mat(1,1) = mat(1,1)+(g1*g2*v**2*Cos(TW)*Cos(TWp)**2*RXiZ*Sin(TW))/2._dp
mat(1,1) = mat(1,1)+(g1**2*v**2*Cos(TWp)**2*RXiZ*Sin(TW)**2)/4._dp
mat(1,1) = mat(1,1)-(gBY*g2*v**2*Cos(TW)*Cos(TWp)*RXiZ*Sin(TWp))/2._dp
mat(1,1) = mat(1,1)+(gBY*g2*v**2*Cos(TW)*Cos(TWp)*RXiZp*Sin(TWp))/2._dp
mat(1,1) = mat(1,1)-(g1*gBY*v**2*Cos(TWp)*RXiZ*Sin(TW)*Sin(TWp))/2._dp
mat(1,1) = mat(1,1)+(g1*gBY*v**2*Cos(TWp)*RXiZp*Sin(TW)*Sin(TWp))/2._dp
mat(1,1) = mat(1,1)+(gBY**2*v**2*RXiZ*Sin(TWp)**2)/4._dp
mat(1,1) = mat(1,1)+(g2**2*v**2*Cos(TW)**2*RXiZp*Sin(TWp)**2)/4._dp
mat(1,1) = mat(1,1)+(g1*g2*v**2*Cos(TW)*RXiZp*Sin(TW)*Sin(TWp)**2)/2._dp
mat(1,1) = mat(1,1)+(g1**2*v**2*RXiZp*Sin(TW)**2*Sin(TWp)**2)/4._dp
mat(1,2) = 0._dp 
mat(1,2) = mat(1,2)+(5*gBY*g1p*v*vX*Cos(TWp)**2*RXiZp)/2._dp
mat(1,2) = mat(1,2)+(5*gYB*g2*v*vX*Cos(TW)*Cos(TWp)**2*RXiZ*Sin(TW))/2._dp
mat(1,2) = mat(1,2)+(5*g1*gYB*v*vX*Cos(TWp)**2*RXiZ*Sin(TW)**2)/2._dp
mat(1,2) = mat(1,2)+(-5*g1p*g2*v*vX*Cos(TW)*Cos(TWp)*RXiZ*Sin(TWp))/2._dp
mat(1,2) = mat(1,2)+(5*g1p*g2*v*vX*Cos(TW)*Cos(TWp)*RXiZp*Sin(TWp))/2._dp
mat(1,2) = mat(1,2)+(-5*g1*g1p*v*vX*Cos(TWp)*RXiZ*Sin(TW)*Sin(TWp))/2._dp
mat(1,2) = mat(1,2)+(-5*gBY*gYB*v*vX*Cos(TWp)*RXiZ*Sin(TW)*Sin(TWp))/2._dp
mat(1,2) = mat(1,2)+(5*g1*g1p*v*vX*Cos(TWp)*RXiZp*Sin(TW)*Sin(TWp))/2._dp
mat(1,2) = mat(1,2)+(5*gBY*gYB*v*vX*Cos(TWp)*RXiZp*Sin(TW)*Sin(TWp))/2._dp
mat(1,2) = mat(1,2)+(5*gBY*g1p*v*vX*RXiZ*Sin(TWp)**2)/2._dp
mat(1,2) = mat(1,2)+(5*gYB*g2*v*vX*Cos(TW)*RXiZp*Sin(TW)*Sin(TWp)**2)/2._dp
mat(1,2) = mat(1,2)+(5*g1*gYB*v*vX*RXiZp*Sin(TW)**2*Sin(TWp)**2)/2._dp
mat(2,2) = 0._dp 
mat(2,2) = mat(2,2)+mui
mat(2,2) = mat(2,2)-(lam3*v**2)/2._dp
mat(2,2) = mat(2,2)-(lam2*vX**2)
mat(2,2) = mat(2,2)+25*g1p**2*vX**2*Cos(TWp)**2*RXiZp
mat(2,2) = mat(2,2)+25*gYB**2*vX**2*Cos(TWp)**2*RXiZ*Sin(TW)**2
mat(2,2) = mat(2,2)-50*g1p*gYB*vX**2*Cos(TWp)*RXiZ*Sin(TW)*Sin(TWp)
mat(2,2) = mat(2,2)+50*g1p*gYB*vX**2*Cos(TWp)*RXiZp*Sin(TW)*Sin(TWp)
mat(2,2) = mat(2,2)+25*g1p**2*vX**2*RXiZ*Sin(TWp)**2
mat(2,2) = mat(2,2)+25*gYB**2*vX**2*RXiZp*Sin(TW)**2*Sin(TWp)**2

 
 Do i1=2,2
  Do i2 = 1, i1-1 
  mat(i1,i2) = mat(i2,i1) 
  End do 
End do 

 
Call EigenSystem(mat,MAh2,ZA,ierr,test) 
 
 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
  End If 
  ierr = 0 
End If 
 
If ((ierr.Ne.0.).And.(ErrorLevel.Ge.-1)) Then 
  Write(10,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
  Write(10,*) 'Diagonalization failed, ierr : ',ierr 
  kont = ierr 
  Iname = Iname - 1 
  Return 
End If 


Do i1=1,2
  If (Abs(MAh2(i1)).Le.MaxMassNumericalZero) MAh2(i1) = 1.E-10_dp 
  If (MAh2(i1).ne.MAh2(i1)) Then 
      Write(*,*) 'NaN appearing in '//NameOfUnit(Iname) 
      Call TerminateProgram 
    End If 
  If (MAh2(i1).Ge.0._dp) Then 
  MAh(i1)=Sqrt(MAh2(i1) ) 
  Else 
    If (ErrorLevel.Ge.0) Then 
      Write(10,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
      Write(10,*) 'a mass squarred is negative: ',i1,MAh2(i1) 
    End If 
  MAh(i1) = 1._dp 
     Write(ErrCan,*) 'Warning from routine '//NameOfUnit(Iname) 
     Write(ErrCan,*) 'in the calculation of the masses' 
     Write(ErrCan,*) 'occurred a negative mass squared!' 
     Write(ErrCan,*) i1,MAh2(i1) 
     Write(*,*) 'Warning from routine '//NameOfUnit(Iname) 
     Write(*,*) 'in the calculation of the masses' 
     Write(*,*) 'occurred a negative mass squared!' 
     Write(*,*) i1,MAh2(i1) 
  MAh2(i1) = 1._dp 
   SignOfMassChanged = .True. 
! kont = -104 
 End if 
End Do 
Iname = Iname - 1 
 
End Subroutine CalculateMAh 

Subroutine CalculateMHm(g2,muh,lamh,lam3,m1,lam6,m2,lam7,Fi,v,vX,ZP,MHm,              & 
& MHm2,kont)

Real(dp), Intent(in) :: g2,Fi(2,2),v,vX

Complex(dp), Intent(in) :: muh,lamh,lam3,m1(2,2),lam6(2,2),m2(2,2),lam7(2,2)

Integer, Intent(inout) :: kont 
Integer :: i1,i2,i3,i4, ierr 
Integer :: j1,j2,j3,j4, pos 
Real(dp), Intent(out) :: MHm(5), MHm2(5) 
Real(dp), Intent(out) :: ZP(5,5) 
 
Real(dp) :: mat(5,5)  

Real(dp) ::  test(2) 

Iname = Iname + 1 
NameOfUnit(Iname) = 'CalculateMHm'
 
mat(1,1) = 0._dp 
mat(1,1) = mat(1,1)+muh
mat(1,1) = mat(1,1)-(lamh*v**2)
mat(1,1) = mat(1,1)-(lam3*vX**2)/2._dp
mat(1,1) = mat(1,1)+(g2**2*v**2*RXiWm)/4._dp
mat(1,2) = 0._dp 
mat(1,3) = 0._dp 
mat(1,4) = 0._dp 
mat(1,5) = 0._dp 
mat(2,2) = 0._dp 
mat(2,2) = mat(2,2)+(v**2*lam6(1,1))/2._dp
mat(2,2) = mat(2,2)+m1(1,1)
mat(2,3) = 0._dp 
mat(2,3) = mat(2,3)+(v**2*lam6(1,2))/2._dp
mat(2,3) = mat(2,3)+m1(1,2)
mat(2,4) = 0._dp 
mat(2,4) = mat(2,4)+(vX*Fi(1,1))/sqrt(2._dp)
mat(2,5) = 0._dp 
mat(2,5) = mat(2,5)+(vX*Fi(1,2))/sqrt(2._dp)
mat(3,3) = 0._dp 
mat(3,3) = mat(3,3)+(v**2*lam6(2,2))/2._dp
mat(3,3) = mat(3,3)+m1(2,2)
mat(3,4) = 0._dp 
mat(3,4) = mat(3,4)+(vX*Fi(2,1))/sqrt(2._dp)
mat(3,5) = 0._dp 
mat(3,5) = mat(3,5)+(vX*Fi(2,2))/sqrt(2._dp)
mat(4,4) = 0._dp 
mat(4,4) = mat(4,4)+(v**2*lam7(1,1))/2._dp
mat(4,4) = mat(4,4)+m2(1,1)
mat(4,5) = 0._dp 
mat(4,5) = mat(4,5)+(v**2*lam7(1,2))/2._dp
mat(4,5) = mat(4,5)+m2(1,2)
mat(5,5) = 0._dp 
mat(5,5) = mat(5,5)+(v**2*lam7(2,2))/2._dp
mat(5,5) = mat(5,5)+m2(2,2)

 
 Do i1=2,5
  Do i2 = 1, i1-1 
  mat(i1,i2) = mat(i2,i1) 
  End do 
End do 

 
Call EigenSystem(mat,MHm2,ZP,ierr,test) 
 
 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
  End If 
  ierr = 0 
End If 
 
If ((ierr.Ne.0.).And.(ErrorLevel.Ge.-1)) Then 
  Write(10,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
  Write(10,*) 'Diagonalization failed, ierr : ',ierr 
  kont = ierr 
  Iname = Iname - 1 
  Return 
End If 


Do i1=1,5
  If (Abs(MHm2(i1)).Le.MaxMassNumericalZero) MHm2(i1) = 1.E-10_dp 
  If (MHm2(i1).ne.MHm2(i1)) Then 
      Write(*,*) 'NaN appearing in '//NameOfUnit(Iname) 
      Call TerminateProgram 
    End If 
  If (MHm2(i1).Ge.0._dp) Then 
  MHm(i1)=Sqrt(MHm2(i1) ) 
  Else 
    If (ErrorLevel.Ge.0) Then 
      Write(10,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
      Write(10,*) 'a mass squarred is negative: ',i1,MHm2(i1) 
    End If 
  MHm(i1) = 1._dp 
     Write(ErrCan,*) 'Warning from routine '//NameOfUnit(Iname) 
     Write(ErrCan,*) 'in the calculation of the masses' 
     Write(ErrCan,*) 'occurred a negative mass squared!' 
     Write(ErrCan,*) i1,MHm2(i1) 
     Write(*,*) 'Warning from routine '//NameOfUnit(Iname) 
     Write(*,*) 'in the calculation of the masses' 
     Write(*,*) 'occurred a negative mass squared!' 
     Write(*,*) i1,MHm2(i1) 
  MHm2(i1) = 1._dp 
   SignOfMassChanged = .True. 
! kont = -104 
 End if 
End Do 
Iname = Iname - 1 
 
End Subroutine CalculateMHm 

Subroutine CalculateMFd(Yd,v,ZDL,ZDR,MFd,kont)

Real(dp),Intent(in) :: v

Complex(dp),Intent(in) :: Yd(3,3)

Integer, Intent(inout) :: kont 
Integer :: i1,i2,i3,i4, ierr 
Integer :: j1,j2,j3,j4 
Real(dp), Intent(out) :: MFd(3) 
 Complex(dp), Intent(out) :: ZDL(3,3), ZDR(3,3) 
 
 Complex(dp) :: mat(3,3)=0._dp, mat2(3,3)=0._dp, phaseM 

Complex(dp) :: ZDL2(3,3), ZDR2(3,3) 
 
 Real(dp) :: ZDL1(3,3), ZDR1(3,3), test(2), MFd2(3) 
 
 Iname = Iname + 1 
NameOfUnit(Iname) = 'CalculateMFd'
 
MFd = 0._dp 
ZDL = 0._dp 
ZDR = 0._dp 
mat(1,1) = 0._dp 
mat(1,1) = mat(1,1)+(v*Yd(1,1))/sqrt(2._dp)
mat(1,2) = 0._dp 
mat(1,2) = mat(1,2)+(v*Yd(2,1))/sqrt(2._dp)
mat(1,3) = 0._dp 
mat(1,3) = mat(1,3)+(v*Yd(3,1))/sqrt(2._dp)
mat(2,1) = 0._dp 
mat(2,1) = mat(2,1)+(v*Yd(1,2))/sqrt(2._dp)
mat(2,2) = 0._dp 
mat(2,2) = mat(2,2)+(v*Yd(2,2))/sqrt(2._dp)
mat(2,3) = 0._dp 
mat(2,3) = mat(2,3)+(v*Yd(3,2))/sqrt(2._dp)
mat(3,1) = 0._dp 
mat(3,1) = mat(3,1)+(v*Yd(1,3))/sqrt(2._dp)
mat(3,2) = 0._dp 
mat(3,2) = mat(3,2)+(v*Yd(2,3))/sqrt(2._dp)
mat(3,3) = 0._dp 
mat(3,3) = mat(3,3)+(v*Yd(3,3))/sqrt(2._dp)

 
mat2 = Matmul(Transpose(Conjg(mat)),mat) 
If (Maxval(Abs(Aimag(mat2))).Eq.0._dp) Then 
Call EigenSystem(Real(mat2,dp),MFd2,ZDR1,ierr,test) 
ZDR2 = ZDR1 
Else 
Call EigenSystem(mat2,MFd2,ZDR2,ierr,test) 
 End If 
 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
  End If 
  ierr = 0 
End If 
 
mat2 = Matmul(mat,Transpose(Conjg(mat))) 
If (Maxval(Abs(Aimag(mat2))).Eq.0._dp) Then 
Call EigenSystem (Real(mat2,dp),MFd2,ZDL1,ierr,test) 
                  
                  
ZDL2 = ZDL1 
Else 
Call EigenSystem(mat2,MFd2,ZDL2,ierr,test) 
 
 
End If 
ZDL2 = Conjg(ZDL2) 
 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
  End If 
  ierr = 0 
End If 
 
mat2 = Matmul(Matmul( Conjg(ZDL2),mat),Transpose( Conjg(ZDR2))) 
Do i1=1,3
If (Abs(mat2(i1,i1)).gt.0._dp) Then 
phaseM = mat2(i1,i1) / Abs(mat2(i1,i1)) 
ZDR2(i1,:) = phaseM *ZDR2(i1,:) 
 End if 
End Do 
 
Do i1=1,3
If (Abs(ZDR2(i1,i1)).gt.0._dp) Then 
phaseM = ZDR2(i1,i1) / Abs(ZDR2(i1,i1)) 
ZDR2(i1,:) = Conjg(phaseM) *ZDR2(i1,:) 
 ZDL2(i1,:) = phaseM *ZDL2(i1,:) 
 End if 
  If (MFd2(i1).ne.MFd2(i1)) Then 
      Write(*,*) 'NaN appearing in '//NameOfUnit(Iname) 
      Call TerminateProgram 
    End If 
  If (Abs(MFd2(i1)).Le.MaxMassNumericalZero) MFd2(i1) = Abs(MFd2(i1))+1.E-10_dp 
  If (MFd2(i1).Le.0._dp) Then 
    If (ErrorLevel.Ge.0) Then 
      Write(10,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
      Write(10,*) 'a mass squarred is negative: ',i1,MFd2(i1) 
      Write(*,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
      Write(*,*) 'a mass squarred is negative: ',i1,MFd2(i1) 
      Call TerminateProgram 
    End If 
     Write(ErrCan,*) 'Warning from routine '//NameOfUnit(Iname) 
     Write(ErrCan,*) 'in the calculation of the masses' 
     Write(ErrCan,*) 'occurred a negative mass squared!' 
     Write(ErrCan,*) i1,MFd2(i1) 
     Write(*,*) 'Warning from routine '//NameOfUnit(Iname) 
     Write(*,*) 'in the calculation of the masses' 
     Write(*,*) 'occurred a negative mass squared!' 
     Write(*,*) i1,MFd2(i1) 
  MFd2(i1) = 1._dp 
   SignOfMassChanged = .True. 
! kont = -104 
 End if 
End Do 
 
If (ierr.Ne.0.) Then 
  Write(ErrCan,*) 'Warning from Subroutine CalculateMFd, ierr =',ierr 
  kont = ierr 
  Iname = Iname - 1 
  Return 
End If 


MFd = Sqrt( MFd2 ) 
ZDL = ZDL2 
ZDR = ZDR2 
Iname = Iname - 1 
 
End Subroutine CalculateMFd 

Subroutine CalculateMFu(Yu,v,ZUL,ZUR,MFu,kont)

Real(dp),Intent(in) :: v

Complex(dp),Intent(in) :: Yu(3,3)

Integer, Intent(inout) :: kont 
Integer :: i1,i2,i3,i4, ierr 
Integer :: j1,j2,j3,j4 
Real(dp), Intent(out) :: MFu(3) 
 Complex(dp), Intent(out) :: ZUL(3,3), ZUR(3,3) 
 
 Complex(dp) :: mat(3,3)=0._dp, mat2(3,3)=0._dp, phaseM 

Complex(dp) :: ZUL2(3,3), ZUR2(3,3) 
 
 Real(dp) :: ZUL1(3,3), ZUR1(3,3), test(2), MFu2(3) 
 
 Iname = Iname + 1 
NameOfUnit(Iname) = 'CalculateMFu'
 
MFu = 0._dp 
ZUL = 0._dp 
ZUR = 0._dp 
mat(1,1) = 0._dp 
mat(1,1) = mat(1,1)-((v*Yu(1,1))/sqrt(2._dp))
mat(1,2) = 0._dp 
mat(1,2) = mat(1,2)-((v*Yu(2,1))/sqrt(2._dp))
mat(1,3) = 0._dp 
mat(1,3) = mat(1,3)-((v*Yu(3,1))/sqrt(2._dp))
mat(2,1) = 0._dp 
mat(2,1) = mat(2,1)-((v*Yu(1,2))/sqrt(2._dp))
mat(2,2) = 0._dp 
mat(2,2) = mat(2,2)-((v*Yu(2,2))/sqrt(2._dp))
mat(2,3) = 0._dp 
mat(2,3) = mat(2,3)-((v*Yu(3,2))/sqrt(2._dp))
mat(3,1) = 0._dp 
mat(3,1) = mat(3,1)-((v*Yu(1,3))/sqrt(2._dp))
mat(3,2) = 0._dp 
mat(3,2) = mat(3,2)-((v*Yu(2,3))/sqrt(2._dp))
mat(3,3) = 0._dp 
mat(3,3) = mat(3,3)-((v*Yu(3,3))/sqrt(2._dp))

 
mat2 = Matmul(Transpose(Conjg(mat)),mat) 
If (Maxval(Abs(Aimag(mat2))).Eq.0._dp) Then 
Call EigenSystem(Real(mat2,dp),MFu2,ZUR1,ierr,test) 
ZUR2 = ZUR1 
Else 
Call EigenSystem(mat2,MFu2,ZUR2,ierr,test) 
 End If 
 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
  End If 
  ierr = 0 
End If 
 
mat2 = Matmul(mat,Transpose(Conjg(mat))) 
If (Maxval(Abs(Aimag(mat2))).Eq.0._dp) Then 
Call EigenSystem (Real(mat2,dp),MFu2,ZUL1,ierr,test) 
                  
                  
ZUL2 = ZUL1 
Else 
Call EigenSystem(mat2,MFu2,ZUL2,ierr,test) 
 
 
End If 
ZUL2 = Conjg(ZUL2) 
 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
  End If 
  ierr = 0 
End If 
 
mat2 = Matmul(Matmul( Conjg(ZUL2),mat),Transpose( Conjg(ZUR2))) 
Do i1=1,3
If (Abs(mat2(i1,i1)).gt.0._dp) Then 
phaseM = mat2(i1,i1) / Abs(mat2(i1,i1)) 
ZUR2(i1,:) = phaseM *ZUR2(i1,:) 
 End if 
End Do 
 
Do i1=1,3
If (Abs(ZUR2(i1,i1)).gt.0._dp) Then 
phaseM = ZUR2(i1,i1) / Abs(ZUR2(i1,i1)) 
ZUR2(i1,:) = Conjg(phaseM) *ZUR2(i1,:) 
 ZUL2(i1,:) = phaseM *ZUL2(i1,:) 
 End if 
  If (MFu2(i1).ne.MFu2(i1)) Then 
      Write(*,*) 'NaN appearing in '//NameOfUnit(Iname) 
      Call TerminateProgram 
    End If 
  If (Abs(MFu2(i1)).Le.MaxMassNumericalZero) MFu2(i1) = Abs(MFu2(i1))+1.E-10_dp 
  If (MFu2(i1).Le.0._dp) Then 
    If (ErrorLevel.Ge.0) Then 
      Write(10,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
      Write(10,*) 'a mass squarred is negative: ',i1,MFu2(i1) 
      Write(*,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
      Write(*,*) 'a mass squarred is negative: ',i1,MFu2(i1) 
      Call TerminateProgram 
    End If 
     Write(ErrCan,*) 'Warning from routine '//NameOfUnit(Iname) 
     Write(ErrCan,*) 'in the calculation of the masses' 
     Write(ErrCan,*) 'occurred a negative mass squared!' 
     Write(ErrCan,*) i1,MFu2(i1) 
     Write(*,*) 'Warning from routine '//NameOfUnit(Iname) 
     Write(*,*) 'in the calculation of the masses' 
     Write(*,*) 'occurred a negative mass squared!' 
     Write(*,*) i1,MFu2(i1) 
  MFu2(i1) = 1._dp 
   SignOfMassChanged = .True. 
! kont = -104 
 End if 
End Do 
 
If (ierr.Ne.0.) Then 
  Write(ErrCan,*) 'Warning from Subroutine CalculateMFu, ierr =',ierr 
  kont = ierr 
  Iname = Iname - 1 
  Return 
End If 


MFu = Sqrt( MFu2 ) 
ZUL = ZUL2 
ZUR = ZUR2 
Iname = Iname - 1 
 
End Subroutine CalculateMFu 

Subroutine CalculateMFe(Ye,v,ZEL,ZER,MFe,kont)

Real(dp),Intent(in) :: v

Complex(dp),Intent(in) :: Ye(3,3)

Integer, Intent(inout) :: kont 
Integer :: i1,i2,i3,i4, ierr 
Integer :: j1,j2,j3,j4 
Real(dp), Intent(out) :: MFe(3) 
 Complex(dp), Intent(out) :: ZEL(3,3), ZER(3,3) 
 
 Complex(dp) :: mat(3,3)=0._dp, mat2(3,3)=0._dp, phaseM 

Complex(dp) :: ZEL2(3,3), ZER2(3,3) 
 
 Real(dp) :: ZEL1(3,3), ZER1(3,3), test(2), MFe2(3) 
 
 Iname = Iname + 1 
NameOfUnit(Iname) = 'CalculateMFe'
 
MFe = 0._dp 
ZEL = 0._dp 
ZER = 0._dp 
mat(1,1) = 0._dp 
mat(1,1) = mat(1,1)+(v*Ye(1,1))/sqrt(2._dp)
mat(1,2) = 0._dp 
mat(1,2) = mat(1,2)+(v*Ye(2,1))/sqrt(2._dp)
mat(1,3) = 0._dp 
mat(1,3) = mat(1,3)+(v*Ye(3,1))/sqrt(2._dp)
mat(2,1) = 0._dp 
mat(2,1) = mat(2,1)+(v*Ye(1,2))/sqrt(2._dp)
mat(2,2) = 0._dp 
mat(2,2) = mat(2,2)+(v*Ye(2,2))/sqrt(2._dp)
mat(2,3) = 0._dp 
mat(2,3) = mat(2,3)+(v*Ye(3,2))/sqrt(2._dp)
mat(3,1) = 0._dp 
mat(3,1) = mat(3,1)+(v*Ye(1,3))/sqrt(2._dp)
mat(3,2) = 0._dp 
mat(3,2) = mat(3,2)+(v*Ye(2,3))/sqrt(2._dp)
mat(3,3) = 0._dp 
mat(3,3) = mat(3,3)+(v*Ye(3,3))/sqrt(2._dp)

 
mat2 = Matmul(Transpose(Conjg(mat)),mat) 
If (Maxval(Abs(Aimag(mat2))).Eq.0._dp) Then 
Call EigenSystem(Real(mat2,dp),MFe2,ZER1,ierr,test) 
ZER2 = ZER1 
Else 
Call EigenSystem(mat2,MFe2,ZER2,ierr,test) 
 End If 
 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
  End If 
  ierr = 0 
End If 
 
mat2 = Matmul(mat,Transpose(Conjg(mat))) 
If (Maxval(Abs(Aimag(mat2))).Eq.0._dp) Then 
Call EigenSystem (Real(mat2,dp),MFe2,ZEL1,ierr,test) 
                  
                  
ZEL2 = ZEL1 
Else 
Call EigenSystem(mat2,MFe2,ZEL2,ierr,test) 
 
 
End If 
ZEL2 = Conjg(ZEL2) 
 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
  End If 
  ierr = 0 
End If 
 
mat2 = Matmul(Matmul( Conjg(ZEL2),mat),Transpose( Conjg(ZER2))) 
Do i1=1,3
If (Abs(mat2(i1,i1)).gt.0._dp) Then 
phaseM = mat2(i1,i1) / Abs(mat2(i1,i1)) 
ZER2(i1,:) = phaseM *ZER2(i1,:) 
 End if 
End Do 
 
Do i1=1,3
If (Abs(ZER2(i1,i1)).gt.0._dp) Then 
phaseM = ZER2(i1,i1) / Abs(ZER2(i1,i1)) 
ZER2(i1,:) = Conjg(phaseM) *ZER2(i1,:) 
 ZEL2(i1,:) = phaseM *ZEL2(i1,:) 
 End if 
  If (MFe2(i1).ne.MFe2(i1)) Then 
      Write(*,*) 'NaN appearing in '//NameOfUnit(Iname) 
      Call TerminateProgram 
    End If 
  If (Abs(MFe2(i1)).Le.MaxMassNumericalZero) MFe2(i1) = Abs(MFe2(i1))+1.E-10_dp 
  If (MFe2(i1).Le.0._dp) Then 
    If (ErrorLevel.Ge.0) Then 
      Write(10,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
      Write(10,*) 'a mass squarred is negative: ',i1,MFe2(i1) 
      Write(*,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
      Write(*,*) 'a mass squarred is negative: ',i1,MFe2(i1) 
      Call TerminateProgram 
    End If 
     Write(ErrCan,*) 'Warning from routine '//NameOfUnit(Iname) 
     Write(ErrCan,*) 'in the calculation of the masses' 
     Write(ErrCan,*) 'occurred a negative mass squared!' 
     Write(ErrCan,*) i1,MFe2(i1) 
     Write(*,*) 'Warning from routine '//NameOfUnit(Iname) 
     Write(*,*) 'in the calculation of the masses' 
     Write(*,*) 'occurred a negative mass squared!' 
     Write(*,*) i1,MFe2(i1) 
  MFe2(i1) = 1._dp 
   SignOfMassChanged = .True. 
! kont = -104 
 End if 
End Do 
 
If (ierr.Ne.0.) Then 
  Write(ErrCan,*) 'Warning from Subroutine CalculateMFe, ierr =',ierr 
  kont = ierr 
  Iname = Iname - 1 
  Return 
End If 


MFe = Sqrt( MFe2 ) 
ZEL = ZEL2 
ZER = ZER2 
Iname = Iname - 1 
 
End Subroutine CalculateMFe 

Subroutine CalculateMFv(UV,MFv,kont)

Integer, Intent(inout) :: kont 
Integer :: i1,i2,i3,i4, ierr, pos 
Integer :: j1,j2,j3,j4 
Logical :: SecondDiagonalisationNeeded 
Real(dp), Intent(out) :: MFv(5) 
Complex(dp), Intent(out) ::  UV(5,5) 
                              
Complex(dp) :: mat(5,5), mat2(5,5), phaseM, E5(5) 

Real(dp) :: UVa(5,5), test(2), eig(5) 

Iname = Iname + 1 
NameOfUnit(Iname) = 'CalculateMFv'
 
mat(1,1) = 0._dp 
mat(1,2) = 0._dp 
mat(1,3) = 0._dp 
mat(1,4) = 0._dp 
mat(1,5) = 0._dp 
mat(2,2) = 0._dp 
mat(2,3) = 0._dp 
mat(2,4) = 0._dp 
mat(2,5) = 0._dp 
mat(3,3) = 0._dp 
mat(3,4) = 0._dp 
mat(3,5) = 0._dp 
mat(4,4) = 0._dp 
mat(4,5) = 0._dp 
mat(5,5) = 0._dp 

 
 Do i1=2,5
  Do i2 = 1, i1-1 
  mat(i1,i2) = mat(i2,i1) 
  End do 
End do 

 
If (Maxval(Abs(Aimag(mat))).Eq.0._dp) Then 
Call EigenSystem(Real(mat,dp),Eig,UVa,ierr,test) 
 
   Do i1=1,5
   If ((Eig(i1).Lt.0._dp).or.(Abs(eig(i1)).lt.1E-15)) Then 
    MFv(i1) = - Eig(i1) 
    UV(i1,:) = (0._dp,1._dp)*UVa(i1,:) 
   Else 
    MFv(i1) = Eig(i1) 
    UV(i1,:) = UVa(i1,:)
    End If 
   End Do 
 
Do i1=1,4
  Do i2=i1+1,5
    If (MFv(i1).Gt.MFv(i2)) Then 
      Eig(1) = MFv(i1) 
      MFv(i1) = MFv(i2) 
      MFv(i2) =  Eig(1) 
      E5 = UV(i1,:) 
      UV(i1,:) = UV(i2,:) 
      UV(i2,:) = E5
    End If 
   End Do 
End Do 
 
Else 
 
mat2 = Matmul( Transpose(Conjg( mat) ), mat ) 
Call Eigensystem(mat2, Eig, UV, ierr, test) 
mat2 = Matmul( Conjg(UV), Matmul( mat, Transpose( Conjg(UV)))) 
! Special efforts are needed for matrices like the Higgsinos one 
SecondDiagonalisationNeeded = .False. 
Do i1=1,5-1
If (MaxVal(Abs(mat2(i1,(i1+1):5))).gt.Abs(mat2(i1,i1))) SecondDiagonalisationNeeded = .True. 

  If (Eig(i1).ne.Eig(i1)) Then 
      Write(*,*) 'NaN appearing in '//NameOfUnit(Iname) 
      Call TerminateProgram 
    End If 
  If ((Abs(Eig(i1)).Le.MaxMassNumericalZero).and.(Eig(i1).lt.0._dp)) Eig(i1) = Abs(Eig(i1))+1.E-10_dp 
  If (Eig(i1).Le.0._dp) Then 
    If (ErrorLevel.Ge.0) Then 
      Write(10,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
      Write(10,*) 'a mass squarred is negative: ',i1,Eig(i1) 
      Write(*,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
      Write(*,*) 'a mass squarred is negative: ',i1,Eig(i1) 
      Call TerminateProgram 
    End If 
     Write(ErrCan,*) 'Warning from routine '//NameOfUnit(Iname) 
     Write(ErrCan,*) 'in the calculation of the masses' 
     Write(ErrCan,*) 'occurred a negative mass squared!' 
     Write(ErrCan,*) i1,Eig(i1) 
     Write(*,*) 'Warning from routine '//NameOfUnit(Iname) 
     Write(*,*) 'in the calculation of the masses' 
     Write(*,*) 'occurred a negative mass squared!' 
     Write(*,*) i1,Eig(i1) 
  Eig(i1) = 1._dp 
   SignOfMassChanged = .True. 
! kont = -104 
 End if 
End do 
If (SecondDiagonalisationNeeded) Then 
Call EigenSystem(Real(mat2,dp),Eig,UVa,ierr,test) 
 
     UV = MatMul(UV,UVa)
  Do i1=1,5
   If ((Eig(i1).Lt.0._dp).or.(Abs(eig(i1)).lt.1E-15)) Then 
    MFv(i1) = - Eig(i1) 
    UV(i1,:) = (0._dp,1._dp)*UVa(i1,:) 
   Else 
    MFv(i1) = Eig(i1) 
    UV(i1,:) = UVa(i1,:)
    End If 
   End Do 
 
Else 
Do i1=1,5
  If (Eig(i1).ne.Eig(i1)) Then 
      Write(*,*) 'NaN appearing in '//NameOfUnit(Iname) 
      Call TerminateProgram 
    End If 
If (Abs(mat2(i1,i1)).gt.0._dp) Then 
  phaseM = Sqrt( mat2(i1,i1) / Abs(mat2(i1,i1))) 
  UV(i1,:)= phaseM * UV(i1,:) 
End if 
  If ((Abs(Eig(i1)).Le.MaxMassNumericalZero).and.(Eig(i1).lt.0._dp)) Eig(i1) = Abs(Eig(i1))+1.E-10_dp 
  If (Eig(i1).Le.0._dp) Then 
    If (ErrorLevel.Ge.0) Then 
      Write(10,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
      Write(10,*) 'a mass squarred is negative: ',i1,Eig(i1) 
      Write(*,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
      Write(*,*) 'a mass squarred is negative: ',i1,Eig(i1) 
      Call TerminateProgram 
    End If 
     Write(ErrCan,*) 'Warning from routine '//NameOfUnit(Iname) 
     Write(ErrCan,*) 'in the calculation of the masses' 
     Write(ErrCan,*) 'occurred a negative mass squared!' 
     Write(ErrCan,*) i1,Eig(i1) 
     Write(*,*) 'Warning from routine '//NameOfUnit(Iname) 
     Write(*,*) 'in the calculation of the masses' 
     Write(*,*) 'occurred a negative mass squared!' 
     Write(*,*) i1,Eig(i1) 
  Eig(i1) = 1._dp 
   SignOfMassChanged = .True. 
! kont = -104 
 End if 
End Do 
MFv = Sqrt( Eig ) 
 
End if ! Second diagonalisation 
End If 
 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
  End If 
  ierr = 0 
End If 
 
If (ierr.Ne.0.) Then 
  Write(ErrCan,*) 'Warning from Subroutine CalculateMFv, ierr =',ierr 
  kont = ierr 
  Iname = Iname - 1 
  Return 
End If 


Iname = Iname - 1 
 
End Subroutine CalculateMFv 

Subroutine CalculateMFeD(YB1,YC1,YRH,YRG,v,vX,ZELD,ZERD,MFeD,kont)

Real(dp),Intent(in) :: YB1,YC1,YRH,v,vX

Complex(dp),Intent(in) :: YRG

Integer, Intent(inout) :: kont 
Integer :: i1,i2,i3,i4, ierr 
Integer :: j1,j2,j3,j4 
Real(dp), Intent(out) :: MFeD(2) 
 Complex(dp), Intent(out) :: ZELD(2,2), ZERD(2,2) 
 
 Complex(dp) :: mat(2,2)=0._dp, mat2(2,2)=0._dp, phaseM 

Complex(dp) :: ZELD2(2,2), ZERD2(2,2) 
 
 Real(dp) :: ZELD1(2,2), ZERD1(2,2), test(2), MFeD2(2) 
 
 Iname = Iname + 1 
NameOfUnit(Iname) = 'CalculateMFeD'
 
MFeD = 0._dp 
ZELD = 0._dp 
ZERD = 0._dp 
mat(1,1) = 0._dp 
mat(1,1) = mat(1,1)+(v*YRG)/sqrt(2._dp)
mat(1,2) = 0._dp 
mat(1,2) = mat(1,2)+(vX*YC1)/sqrt(2._dp)
mat(2,1) = 0._dp 
mat(2,1) = mat(2,1)+(vX*YB1)/sqrt(2._dp)
mat(2,2) = 0._dp 
mat(2,2) = mat(2,2)+(v*YRH)/sqrt(2._dp)

 
mat2 = Matmul(Transpose(Conjg(mat)),mat) 
If (Maxval(Abs(Aimag(mat2))).Eq.0._dp) Then 
Call EigenSystem(Real(mat2,dp),MFeD2,ZERD1,ierr,test) 
ZERD2 = ZERD1 
Else 
Call EigenSystem(mat2,MFeD2,ZERD2,ierr,test) 
 End If 
 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
  End If 
  ierr = 0 
End If 
 
mat2 = Matmul(mat,Transpose(Conjg(mat))) 
If (Maxval(Abs(Aimag(mat2))).Eq.0._dp) Then 
Call EigenSystem (Real(mat2,dp),MFeD2,ZELD1,ierr,test) 
                  
                  
ZELD2 = ZELD1 
Else 
Call EigenSystem(mat2,MFeD2,ZELD2,ierr,test) 
 
 
End If 
ZELD2 = Conjg(ZELD2) 
 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
  End If 
  ierr = 0 
End If 
 
mat2 = Matmul(Matmul( Conjg(ZELD2),mat),Transpose( Conjg(ZERD2))) 
Do i1=1,2
If (Abs(mat2(i1,i1)).gt.0._dp) Then 
phaseM = mat2(i1,i1) / Abs(mat2(i1,i1)) 
ZERD2(i1,:) = phaseM *ZERD2(i1,:) 
 End if 
End Do 
 
Do i1=1,2
If (Abs(ZERD2(i1,i1)).gt.0._dp) Then 
phaseM = ZERD2(i1,i1) / Abs(ZERD2(i1,i1)) 
ZERD2(i1,:) = Conjg(phaseM) *ZERD2(i1,:) 
 ZELD2(i1,:) = phaseM *ZELD2(i1,:) 
 End if 
  If (MFeD2(i1).ne.MFeD2(i1)) Then 
      Write(*,*) 'NaN appearing in '//NameOfUnit(Iname) 
      Call TerminateProgram 
    End If 
  If (Abs(MFeD2(i1)).Le.MaxMassNumericalZero) MFeD2(i1) = Abs(MFeD2(i1))+1.E-10_dp 
  If (MFeD2(i1).Le.0._dp) Then 
    If (ErrorLevel.Ge.0) Then 
      Write(10,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
      Write(10,*) 'a mass squarred is negative: ',i1,MFeD2(i1) 
      Write(*,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
      Write(*,*) 'a mass squarred is negative: ',i1,MFeD2(i1) 
      Call TerminateProgram 
    End If 
     Write(ErrCan,*) 'Warning from routine '//NameOfUnit(Iname) 
     Write(ErrCan,*) 'in the calculation of the masses' 
     Write(ErrCan,*) 'occurred a negative mass squared!' 
     Write(ErrCan,*) i1,MFeD2(i1) 
     Write(*,*) 'Warning from routine '//NameOfUnit(Iname) 
     Write(*,*) 'in the calculation of the masses' 
     Write(*,*) 'occurred a negative mass squared!' 
     Write(*,*) i1,MFeD2(i1) 
  MFeD2(i1) = 1._dp 
   SignOfMassChanged = .True. 
! kont = -104 
 End if 
End Do 
 
If (ierr.Ne.0.) Then 
  Write(ErrCan,*) 'Warning from Subroutine CalculateMFeD, ierr =',ierr 
  kont = ierr 
  Iname = Iname - 1 
  Return 
End If 


MFeD = Sqrt( MFeD2 ) 
ZELD = ZELD2 
ZERD = ZERD2 
Iname = Iname - 1 
 
End Subroutine CalculateMFeD 

Subroutine CalculateVPVZVZp(g1,g2,g1p,gBY,gYB,v,vX,ZZ,MVZ,MVZp,MVZ2,MVZp2,            & 
& TW,TWp,kont)

Real(dp), Intent(in) :: g1,g2,g1p,gBY,gYB,v,vX

Real(dp), Intent(out) :: TW,TWp

Integer, Intent(inout) :: kont 
Integer :: i1,i2,i3,i4, ierr 
Integer :: j1,j2,j3,j4 
Real(dp), Intent(out) :: MVZ, MVZ2
Real(dp), Intent(out) :: MVZp, MVZp2
Real(dp) :: VPVZVZp2(3),VPVZVZp(3)  

Real(dp), Intent(out) :: ZZ(3,3) 
 
Real(dp) :: mat(3,3)  

Real(dp) ::  test(2) 

Iname = Iname + 1 
NameOfUnit(Iname) = 'CalculateVPVZVZp'
 
mat(1,1) = 0._dp 
mat(1,1) = mat(1,1)+(g1**2*v**2)/4._dp
mat(1,1) = mat(1,1)+25*gYB**2*vX**2
mat(1,2) = 0._dp 
mat(1,2) = mat(1,2)-(g1*g2*v**2)/4._dp
mat(1,3) = 0._dp 
mat(1,3) = mat(1,3)+(g1*gBY*v**2)/4._dp
mat(1,3) = mat(1,3)+25*g1p*gYB*vX**2
mat(2,2) = 0._dp 
mat(2,2) = mat(2,2)+(g2**2*v**2)/4._dp
mat(2,3) = 0._dp 
mat(2,3) = mat(2,3)-(gBY*g2*v**2)/4._dp
mat(3,3) = 0._dp 
mat(3,3) = mat(3,3)+(gBY**2*v**2)/4._dp
mat(3,3) = mat(3,3)+25*g1p**2*vX**2

 
 Do i1=2,3
  Do i2 = 1, i1-1 
  mat(i1,i2) = mat(i2,i1) 
  End do 
End do 

 
Call EigenSystem(mat,VPVZVZp2,ZZ,ierr,test) 
 
 
ZZ = Transpose(ZZ) 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
  End If 
  ierr = 0 
End If 
 
If ((ierr.Ne.0.).And.(ErrorLevel.Ge.-1)) Then 
  Write(10,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
  Write(10,*) 'Diagonalization failed, ierr : ',ierr 
  kont = ierr 
  Iname = Iname - 1 
  Return 
End If 


Do i1=1,3
  If (Abs(VPVZVZp2(i1)).Le.1.E-10_dp*(Maxval(VPVZVZp2))) VPVZVZp2(i1) = 1.E-10_dp 
  If (VPVZVZp2(i1).ne.VPVZVZp2(i1)) Then 
      Write(*,*) 'NaN appearing in '//NameOfUnit(Iname) 
      Call TerminateProgram 
    End If 
  If (VPVZVZp2(i1).Ge.0._dp) Then 
  VPVZVZp(i1) =Sqrt(VPVZVZp2(i1) ) 
  Else 
    If (ErrorLevel.Ge.0) Then 
      Write(10,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
      Write(10,*) 'a mass squarred is negative: ',i1,VPVZVZp2(i1) 
    End If 
  VPVZVZp(i1)= 1._dp 
  VPVZVZp2(i1)= 1._dp  
     Write(ErrCan,*) 'Warning from routine '//NameOfUnit(Iname) 
     Write(ErrCan,*) 'in the calculation of the masses' 
     Write(ErrCan,*) 'occurred a negative mass squared!' 
     Write(ErrCan,*) i1,VPVZVZp2(i1) 
     Write(*,*) 'Warning from routine '//NameOfUnit(Iname) 
     Write(*,*) 'in the calculation of the masses' 
     Write(*,*) 'occurred a negative mass squared!' 
     Write(*,*) i1,VPVZVZp2(i1) 
  VPVZVZp(i1)= 1._dp 
  VPVZVZp2(i1) = 1._dp  
   SignOfMassChanged = .True. 
! kont = -104 
 End if 
End Do 
 
MVZ = VPVZVZp(2) 
MVZ2 = VPVZVZp2(2) 
MVZp = VPVZVZp(3) 
MVZp2 = VPVZVZp2(3) 
TW = ACos(Abs(ZZ(1,1)))
TWp = ACos(Abs(ZZ(3,3)))

 Iname = Iname - 1 
 
End Subroutine CalculateVPVZVZp 

Subroutine CalculateVWm(g2,v,ZW,MVWm,MVWm2,kont)

Real(dp), Intent(in) :: g2,v

Integer, Intent(inout) :: kont 
Integer :: i1,i2,i3,i4, ierr 
Integer :: j1,j2,j3,j4 
Real(dp), Intent(out) :: MVWm, MVWm2
Real(dp) :: VWm2(2),VWm(2)  

Complex(dp), Intent(out) :: ZW(2,2) 
 
Complex(dp) :: mat(2,2)  

Real(dp) ::  test(2) 

Iname = Iname + 1 
NameOfUnit(Iname) = 'CalculateVWm'
 
mat(1,1) = 0._dp 
mat(1,1) = mat(1,1)+(g2**2*v**2)/4._dp
mat(1,2) = 0._dp 
mat(2,2) = 0._dp 
mat(2,2) = mat(2,2)+(g2**2*v**2)/4._dp

 
 Do i1=2,2
  Do i2 = 1, i1-1 
  mat(i1,i2) = Conjg(mat(i2,i1)) 
  End do 
End do 

 
Call EigenSystem(mat,VWm2,ZW,ierr,test) 
 
 
ZW = Transpose(ZW) 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
  End If 
  ierr = 0 
End If 
 
If ((ierr.Ne.0.).And.(ErrorLevel.Ge.-1)) Then 
  Write(10,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
  Write(10,*) 'Diagonalization failed, ierr : ',ierr 
  kont = ierr 
  Iname = Iname - 1 
  Return 
End If 


Do i1=1,2
  If (Abs(VWm2(i1)).Le.1.E-10_dp*(Maxval(VWm2))) VWm2(i1) = 1.E-10_dp 
  If (VWm2(i1).ne.VWm2(i1)) Then 
      Write(*,*) 'NaN appearing in '//NameOfUnit(Iname) 
      Call TerminateProgram 
    End If 
  If (VWm2(i1).Ge.0._dp) Then 
  VWm(i1) =Sqrt(VWm2(i1) ) 
  Else 
    If (ErrorLevel.Ge.0) Then 
      Write(10,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
      Write(10,*) 'a mass squarred is negative: ',i1,VWm2(i1) 
    End If 
  VWm(i1)= 1._dp 
  VWm2(i1)= 1._dp  
     Write(ErrCan,*) 'Warning from routine '//NameOfUnit(Iname) 
     Write(ErrCan,*) 'in the calculation of the masses' 
     Write(ErrCan,*) 'occurred a negative mass squared!' 
     Write(ErrCan,*) i1,VWm2(i1) 
     Write(*,*) 'Warning from routine '//NameOfUnit(Iname) 
     Write(*,*) 'in the calculation of the masses' 
     Write(*,*) 'occurred a negative mass squared!' 
     Write(*,*) i1,VWm2(i1) 
  VWm(i1)= 1._dp 
  VWm2(i1) = 1._dp  
   SignOfMassChanged = .True. 
! kont = -104 
 End if 
End Do 
 
MVWm = VWm(1) 
MVWm2 = VWm2(1) 

 Iname = Iname - 1 
 
End Subroutine CalculateVWm 

Subroutine CalculateMhhEffPot(muh,lamh,mui,lam2,lam3,v,vX,ZH,Mhh,Mhh2,kont)

Real(dp), Intent(in) :: v,vX

Complex(dp), Intent(in) :: muh,lamh,mui,lam2,lam3

Integer, Intent(inout) :: kont 
Integer :: i1,i2,i3,i4, ierr 
Integer :: j1,j2,j3,j4, pos 
Real(dp), Intent(out) :: Mhh(2), Mhh2(2) 
Real(dp), Intent(out) :: ZH(2,2) 
 
Real(dp) :: mat(2,2)  

Real(dp) :: Mhh2temp(2), Q2 
Real(dp) :: ZHtemp(2,2),ZHtemp2(2,2) 
 
Real(dp) ::  test(2) 

Iname = Iname + 1 
NameOfUnit(Iname) = 'CalculateMhh'
 
mat(1,1) = 0._dp 
mat(1,1) = mat(1,1)+muh
mat(1,1) = mat(1,1)-3*lamh*v**2
mat(1,1) = mat(1,1)-(lam3*vX**2)/2._dp
mat(1,2) = 0._dp 
mat(1,2) = mat(1,2)-(lam3*v*vX)
mat(2,2) = 0._dp 
mat(2,2) = mat(2,2)+mui
mat(2,2) = mat(2,2)-(lam3*v**2)/2._dp
mat(2,2) = mat(2,2)-3*lam2*vX**2

 
 Do i1=2,2
  Do i2 = 1, i1-1 
  mat(i1,i2) = mat(i2,i1) 
  End do 
End do 

 
Call EigenSystem(mat,Mhh2,ZH,ierr,test) 
 
 
! Fix phases
Do i1=1,2
  pos=Maxloc(Abs(ZH(i1,:)),1)
  If (Real(ZH(i1,pos),dp).lt.0._dp) Then
    ZH(i1,:)=-ZH(i1,:)
  End if
End do
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
  End If 
  ierr = 0 
End If 
 
If ((ierr.Ne.0.).And.(ErrorLevel.Ge.-1)) Then 
  Write(10,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
  Write(10,*) 'Diagonalization failed, ierr : ',ierr 
  kont = ierr 
  Iname = Iname - 1 
  Return 
End If 


Do i1=1,2
  If (Mhh2(i1).ne.Mhh2(i1)) Then 
      Write(*,*) 'NaN appearing in '//NameOfUnit(Iname) 
      Call TerminateProgram 
    End If 
  If (Mhh2(i1).Ge.0._dp) Then 
  Mhh(i1)=Sqrt(Mhh2(i1) ) 
  Else 
  Mhh(i1) = 1._dp 
! kont = -104 
 End if 
End Do 
Iname = Iname - 1 
 
End Subroutine CalculateMhhEffPot 

Subroutine CalculateMAhEffPot(g1,g2,g1p,gBY,gYB,muh,lamh,mui,lam2,lam3,               & 
& v,vX,TW,TWp,ZA,MAh,MAh2,kont)

Real(dp), Intent(in) :: g1,g2,g1p,gBY,gYB,v,vX,TW,TWp

Complex(dp), Intent(in) :: muh,lamh,mui,lam2,lam3

Integer, Intent(inout) :: kont 
Integer :: i1,i2,i3,i4, ierr, pos 
Integer :: j1,j2,j3,j4 
Real(dp), Intent(out) :: MAh(2), MAh2(2) 
Real(dp), Intent(out) :: ZA(2,2) 
 
Real(dp) :: ZAFIX(2,2) 
 
Real(dp) :: mat(2,2)  

Real(dp) ::  test(2), Q2 

Iname = Iname + 1 
NameOfUnit(Iname) = 'CalculateMAh'
 
mat(1,1) = 0._dp 
mat(1,1) = mat(1,1)+muh
mat(1,1) = mat(1,1)-(lamh*v**2)
mat(1,1) = mat(1,1)-(lam3*vX**2)/2._dp
mat(1,2) = 0._dp 
mat(2,2) = 0._dp 
mat(2,2) = mat(2,2)+mui
mat(2,2) = mat(2,2)-(lam3*v**2)/2._dp
mat(2,2) = mat(2,2)-(lam2*vX**2)

 
 Do i1=2,2
  Do i2 = 1, i1-1 
  mat(i1,i2) = mat(i2,i1) 
  End do 
End do 

 
Call EigenSystem(mat,MAh2,ZA,ierr,test) 
 
 
! Fix phases
Do i1=1,2
  pos=Maxloc(Abs(ZA(i1,:)),1)
  If (Real(ZA(i1,pos),dp).lt.0._dp) Then
    ZA(i1,:)=-ZA(i1,:)
  End if
End do
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
  End If 
  ierr = 0 
End If 
 
If ((ierr.Ne.0.).And.(ErrorLevel.Ge.-1)) Then 
  Write(10,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
  Write(10,*) 'Diagonalization failed, ierr : ',ierr 
  kont = ierr 
  Iname = Iname - 1 
  Return 
End If 


Do i1=1,2
  If (MAh2(i1).ne.MAh2(i1)) Then 
      Write(*,*) 'NaN appearing in '//NameOfUnit(Iname) 
      Call TerminateProgram 
    End If 
  If (MAh2(i1).Ge.0._dp) Then 
  MAh(i1)=Sqrt(MAh2(i1) ) 
  Else 
  MAh(i1) = 1._dp 
! kont = -104 
 End if 
End Do 
Iname = Iname - 1 
 
End Subroutine CalculateMAhEffPot 

Subroutine CalculateMHmEffPot(g2,muh,lamh,lam3,m1,lam6,m2,lam7,Fi,v,vX,               & 
& ZP,MHm,MHm2,kont)

Real(dp), Intent(in) :: g2,Fi(2,2),v,vX

Complex(dp), Intent(in) :: muh,lamh,lam3,m1(2,2),lam6(2,2),m2(2,2),lam7(2,2)

Integer, Intent(inout) :: kont 
Integer :: i1,i2,i3,i4, ierr, pos 
Integer :: j1,j2,j3,j4 
Real(dp), Intent(out) :: MHm(5), MHm2(5) 
Real(dp), Intent(out) :: ZP(5,5) 
 
Real(dp) :: ZPFIX(5,5) 
 
Real(dp) :: mat(5,5)  

Real(dp) ::  test(2), Q2 

Iname = Iname + 1 
NameOfUnit(Iname) = 'CalculateMHm'
 
mat(1,1) = 0._dp 
mat(1,1) = mat(1,1)+muh
mat(1,1) = mat(1,1)-(lamh*v**2)
mat(1,1) = mat(1,1)-(lam3*vX**2)/2._dp
mat(1,2) = 0._dp 
mat(1,3) = 0._dp 
mat(1,4) = 0._dp 
mat(1,5) = 0._dp 
mat(2,2) = 0._dp 
mat(2,2) = mat(2,2)+(v**2*lam6(1,1))/2._dp
mat(2,2) = mat(2,2)+m1(1,1)
mat(2,3) = 0._dp 
mat(2,3) = mat(2,3)+(v**2*lam6(1,2))/2._dp
mat(2,3) = mat(2,3)+m1(1,2)
mat(2,4) = 0._dp 
mat(2,4) = mat(2,4)+(vX*Fi(1,1))/sqrt(2._dp)
mat(2,5) = 0._dp 
mat(2,5) = mat(2,5)+(vX*Fi(1,2))/sqrt(2._dp)
mat(3,3) = 0._dp 
mat(3,3) = mat(3,3)+(v**2*lam6(2,2))/2._dp
mat(3,3) = mat(3,3)+m1(2,2)
mat(3,4) = 0._dp 
mat(3,4) = mat(3,4)+(vX*Fi(2,1))/sqrt(2._dp)
mat(3,5) = 0._dp 
mat(3,5) = mat(3,5)+(vX*Fi(2,2))/sqrt(2._dp)
mat(4,4) = 0._dp 
mat(4,4) = mat(4,4)+(v**2*lam7(1,1))/2._dp
mat(4,4) = mat(4,4)+m2(1,1)
mat(4,5) = 0._dp 
mat(4,5) = mat(4,5)+(v**2*lam7(1,2))/2._dp
mat(4,5) = mat(4,5)+m2(1,2)
mat(5,5) = 0._dp 
mat(5,5) = mat(5,5)+(v**2*lam7(2,2))/2._dp
mat(5,5) = mat(5,5)+m2(2,2)

 
 Do i1=2,5
  Do i2 = 1, i1-1 
  mat(i1,i2) = mat(i2,i1) 
  End do 
End do 

 
Call EigenSystem(mat,MHm2,ZP,ierr,test) 
 
 
! Fix phases
Do i1=1,5
  pos=Maxloc(Abs(ZP(i1,:)),1)
  If (Real(ZP(i1,pos),dp).lt.0._dp) Then
    ZP(i1,:)=-ZP(i1,:)
  End if
End do
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
  End If 
  ierr = 0 
End If 
 
If ((ierr.Ne.0.).And.(ErrorLevel.Ge.-1)) Then 
  Write(10,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
  Write(10,*) 'Diagonalization failed, ierr : ',ierr 
  kont = ierr 
  Iname = Iname - 1 
  Return 
End If 


Do i1=1,5
  If (MHm2(i1).ne.MHm2(i1)) Then 
      Write(*,*) 'NaN appearing in '//NameOfUnit(Iname) 
      Call TerminateProgram 
    End If 
  If (MHm2(i1).Ge.0._dp) Then 
  MHm(i1)=Sqrt(MHm2(i1) ) 
  Else 
  MHm(i1) = 1._dp 
! kont = -104 
 End if 
End Do 
Iname = Iname - 1 
 
End Subroutine CalculateMHmEffPot 

Subroutine CalculateMFdEffPot(Yd,v,ZDL,ZDR,MFd,kont)

Real(dp),Intent(in) :: v

Complex(dp),Intent(in) :: Yd(3,3)

Integer, Intent(inout) :: kont 
Integer :: i1,i2,i3,i4, ierr 
Integer :: j1,j2,j3,j4 
Real(dp), Intent(out) :: MFd(3) 
 Complex(dp), Intent(out) :: ZDL(3,3), ZDR(3,3) 
 
 Complex(dp) :: mat(3,3)=0._dp, mat2(3,3)=0._dp, phaseM 

Complex(dp) :: ZDL2(3,3), ZDR2(3,3) 
 
 Real(dp) :: ZDL1(3,3), ZDR1(3,3), test(2), MFd2(3) 
 
 Iname = Iname + 1 
NameOfUnit(Iname) = 'CalculateMFd'
 
MFd = 0._dp 
ZDL = 0._dp 
ZDR = 0._dp 
mat(1,1) = 0._dp 
mat(1,1) = mat(1,1)+(v*Yd(1,1))/sqrt(2._dp)
mat(1,2) = 0._dp 
mat(1,2) = mat(1,2)+(v*Yd(2,1))/sqrt(2._dp)
mat(1,3) = 0._dp 
mat(1,3) = mat(1,3)+(v*Yd(3,1))/sqrt(2._dp)
mat(2,1) = 0._dp 
mat(2,1) = mat(2,1)+(v*Yd(1,2))/sqrt(2._dp)
mat(2,2) = 0._dp 
mat(2,2) = mat(2,2)+(v*Yd(2,2))/sqrt(2._dp)
mat(2,3) = 0._dp 
mat(2,3) = mat(2,3)+(v*Yd(3,2))/sqrt(2._dp)
mat(3,1) = 0._dp 
mat(3,1) = mat(3,1)+(v*Yd(1,3))/sqrt(2._dp)
mat(3,2) = 0._dp 
mat(3,2) = mat(3,2)+(v*Yd(2,3))/sqrt(2._dp)
mat(3,3) = 0._dp 
mat(3,3) = mat(3,3)+(v*Yd(3,3))/sqrt(2._dp)

 
mat2 = Matmul(Transpose(Conjg(mat)),mat) 
If (Maxval(Abs(Aimag(mat2))).Eq.0._dp) Then 
Call EigenSystem(Real(mat2,dp),MFd2,ZDR1,ierr,test) 
ZDR2 = ZDR1 
Else 
Call EigenSystem(mat2,MFd2,ZDR2,ierr,test) 
 End If 
 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
  End If 
  ierr = 0 
End If 
 
mat2 = Matmul(mat,Transpose(Conjg(mat))) 
If (Maxval(Abs(Aimag(mat2))).Eq.0._dp) Then 
Call EigenSystem (Real(mat2,dp),MFd2,ZDL1,ierr,test) 
                  
                  
ZDL2 = ZDL1 
Else 
Call EigenSystem(mat2,MFd2,ZDL2,ierr,test) 
 
 
End If 
ZDL2 = Conjg(ZDL2) 
 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
  End If 
  ierr = 0 
End If 
 
mat2 = Matmul(Matmul( Conjg(ZDL2),mat),Transpose( Conjg(ZDR2))) 
Do i1=1,3
If (Abs(mat2(i1,i1)).gt.0._dp) Then 
phaseM = mat2(i1,i1) / Abs(mat2(i1,i1)) 
ZDR2(i1,:) = phaseM *ZDR2(i1,:) 
 End if 
End Do 
 
Do i1=1,3
If (Abs(ZDR2(i1,i1)).gt.0._dp) Then 
phaseM = ZDR2(i1,i1) / Abs(ZDR2(i1,i1)) 
ZDR2(i1,:) = Conjg(phaseM) *ZDR2(i1,:) 
 ZDL2(i1,:) = phaseM *ZDL2(i1,:) 
 End if 
  If (MFd2(i1).ne.MFd2(i1)) Then 
      Write(*,*) 'NaN appearing in '//NameOfUnit(Iname) 
      Call TerminateProgram 
    End If 
  If (Abs(MFd2(i1)).Le.MaxMassNumericalZero) MFd2(i1) = Abs(MFd2(i1))+1.E-10_dp 
  If (MFd2(i1).Le.0._dp) Then 
! kont = -104 
 End if 
End Do 
 
If (ierr.Ne.0.) Then 
  Write(ErrCan,*) 'Warning from Subroutine CalculateMFd, ierr =',ierr 
  kont = ierr 
  Iname = Iname - 1 
  Return 
End If 


MFd = Sqrt( MFd2 ) 
ZDL = ZDL2 
ZDR = ZDR2 
Iname = Iname - 1 
 
End Subroutine CalculateMFdEffPot 

Subroutine CalculateMFuEffPot(Yu,v,ZUL,ZUR,MFu,kont)

Real(dp),Intent(in) :: v

Complex(dp),Intent(in) :: Yu(3,3)

Integer, Intent(inout) :: kont 
Integer :: i1,i2,i3,i4, ierr 
Integer :: j1,j2,j3,j4 
Real(dp), Intent(out) :: MFu(3) 
 Complex(dp), Intent(out) :: ZUL(3,3), ZUR(3,3) 
 
 Complex(dp) :: mat(3,3)=0._dp, mat2(3,3)=0._dp, phaseM 

Complex(dp) :: ZUL2(3,3), ZUR2(3,3) 
 
 Real(dp) :: ZUL1(3,3), ZUR1(3,3), test(2), MFu2(3) 
 
 Iname = Iname + 1 
NameOfUnit(Iname) = 'CalculateMFu'
 
MFu = 0._dp 
ZUL = 0._dp 
ZUR = 0._dp 
mat(1,1) = 0._dp 
mat(1,1) = mat(1,1)-((v*Yu(1,1))/sqrt(2._dp))
mat(1,2) = 0._dp 
mat(1,2) = mat(1,2)-((v*Yu(2,1))/sqrt(2._dp))
mat(1,3) = 0._dp 
mat(1,3) = mat(1,3)-((v*Yu(3,1))/sqrt(2._dp))
mat(2,1) = 0._dp 
mat(2,1) = mat(2,1)-((v*Yu(1,2))/sqrt(2._dp))
mat(2,2) = 0._dp 
mat(2,2) = mat(2,2)-((v*Yu(2,2))/sqrt(2._dp))
mat(2,3) = 0._dp 
mat(2,3) = mat(2,3)-((v*Yu(3,2))/sqrt(2._dp))
mat(3,1) = 0._dp 
mat(3,1) = mat(3,1)-((v*Yu(1,3))/sqrt(2._dp))
mat(3,2) = 0._dp 
mat(3,2) = mat(3,2)-((v*Yu(2,3))/sqrt(2._dp))
mat(3,3) = 0._dp 
mat(3,3) = mat(3,3)-((v*Yu(3,3))/sqrt(2._dp))

 
mat2 = Matmul(Transpose(Conjg(mat)),mat) 
If (Maxval(Abs(Aimag(mat2))).Eq.0._dp) Then 
Call EigenSystem(Real(mat2,dp),MFu2,ZUR1,ierr,test) 
ZUR2 = ZUR1 
Else 
Call EigenSystem(mat2,MFu2,ZUR2,ierr,test) 
 End If 
 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
  End If 
  ierr = 0 
End If 
 
mat2 = Matmul(mat,Transpose(Conjg(mat))) 
If (Maxval(Abs(Aimag(mat2))).Eq.0._dp) Then 
Call EigenSystem (Real(mat2,dp),MFu2,ZUL1,ierr,test) 
                  
                  
ZUL2 = ZUL1 
Else 
Call EigenSystem(mat2,MFu2,ZUL2,ierr,test) 
 
 
End If 
ZUL2 = Conjg(ZUL2) 
 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
  End If 
  ierr = 0 
End If 
 
mat2 = Matmul(Matmul( Conjg(ZUL2),mat),Transpose( Conjg(ZUR2))) 
Do i1=1,3
If (Abs(mat2(i1,i1)).gt.0._dp) Then 
phaseM = mat2(i1,i1) / Abs(mat2(i1,i1)) 
ZUR2(i1,:) = phaseM *ZUR2(i1,:) 
 End if 
End Do 
 
Do i1=1,3
If (Abs(ZUR2(i1,i1)).gt.0._dp) Then 
phaseM = ZUR2(i1,i1) / Abs(ZUR2(i1,i1)) 
ZUR2(i1,:) = Conjg(phaseM) *ZUR2(i1,:) 
 ZUL2(i1,:) = phaseM *ZUL2(i1,:) 
 End if 
  If (MFu2(i1).ne.MFu2(i1)) Then 
      Write(*,*) 'NaN appearing in '//NameOfUnit(Iname) 
      Call TerminateProgram 
    End If 
  If (Abs(MFu2(i1)).Le.MaxMassNumericalZero) MFu2(i1) = Abs(MFu2(i1))+1.E-10_dp 
  If (MFu2(i1).Le.0._dp) Then 
! kont = -104 
 End if 
End Do 
 
If (ierr.Ne.0.) Then 
  Write(ErrCan,*) 'Warning from Subroutine CalculateMFu, ierr =',ierr 
  kont = ierr 
  Iname = Iname - 1 
  Return 
End If 


MFu = Sqrt( MFu2 ) 
ZUL = ZUL2 
ZUR = ZUR2 
Iname = Iname - 1 
 
End Subroutine CalculateMFuEffPot 

Subroutine CalculateMFeEffPot(Ye,v,ZEL,ZER,MFe,kont)

Real(dp),Intent(in) :: v

Complex(dp),Intent(in) :: Ye(3,3)

Integer, Intent(inout) :: kont 
Integer :: i1,i2,i3,i4, ierr 
Integer :: j1,j2,j3,j4 
Real(dp), Intent(out) :: MFe(3) 
 Complex(dp), Intent(out) :: ZEL(3,3), ZER(3,3) 
 
 Complex(dp) :: mat(3,3)=0._dp, mat2(3,3)=0._dp, phaseM 

Complex(dp) :: ZEL2(3,3), ZER2(3,3) 
 
 Real(dp) :: ZEL1(3,3), ZER1(3,3), test(2), MFe2(3) 
 
 Iname = Iname + 1 
NameOfUnit(Iname) = 'CalculateMFe'
 
MFe = 0._dp 
ZEL = 0._dp 
ZER = 0._dp 
mat(1,1) = 0._dp 
mat(1,1) = mat(1,1)+(v*Ye(1,1))/sqrt(2._dp)
mat(1,2) = 0._dp 
mat(1,2) = mat(1,2)+(v*Ye(2,1))/sqrt(2._dp)
mat(1,3) = 0._dp 
mat(1,3) = mat(1,3)+(v*Ye(3,1))/sqrt(2._dp)
mat(2,1) = 0._dp 
mat(2,1) = mat(2,1)+(v*Ye(1,2))/sqrt(2._dp)
mat(2,2) = 0._dp 
mat(2,2) = mat(2,2)+(v*Ye(2,2))/sqrt(2._dp)
mat(2,3) = 0._dp 
mat(2,3) = mat(2,3)+(v*Ye(3,2))/sqrt(2._dp)
mat(3,1) = 0._dp 
mat(3,1) = mat(3,1)+(v*Ye(1,3))/sqrt(2._dp)
mat(3,2) = 0._dp 
mat(3,2) = mat(3,2)+(v*Ye(2,3))/sqrt(2._dp)
mat(3,3) = 0._dp 
mat(3,3) = mat(3,3)+(v*Ye(3,3))/sqrt(2._dp)

 
mat2 = Matmul(Transpose(Conjg(mat)),mat) 
If (Maxval(Abs(Aimag(mat2))).Eq.0._dp) Then 
Call EigenSystem(Real(mat2,dp),MFe2,ZER1,ierr,test) 
ZER2 = ZER1 
Else 
Call EigenSystem(mat2,MFe2,ZER2,ierr,test) 
 End If 
 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
  End If 
  ierr = 0 
End If 
 
mat2 = Matmul(mat,Transpose(Conjg(mat))) 
If (Maxval(Abs(Aimag(mat2))).Eq.0._dp) Then 
Call EigenSystem (Real(mat2,dp),MFe2,ZEL1,ierr,test) 
                  
                  
ZEL2 = ZEL1 
Else 
Call EigenSystem(mat2,MFe2,ZEL2,ierr,test) 
 
 
End If 
ZEL2 = Conjg(ZEL2) 
 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
  End If 
  ierr = 0 
End If 
 
mat2 = Matmul(Matmul( Conjg(ZEL2),mat),Transpose( Conjg(ZER2))) 
Do i1=1,3
If (Abs(mat2(i1,i1)).gt.0._dp) Then 
phaseM = mat2(i1,i1) / Abs(mat2(i1,i1)) 
ZER2(i1,:) = phaseM *ZER2(i1,:) 
 End if 
End Do 
 
Do i1=1,3
If (Abs(ZER2(i1,i1)).gt.0._dp) Then 
phaseM = ZER2(i1,i1) / Abs(ZER2(i1,i1)) 
ZER2(i1,:) = Conjg(phaseM) *ZER2(i1,:) 
 ZEL2(i1,:) = phaseM *ZEL2(i1,:) 
 End if 
  If (MFe2(i1).ne.MFe2(i1)) Then 
      Write(*,*) 'NaN appearing in '//NameOfUnit(Iname) 
      Call TerminateProgram 
    End If 
  If (Abs(MFe2(i1)).Le.MaxMassNumericalZero) MFe2(i1) = Abs(MFe2(i1))+1.E-10_dp 
  If (MFe2(i1).Le.0._dp) Then 
! kont = -104 
 End if 
End Do 
 
If (ierr.Ne.0.) Then 
  Write(ErrCan,*) 'Warning from Subroutine CalculateMFe, ierr =',ierr 
  kont = ierr 
  Iname = Iname - 1 
  Return 
End If 


MFe = Sqrt( MFe2 ) 
ZEL = ZEL2 
ZER = ZER2 
Iname = Iname - 1 
 
End Subroutine CalculateMFeEffPot 

Subroutine CalculateMFvEffPot(UV,MFv,kont)

Integer, Intent(inout) :: kont 
Integer :: i1,i2,i3,i4, ierr, pos 
Integer :: j1,j2,j3,j4 
Logical :: SecondDiagonalisationNeeded 
Real(dp), Intent(out) :: MFv(5) 
Complex(dp), Intent(out) ::  UV(5,5) 
                              
Complex(dp) :: mat(5,5), mat2(5,5), phaseM, E5(5) 

Real(dp) :: UVa(5,5), test(2), eig(5) 

Iname = Iname + 1 
NameOfUnit(Iname) = 'CalculateMFv'
 
mat(1,1) = 0._dp 
mat(1,2) = 0._dp 
mat(1,3) = 0._dp 
mat(1,4) = 0._dp 
mat(1,5) = 0._dp 
mat(2,2) = 0._dp 
mat(2,3) = 0._dp 
mat(2,4) = 0._dp 
mat(2,5) = 0._dp 
mat(3,3) = 0._dp 
mat(3,4) = 0._dp 
mat(3,5) = 0._dp 
mat(4,4) = 0._dp 
mat(4,5) = 0._dp 
mat(5,5) = 0._dp 

 
 Do i1=2,5
  Do i2 = 1, i1-1 
  mat(i1,i2) = mat(i2,i1) 
  End do 
End do 

 
If (Maxval(Abs(Aimag(mat))).Eq.0._dp) Then 
Call EigenSystem(Real(mat,dp),Eig,UVa,ierr,test) 
 
   Do i1=1,5
   If ((Eig(i1).Lt.0._dp).or.(Abs(eig(i1)).lt.1E-15)) Then 
    MFv(i1) = - Eig(i1) 
    UV(i1,:) = (0._dp,1._dp)*UVa(i1,:) 
   Else 
    MFv(i1) = Eig(i1) 
    UV(i1,:) = UVa(i1,:)
    End If 
   End Do 
 
  Do i1=1,5
   pos=Maxloc(Abs(UV(i1,:)),1) 
   If (Abs(Real(UV(i1,pos),dp)).gt.Abs(Aimag(UV(i1,pos)))) Then 
      If (Real(UV(i1,pos),dp).lt.0._dp) Then 
        UV(i1,:)=-UV(i1,:) 
       End If 
    Else 
      If (Aimag(UV(i1,pos)).lt.0._dp) Then 
        UV(i1,:)=-UV(i1,:) 
      End If 
    End If 
 End Do 
 
Do i1=1,4
  Do i2=i1+1,5
    If (MFv(i1).Gt.MFv(i2)) Then 
      Eig(1) = MFv(i1) 
      MFv(i1) = MFv(i2) 
      MFv(i2) =  Eig(1) 
      E5 = UV(i1,:) 
      UV(i1,:) = UV(i2,:) 
      UV(i2,:) = E5
    End If 
   End Do 
End Do 
 
Else 
 
mat2 = Matmul( Transpose(Conjg( mat) ), mat ) 
Call Eigensystem(mat2, Eig, UV, ierr, test) 
mat2 = Matmul( Conjg(UV), Matmul( mat, Transpose( Conjg(UV)))) 
! Special efforts are needed for matrices like the Higgsinos one 
SecondDiagonalisationNeeded = .False. 
Do i1=1,5-1
If (MaxVal(Abs(mat2(i1,(i1+1):5))).gt.Abs(mat2(i1,i1))) SecondDiagonalisationNeeded = .True. 

  If (Eig(i1).ne.Eig(i1)) Then 
      Write(*,*) 'NaN appearing in '//NameOfUnit(Iname) 
      Call TerminateProgram 
    End If 
  If ((Abs(Eig(i1)).Le.MaxMassNumericalZero).and.(Eig(i1).lt.0._dp)) Eig(i1) = Abs(Eig(i1))+1.E-10_dp 
  If (Eig(i1).Le.0._dp) Then 
! kont = -104 
 End if 
End do 
If (SecondDiagonalisationNeeded) Then 
Call EigenSystem(Real(mat2,dp),Eig,UVa,ierr,test) 
 
     UV = MatMul(UV,UVa)
  Do i1=1,5
   If ((Eig(i1).Lt.0._dp).or.(Abs(eig(i1)).lt.1E-15)) Then 
    MFv(i1) = - Eig(i1) 
    UV(i1,:) = (0._dp,1._dp)*UVa(i1,:) 
   Else 
    MFv(i1) = Eig(i1) 
    UV(i1,:) = UVa(i1,:)
    End If 
   End Do 
 
Else 
Do i1=1,5
  If (Eig(i1).ne.Eig(i1)) Then 
      Write(*,*) 'NaN appearing in '//NameOfUnit(Iname) 
      Call TerminateProgram 
    End If 
If (Abs(mat2(i1,i1)).gt.0._dp) Then 
  phaseM = Sqrt( mat2(i1,i1) / Abs(mat2(i1,i1))) 
  UV(i1,:)= phaseM * UV(i1,:) 
End if 
  If ((Abs(Eig(i1)).Le.MaxMassNumericalZero).and.(Eig(i1).lt.0._dp)) Eig(i1) = Abs(Eig(i1))+1.E-10_dp 
  If (Eig(i1).Le.0._dp) Then 
! kont = -104 
 End if 
End Do 
MFv = Sqrt( Eig ) 
 
End if ! Second diagonalisation 
End If 
 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
  End If 
  ierr = 0 
End If 
 
If (ierr.Ne.0.) Then 
  Write(ErrCan,*) 'Warning from Subroutine CalculateMFv, ierr =',ierr 
  kont = ierr 
  Iname = Iname - 1 
  Return 
End If 


Iname = Iname - 1 
 
End Subroutine CalculateMFvEffPot 

Subroutine CalculateMFeDEffPot(YB1,YC1,YRH,YRG,v,vX,ZELD,ZERD,MFeD,kont)

Real(dp),Intent(in) :: YB1,YC1,YRH,v,vX

Complex(dp),Intent(in) :: YRG

Integer, Intent(inout) :: kont 
Integer :: i1,i2,i3,i4, ierr 
Integer :: j1,j2,j3,j4 
Real(dp), Intent(out) :: MFeD(2) 
 Complex(dp), Intent(out) :: ZELD(2,2), ZERD(2,2) 
 
 Complex(dp) :: mat(2,2)=0._dp, mat2(2,2)=0._dp, phaseM 

Complex(dp) :: ZELD2(2,2), ZERD2(2,2) 
 
 Real(dp) :: ZELD1(2,2), ZERD1(2,2), test(2), MFeD2(2) 
 
 Iname = Iname + 1 
NameOfUnit(Iname) = 'CalculateMFeD'
 
MFeD = 0._dp 
ZELD = 0._dp 
ZERD = 0._dp 
mat(1,1) = 0._dp 
mat(1,1) = mat(1,1)+(v*YRG)/sqrt(2._dp)
mat(1,2) = 0._dp 
mat(1,2) = mat(1,2)+(vX*YC1)/sqrt(2._dp)
mat(2,1) = 0._dp 
mat(2,1) = mat(2,1)+(vX*YB1)/sqrt(2._dp)
mat(2,2) = 0._dp 
mat(2,2) = mat(2,2)+(v*YRH)/sqrt(2._dp)

 
mat2 = Matmul(Transpose(Conjg(mat)),mat) 
If (Maxval(Abs(Aimag(mat2))).Eq.0._dp) Then 
Call EigenSystem(Real(mat2,dp),MFeD2,ZERD1,ierr,test) 
ZERD2 = ZERD1 
Else 
Call EigenSystem(mat2,MFeD2,ZERD2,ierr,test) 
 End If 
 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
  End If 
  ierr = 0 
End If 
 
mat2 = Matmul(mat,Transpose(Conjg(mat))) 
If (Maxval(Abs(Aimag(mat2))).Eq.0._dp) Then 
Call EigenSystem (Real(mat2,dp),MFeD2,ZELD1,ierr,test) 
                  
                  
ZELD2 = ZELD1 
Else 
Call EigenSystem(mat2,MFeD2,ZELD2,ierr,test) 
 
 
End If 
ZELD2 = Conjg(ZELD2) 
 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
  End If 
  ierr = 0 
End If 
 
mat2 = Matmul(Matmul( Conjg(ZELD2),mat),Transpose( Conjg(ZERD2))) 
Do i1=1,2
If (Abs(mat2(i1,i1)).gt.0._dp) Then 
phaseM = mat2(i1,i1) / Abs(mat2(i1,i1)) 
ZERD2(i1,:) = phaseM *ZERD2(i1,:) 
 End if 
End Do 
 
Do i1=1,2
If (Abs(ZERD2(i1,i1)).gt.0._dp) Then 
phaseM = ZERD2(i1,i1) / Abs(ZERD2(i1,i1)) 
ZERD2(i1,:) = Conjg(phaseM) *ZERD2(i1,:) 
 ZELD2(i1,:) = phaseM *ZELD2(i1,:) 
 End if 
  If (MFeD2(i1).ne.MFeD2(i1)) Then 
      Write(*,*) 'NaN appearing in '//NameOfUnit(Iname) 
      Call TerminateProgram 
    End If 
  If (Abs(MFeD2(i1)).Le.MaxMassNumericalZero) MFeD2(i1) = Abs(MFeD2(i1))+1.E-10_dp 
  If (MFeD2(i1).Le.0._dp) Then 
! kont = -104 
 End if 
End Do 
 
If (ierr.Ne.0.) Then 
  Write(ErrCan,*) 'Warning from Subroutine CalculateMFeD, ierr =',ierr 
  kont = ierr 
  Iname = Iname - 1 
  Return 
End If 


MFeD = Sqrt( MFeD2 ) 
ZELD = ZELD2 
ZERD = ZERD2 
Iname = Iname - 1 
 
End Subroutine CalculateMFeDEffPot 

Subroutine CalculateVPVZVZpEffPot(g1,g2,g1p,gBY,gYB,v,vX,ZZ,MVZ,MVZp,MVZ2,            & 
& MVZp2,TW,TWp,kont)

Real(dp), Intent(in) :: g1,g2,g1p,gBY,gYB,v,vX

Real(dp), Intent(out) :: TW,TWp

Integer, Intent(inout) :: kont 
Integer :: i1,i2,i3,i4, ierr 
Integer :: j1,j2,j3,j4 
Real(dp), Intent(out) :: MVZ, MVZ2
Real(dp), Intent(out) :: MVZp, MVZp2
Real(dp) :: VPVZVZp2(3),VPVZVZp(3)  

Real(dp), Intent(out) :: ZZ(3,3) 
 
Real(dp) :: mat(3,3)  

Real(dp) ::  test(2) 

Iname = Iname + 1 
NameOfUnit(Iname) = 'CalculateVPVZVZp'
 
mat(1,1) = 0._dp 
mat(1,1) = mat(1,1)+(g1**2*vFix**2)/4._dp
mat(1,1) = mat(1,1)+25*gYB**2*vXFix**2
mat(1,2) = 0._dp 
mat(1,2) = mat(1,2)-(g1*g2*vFix**2)/4._dp
mat(1,3) = 0._dp 
mat(1,3) = mat(1,3)+(g1*gBY*vFix**2)/4._dp
mat(1,3) = mat(1,3)+25*g1p*gYB*vXFix**2
mat(2,2) = 0._dp 
mat(2,2) = mat(2,2)+(g2**2*vFix**2)/4._dp
mat(2,3) = 0._dp 
mat(2,3) = mat(2,3)-(gBY*g2*vFix**2)/4._dp
mat(3,3) = 0._dp 
mat(3,3) = mat(3,3)+(gBY**2*vFix**2)/4._dp
mat(3,3) = mat(3,3)+25*g1p**2*vXFix**2

 
 Do i1=2,3
  Do i2 = 1, i1-1 
  mat(i1,i2) = mat(i2,i1) 
  End do 
End do 

 
Call EigenSystem(mat,VPVZVZp2,ZZ,ierr,test) 
 
 
ZZ = Transpose(ZZ) 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
  End If 
  ierr = 0 
End If 
 
If ((ierr.Ne.0.).And.(ErrorLevel.Ge.-1)) Then 
  Write(10,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
  Write(10,*) 'Diagonalization failed, ierr : ',ierr 
  kont = ierr 
  Iname = Iname - 1 
  Return 
End If 


Do i1=1,3
  If (Abs(VPVZVZp2(i1)).Le.1.E-10_dp*(Maxval(VPVZVZp2))) VPVZVZp2(i1) = 1.E-10_dp 
  If (VPVZVZp2(i1).ne.VPVZVZp2(i1)) Then 
      Write(*,*) 'NaN appearing in '//NameOfUnit(Iname) 
      Call TerminateProgram 
    End If 
  If (VPVZVZp2(i1).Ge.0._dp) Then 
  VPVZVZp(i1) =Sqrt(VPVZVZp2(i1) ) 
  Else 
  VPVZVZp(i1)= 1._dp 
  VPVZVZp(i1)= 1._dp 
! kont = -104 
 End if 
End Do 
 
MVZ = VPVZVZp(2) 
MVZ2 = VPVZVZp2(2) 
MVZp = VPVZVZp(3) 
MVZp2 = VPVZVZp2(3) 

 Iname = Iname - 1 
 
End Subroutine CalculateVPVZVZpEffPot 

Subroutine CalculateVWmEffPot(g2,v,ZW,MVWm,MVWm2,kont)

Real(dp), Intent(in) :: g2,v

Integer, Intent(inout) :: kont 
Integer :: i1,i2,i3,i4, ierr 
Integer :: j1,j2,j3,j4 
Real(dp), Intent(out) :: MVWm, MVWm2
Real(dp) :: VWm2(2),VWm(2)  

Complex(dp), Intent(out) :: ZW(2,2) 
 
Complex(dp) :: mat(2,2)  

Real(dp) ::  test(2) 

Iname = Iname + 1 
NameOfUnit(Iname) = 'CalculateVWm'
 
mat(1,1) = 0._dp 
mat(1,1) = mat(1,1)+(g2**2*vFix**2)/4._dp
mat(1,2) = 0._dp 
mat(2,2) = 0._dp 
mat(2,2) = mat(2,2)+(g2**2*vFix**2)/4._dp

 
 Do i1=2,2
  Do i2 = 1, i1-1 
  mat(i1,i2) = Conjg(mat(i2,i1)) 
  End do 
End do 

 
Call EigenSystem(mat,VWm2,ZW,ierr,test) 
 
 
ZW = Transpose(ZW) 
If ((ierr.Eq.-8).Or.(ierr.Eq.-9)) Then 
  Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname) 
  If (ErrorLevel.Eq.2) Then 
  Write(*,*) "Possible numerical problem in "//NameOfUnit(Iname) 
    Call TerminateProgram 
  End If 
  ierr = 0 
End If 
 
If ((ierr.Ne.0.).And.(ErrorLevel.Ge.-1)) Then 
  Write(10,*) 'Warning from Subroutine '//NameOfUnit(Iname) 
  Write(10,*) 'Diagonalization failed, ierr : ',ierr 
  kont = ierr 
  Iname = Iname - 1 
  Return 
End If 


Do i1=1,2
  If (Abs(VWm2(i1)).Le.1.E-10_dp*(Maxval(VWm2))) VWm2(i1) = 1.E-10_dp 
  If (VWm2(i1).ne.VWm2(i1)) Then 
      Write(*,*) 'NaN appearing in '//NameOfUnit(Iname) 
      Call TerminateProgram 
    End If 
  If (VWm2(i1).Ge.0._dp) Then 
  VWm(i1) =Sqrt(VWm2(i1) ) 
  Else 
  VWm(i1)= 1._dp 
  VWm(i1)= 1._dp 
! kont = -104 
 End if 
End Do 
 
MVWm = VWm(1) 
MVWm2 = VWm2(1) 

 Iname = Iname - 1 
 
End Subroutine CalculateVWmEffPot 

Subroutine TreeMassesSM(MFd,MFd2,MFe,MFe2,MFu,MFu2,MFv,MFv2,MVWm,MVWm2,               & 
& MVZ,MVZ2,MVZp,MVZp2,TW,TWp,ZDR,ZER,ZUR,ZDL,ZEL,ZUL,UV,ZW,ZZ,v,vX,g1,gBY,               & 
& g2,g3,g1p,gYB,lamh,lam2,lam3,lam4,lam5,lam6,lam7,Yu,Yd,Ye,YB1,YC1,YA1,YB2,             & 
& YC2,YA2,YRH,LE,YRG,LD,Fi,Fj,muh,mui,muj,m1,m2,GenerationMixing,kont)

Implicit None 
 
Real(dp),Intent(in) :: g1,gBY,g2,g3,g1p,gYB,YB1,YC1,YA1,YB2,YC2,YA2,YRH,Fi(2,2),Fj(2,2)

Complex(dp),Intent(in) :: lamh,lam2,lam3,lam4,lam5,lam6(2,2),lam7(2,2),Yu(3,3),Yd(3,3),Ye(3,3),LE(2,2),         & 
& YRG,LD(2,3),muh,mui,muj,m1(2,2),m2(2,2)

Real(dp),Intent(out) :: MFd(3),MFd2(3),MFe(3),MFe2(3),MFu(3),MFu2(3),MFv(5),MFv2(5),MVWm,MVWm2,               & 
& MVZ,MVZ2,MVZp,MVZp2,TW,TWp,ZZ(3,3)

Complex(dp),Intent(out) :: ZDR(3,3),ZER(3,3),ZUR(3,3),ZDL(3,3),ZEL(3,3),ZUL(3,3),UV(5,5),ZW(2,2)

Real(dp),Intent(in) :: v,vX

Logical, Intent(in) :: GenerationMixing 
Integer, Intent(inout) :: kont 
Integer :: i1,i2,i3,i4,j1,j2,j3,kontSave 
Iname = Iname + 1 
NameOfUnit(Iname) = 'TreeMassesU1B'
 
kont = 0 
Call CalculateMFd(Yd,v,ZDL,ZDR,MFd,kont)

MFd2 = MFd**2 
Call CalculateMFu(Yu,v,ZUL,ZUR,MFu,kont)

MFu2 = MFu**2 
Call CalculateMFe(Ye,v,ZEL,ZER,MFe,kont)

MFe2 = MFe**2 
Call CalculateMFv(UV,MFv,kont)

MFv2 = MFv**2 

 
 Call CalculateVPVZVZp(g1,g2,g1p,gBY,gYB,v,vX,ZZ,MVZ,MVZp,MVZ2,MVZp2,TW,TWp,kont)

Call CalculateVWm(g2,v,ZW,MVWm,MVWm2,kont)

Iname = Iname - 1 
 
End Subroutine  TreeMassesSM 
 
 
Subroutine SortGoldstones(MAh,MAh2,MFd,MFd2,MFe,MFe2,MFeD,MFeD2,MFu,MFu2,             & 
& MFv,MFv2,MFvv,MFvv2,MFx,MFx2,Mhh,Mhh2,MHm,MHm2,Mss,Mss2,MVWm,MVWm2,MVZ,MVZ2,           & 
& MVZp,MVZp2,TW,TWp,ZDR,ZER,ZERD,ZUR,ZDL,ZEL,ZELD,ZUL,ZA,ZH,UV,ZP,ZW,ZZ,kont)

Real(dp),Intent(inout) :: MAh(2),MAh2(2),MFd(3),MFd2(3),MFe(3),MFe2(3),MFeD(2),MFeD2(2),MFu(3),MFu2(3),         & 
& MFv(5),MFv2(5),MFvv,MFvv2,MFx,MFx2,Mhh(2),Mhh2(2),MHm(5),MHm2(5),Mss,Mss2,             & 
& MVWm,MVWm2,MVZ,MVZ2,MVZp,MVZp2,TW,TWp,ZA(2,2),ZH(2,2),ZP(5,5),ZZ(3,3)

Complex(dp),Intent(inout) :: ZDR(3,3),ZER(3,3),ZERD(2,2),ZUR(3,3),ZDL(3,3),ZEL(3,3),ZELD(2,2),ZUL(3,3),            & 
& UV(5,5),ZW(2,2)

Integer, Intent(inout) :: kont 
Integer :: i1, i2, pos 
Real(dp) :: MAhtemp(2) 
Complex(dp) :: ZAhtemp(2,2) 
Real(dp) :: MHmtemp(5) 
Complex(dp) :: ZHmtemp(5,5) 


pos = MinLoc(Abs(MAh2-MVZ2*RXiZ),1) 
If (pos.ne.1) Then 
  MAhtemp = MAh2 
  ZAhtemp = ZA 
  MAh2(1) = MAhtemp(pos) 
  ZA(1,:) = ZAhtemp(pos,:) 
  MAh2(pos) = MAhtemp(1) 
  ZA(pos,:) = ZAhtemp(1,:) 
End if 
pos = MinLoc(Abs(MAh2-MVZp2*RXiZp),1) 
If (pos.ne.2) Then 
  MAhtemp = MAh2 
  ZAhtemp = ZA 
  MAh2(2) = MAhtemp(pos) 
  ZA(2,:) = ZAhtemp(pos,:) 
  MAh2(pos) = MAhtemp(2) 
  ZA(pos,:) = ZAhtemp(2,:) 
End if 

 ! Reorder the physical states by their mass 
Do i1=3,2
 pos = Minloc(MAh2(i1:2),1) + i1 -1  
If (pos.ne.i1) Then 
  MAhtemp = MAh2 
  ZAhtemp = ZA 
  MAh2(i1) = MAhtemp(pos) 
  ZA(i1,:) = ZAhtemp(pos,:) 
  MAh2(pos) = MAhtemp(i1) 
  ZA(pos,:) = ZAhtemp(i1,:) 
End if 
End do 
MAh = sqrt(MAh2) 

 
 
pos = MinLoc(Abs(MHm2-MVWm2*RXiWm),1) 
If (pos.ne.1) Then 
  MHmtemp = MHm2 
  ZHmtemp = ZP 
  MHm2(1) = MHmtemp(pos) 
  ZP(1,:) = ZHmtemp(pos,:) 
  MHm2(pos) = MHmtemp(1) 
  ZP(pos,:) = ZHmtemp(1,:) 
End if 

 ! Reorder the physical states by their mass 
Do i1=2,5
 pos = Minloc(MHm2(i1:5),1) + i1 -1  
If (pos.ne.i1) Then 
  MHmtemp = MHm2 
  ZHmtemp = ZP 
  MHm2(i1) = MHmtemp(pos) 
  ZP(i1,:) = ZHmtemp(pos,:) 
  MHm2(pos) = MHmtemp(i1) 
  ZP(pos,:) = ZHmtemp(i1,:) 
End if 
End do 
MHm = sqrt(MHm2) 

 
 
End subroutine SortGoldstones 


End Module TreeLevelMasses_U1B 
 
