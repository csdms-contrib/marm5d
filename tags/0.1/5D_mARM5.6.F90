! The mARM5D model - Simulating 3D soil evolution. 
!Copyright (C) 2014 Sagy Cohen
!Developer can be contacted by sagy.cohen@ua.edu and Box 870322, Tuscaloosa, AL 35487-0322 U.S.A.
!This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
!This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA. 

! First version: 1/7/2008
! Last modified: 18/6/10 (5D version 5.6)

 PROGRAM mARM4D
 USE IFPORT
 IMPLICIT NONE

    Integer NumGrade,NoNodes,MaxGrads,iii,ig,i,ii,j,n,xii,ix,unit1,unit2,NCols,NRows,r,c,Numout,dt,Numlayers,l,temp,Mode,NumIter,GT2grade,count,tmpj,counter,indx,cj,rj,dir
    Integer to_r,to_c ,ND(8,2),cr,cc,rr,NumOutputs,OUTcount,WeatherFront,ChangeCount
    Character*80 D50aFile,TEFile,tmp,FAccName,SlopeName ,CurvName,DirName,DFile,FEFile,InputName,ClimateName,GradingName,AeolGradName
    Real(8) ErosionSum,NetErosion,CGSum,Etmp,a,TmpErosion,SplitProp,WeatherAlpha,LayerDepth,Erosion,ArmourDepth,tmpRock,d50,tmpp,DepoProp,SurfLayerVol,E
    Real(8) xllcor, yllcor,CellSize,pre,tNRows,FracSum, TempGrad,rn,Qf,cd50,TmpBR,DepthWeatheringAlpha,StepWeatheringChange, StepErosionChange, StepDischargeChange
    Real(8) DepthWeatheringAlpha2,tmpSlope,tmpQf,kErode,tmpkErode, tmpWeatheringChange, tmpErosionChange, tmpDischargeChange,zero,CreepRate,tmpDeficit 
    Real(8) RandChange,WeatherChange,UpSlopTotal,tmpPC,OutPC,Qs,b,AeolianRate,PushProp,tmpBRprop,tmpAeolianChange,StepAeolianChange,tmpAeolianRate
	Real(8) tmpCreepChange,StepCreepChange, tmpCreepRate 
    REAL*4 ran2
    Real(8),ALLOCATABLE :: IterPrecent(:), WeatheringChange(:), ErosionChange(:),Frac2mm(:),Dmin(:,:),tmpUSGrad(:),AeolDepGrad(:)
    Real(8),ALLOCATABLE :: InitGrad(:), ArmourTrans(:), GradingSize(:),MeanGrading(:),LayerWeatheringM(:,:,:),tmpLayerWeatheringM(:,:,:)
    Real(8),ALLOCATABLE :: NextGrad(:),ResupplyGrad(:),ErosionGrad(:),ErodedGrad(:),tmpGradSize(:),TempSurfGrad(:),TempLayerGrad(:,:),LayerWeatherAlpha(:)
    Real(8),ALLOCATABLE :: NodeTransition(:), WeatheringTrans(:,:), ProfileGrading(:,:,:,:),LayerD50(:,:),Depth(:,:),TempSurfGrad2(:),TempLayerGrad2(:,:)
    Real(8),ALLOCATABLE :: PureErosion(:),FinalErosion(:,:),FinalNetErosion(:,:),DischargeChange(:),tmpLayerWeatherAlpha(:),AeolianChange(:),CreepChange(:)
    Real(8),ALLOCATABLE :: FinSurfGrad(:,:),FlowAcc(:,:),FlowDir(:,:),BRLayerWeatherAlpha(:),tmpBRLayerWeatherAlpha(:),OutPCv(:),DepositedGrad(:),PotenDepo(:)
    Real(8),ALLOCATABLE :: TotalErosion(:,:), Slope(:,:), Qd(:,:),SizeTh(:,:),Curvature(:),PreD50(:),d50map(:,:),Junc(:,:),tmpJuncsions(:,:)
    Real(8),ALLOCATABLE :: SurfaceGrad(:,:,:), UpSlopFlowGrad(:,:,:),QsUS(:,:),DepositionMatrix(:),CreepDepth(:),SurfCreepProp(:,:),SurfCreepRate(:,:)
    INTEGER(4) tmpday, tmpmonth, tmpyear, ClimatFluc, ClimateRec
    INTEGER(4) tmphour, tmpminute, tmpsecond, tmphund, CurrentStep,seed1
    INTEGER,ALLOCATABLE ::  Junc0(:,:),tmpJunc0(:,:),OutputTS(:),InputTS(:)
    CHARACTER(1) mer
    DATA ND  /1,1,0,-1,-1,-1,0,1,0,-1,-1,-1,0,1,1,1/  
    DATA  seed1 / 1 /
            
    WRITE (*,*) 'mARM5D(c) - Simulating soil evolution as a function of weathering and' 
    WRITE (*,*) 'armouring processes'
    WRITE (*,*) 'Written by Sagy Cohen in cooperation with '
    WRITE (*,*) 'Prof. Garry Willgoose and Dr Greg Hancock'
    WRITE (*,*) 'The University of Newcastle, Australia'
    WRITE (*,*) 'Version 5.6- updated 27/5/2010'
    WRITE (*,*) '==================================================================='
    WRITE (*,*) 'Input layers (contributing area, slope and flow direction) must be in'
    WRITE (*,*) 'ArcGIS raster ASCII format.'
    WRITE (*,*) '==================================================================='
    WRITE (*,*) ''
    
    WRITE (*,*) 'What is the input data file name?'
    Read (*,*)  InputName !number of Nodes
    InputName='inputparam.dat'
1   Format(a,a)
2   Format(a30,5F10.5)
3   Format(a30,I0) 
   
   ! Read Parameters file 
    Open (unit=1,file=InputName,status='old')
    Read (1,"(a27,a)") tmp,FAccName
    Read (1,"(a16,a)") tmp,SlopeName
    Read (1,"(a18,a)") tmp,DirName
    Read (1,"(a21,I)") tmp, NumIter
    Read (1,"(a19,I)") tmp,dt
    Read (1,"(a22,5F10.5)") tmp,Qf
    Read (1,"(a23,5F10.5)") tmp,WeatherAlpha
    Read (1,"(a25,5F10.5)") tmp,AeolianRate
    Read (1,"(a25,5F10.5)") tmp,CreepRate 
    Read (1,"(a25,I)") tmp,Numlayers
    Read (1,"(a32,5F10.5)") tmp,LayerDepth
    Read (1,"(a32,5F10.5)") tmp,ArmourDepth
    Read (1,"(a22,5F10.5)") tmp,kErode
    Read (1,"(a35,5F10.5)") tmp,SplitProp
    Read (1,"(a80)") tmp
    Read (1,"(a75,I)") tmp,Mode
    Read (1,"(a44,I)") tmp,ClimatFluc
    If (ClimatFluc.EQ.1) Read (1,"(a31,a)") tmp,ClimateName
    Read (1,"(a23,a)") tmp,GradingName
    Read (1,"(a27,a)") tmp,AeolGradName
    CALL GETDAT(tmpyear, tmpmonth, tmpday)
    CALL GETTIM(tmphour, tmpminute, tmpsecond, tmphund)
22 Format(a,I0,a,I0,a,I0,a,I0,a,I0,a,I0)
    Write (1,22) 'Start: ',tmpday,'/',tmpmonth,'/',tmpyear,' ',tmphour,':',tmpminute,':',tmpsecond

    WeatherAlpha=WeatherAlpha*dt
    AeolianRate=(AeolianRate*dt)
    Open (unit=9,file=FAccName,status='old')
    Read (9,*) tmp,NCols
    Read (9,*) tmp,NRows
    Read (9,*) tmp,xllcor
    Read (9,*) tmp,yllcor
    Read (9,*) tmp,CellSize
    Read (9,*) tmp

    Open (unit=8,file=SlopeName,status='old')
    Read (8,*) tmp,NCols
    Read (8,*) tmp,NRows
    Read (8,*) tmp,xllcor
    Read (8,*) tmp,yllcor
    Read (8,*) tmp,CellSize
    Read (8,*) tmp
    
    Open (unit=6,file=DirName,status='old')
    Read (6,*) tmp,NCols
    Read (6,*) tmp,NRows
    Read (6,*) tmp,xllcor
    Read (6,*) tmp,yllcor
    Read (6,*) tmp,CellSize
    Read (6,*) tmp    
    Allocate (Junc0(NCols*NRows,2),tmpJunc0(NCols*NRows,2))
    Allocate (FlowAcc(NCols,NRows),Slope(NCols,NRows),d50map(NCols,NRows),FlowDir(NCols,NRows),Junc(NCols,NRows),tmpJuncsions(NCols,NRows))
    Do i=1,NRows
        Read (9,*) FlowAcc(:,i) ! Read the flow accumulation values
        Read (8,*) Slope(:,i)
        Read (6,*) FlowDir(:,i) ! Read the slope values 
    End Do 
     
    Open (unit=10,file=GradingName,status='old')
    Open (unit=11,file=AeolGradName,status='old')
    Read (10,*) NumGrade
    Read (11,*) NumGrade
! Read the recoreds of the climate fluctuation file
    If (ClimatFluc.NE.0) Then
        Open (unit=890,file=ClimateName,status='old')
        Read (890,*) tmp
        Read (890,*) ClimateRec, NumOutputs
        Read (890,*) tmp
        Allocate (IterPrecent(ClimateRec), WeatheringChange(ClimateRec), ErosionChange(ClimateRec), DischargeChange(ClimateRec),OutputTS(NumOutputs),InputTS(ClimateRec))
!for 1D! Allocate (CCd50(NumOutputs,NCols),CCdepth(NumOutputs,NCols),CCNetE(NumOutputs,NCols),CCFrac2mm(NumOutputs,NCols),Frac2mm(NCols))
!for 1D! Allocate (formatCCd50(NumOutputs,NRows+1),formatCCdepth(NumOutputs,NRows+1),formatCCNetE(NumOutputs,NRows+1),formatCCFrac2mm(NumOutputs,NRows+1))
        Allocate (OutPCv(NumOutputs),AeolianChange(ClimateRec),CreepChange(ClimateRec))
        Do i=1, ClimateRec
            Read (890,*) IterPrecent(i), WeatheringChange(i), ErosionChange(i), DischargeChange(i), AeolianChange(i), CreepChange(i)
        End Do
    End If

 ! Calculate output and climate input time steps
    tmpPC = (100/(Float(NumOutputs)))      
    Do  i=1,ClimateRec
        InputTS(i)=IterPrecent(i)*(NumIter/100)
    End Do
    Do  i=1, NumOutputs
        OutputTS(i)=(tmpPC/100)*(i-1)*NumIter
    End Do       

! Climate change output files
995  Format (a,F5.2,a)      
    If (ClimatFluc.NE.0) Then   
        Do i=1, NumOutputs
            OutPC=REAL(OutputTS(i))/REAL(NumIter)*100
            write (TEFile,995) 'd50aL0',OutPC,'pc.txt'    
            OPEN(unit=(i+590),file=TEFile,status='unknown')
            Write ((i+590),*) 'ncols ',NCols
            Write ((i+590),*) 'nrows ',NRows
            Write ((i+590),*) 'xllcorner ',xllcor 
            Write ((i+590),*) 'yllcorner ',yllcor 
            Write ((i+590),*) 'cellsize  ',CellSize 
            Write ((i+590),*) 'NODATA_value  -9999' 
            
            write (TEFile,995) 'Depth',OutPC,'pc.txt'    
            OPEN(unit=(i+690),file=TEFile,status='unknown')
            Write ((i+690),*) 'ncols ',NCols
            Write ((i+690),*) 'nrows ',NRows
            Write ((i+690),*) 'xllcorner ',xllcor 
            Write ((i+690),*) 'yllcorner ',yllcor 
            Write ((i+690),*) 'cellsize  ',CellSize 
            Write ((i+690),*) 'NODATA_value  -9999' 
            
!            write (TEFile,995) 'FracGT2mm',OutPC,'pc.txt'    
!            OPEN(unit=(i+1090),file=TEFile,status='unknown')
!            Write ((i+1090),*) 'ncols ',NCols
!            Write ((i+1090),*) 'nrows ',NRows
!            Write ((i+1090),*) 'xllcorner ',xllcor 
!            Write ((i+1090),*) 'yllcorner ',yllcor 
!            Write ((i+1090),*) 'cellsize  ',CellSize 
!            Write ((i+1090),*) 'NODATA_value  -9999' 
            
            write (TEFile,995) 'FinalNetE',OutPC,'pc.txt'    
            OPEN(unit=(i+790),file=TEFile,status='unknown')
            Write ((i+790),*) 'ncols ',NCols
            Write ((i+790),*) 'nrows ',NRows
            Write ((i+790),*) 'xllcorner ',xllcor 
            Write ((i+790),*) 'yllcorner ',yllcor 
            Write ((i+790),*) 'cellsize  ',CellSize 
            Write ((i+790),*) 'NODATA_value  -9999' 
        End Do
    End If

!write final surface grading files
1500 Format (a,I2,a) 
    Do i=1,NumGrade
        write (TEFile,1500) 'FinSurfGrad',i,'.txt'
        OPEN(unit=(i+1500),file=TEFile,status='unknown')
        Write ((i+1500),*) 'ncols ',NCols
        Write ((i+1500),*) 'nrows ',NRows
        Write ((i+1500),*) 'xllcorner ',xllcor 
        Write ((i+1500),*) 'yllcorner ',yllcor 
        Write ((i+1500),*) 'cellsize  ',CellSize 
        Write ((i+1500),*) 'NODATA_value  -9999' 
    End Do
                  
    Allocate (InitGrad(NumGrade), ArmourTrans(NumGrade), GradingSize(NumGrade), NodeTransition(NumGrade),AeolDepGrad(NumGrade) &
    ,MeanGrading(NumGrade),TempSurfGrad(NumGrade),TempLayerGrad(Numlayers,NumGrade),LayerWeatherAlpha(NumLayers+1),tmpLayerWeatherAlpha(NumLayers+1))
    Allocate (NextGrad(NumGrade),ResupplyGrad(NumGrade),ErosionGrad(NumGrade),ErodedGrad(NumGrade),tmpGradSize(NumGrade)&
    , WeatheringTrans(NumGrade,NumGrade), ProfileGrading(NCols,NRows,Numlayers,NumGrade),LayerD50(NCols,NRows),PreD50(NCols))
    Allocate (LayerWeatheringM(NumLayers+1,NumGrade,NumGrade),tmpLayerWeatheringM(NumLayers+1,NumGrade,NumGrade),PotenDepo(NumGrade),tmpBRLayerWeatherAlpha(NumLayers+1)&
    ,FinSurfGrad(NumGrade,NCols),Dmin(NCols,NRows),BRLayerWeatherAlpha(NumLayers+1),DepositionMatrix(NumGrade),tmpUSGrad(NumGrade),DepositedGrad(NumGrade))
    Allocate (CreepDepth(NumLayers+1),SurfCreepProp(NCols,NRows),SurfCreepRate(NCols,NRows),TempSurfGrad2(NumGrade),TempLayerGrad2(Numlayers,NumGrade)) 
    Do i=1, NumGrade
        Read (10,*) GradingSize(i), InitGrad(i)
		Read (11,*) GradingSize(i), AeolDepGrad(i)
    End Do
   
    Allocate (TotalErosion(NCols,NRows),SurfaceGrad(NCols,NRows,NumGrade),Qd(NCols,NRows),SizeTh(NCols,NRows),Curvature(NCols),Depth(NCols,NRows) &
    ,PureErosion(NCols),FinalErosion(NCols,NRows),FinalNetErosion(NCols,NRows),UpSlopFlowGrad(NCols,NRows,NumGrade),QsUS(NCols,NRows))

    unit2 = 200
 999  Format (a)
    write (TEFile,999) 'TotalErosion.txt'
    OPEN(unit=unit2,file=TEFile,status='unknown')
    Write (unit2,*) 'ncols ',NCols
    Write (unit2,*) 'nrows ',NRows
    Write (unit2,*) 'xllcorner ',xllcor 
    Write (unit2,*) 'yllcorner ',yllcor 
    Write (unit2,*) 'cellsize  ',CellSize 
    Write (unit2,*) 'NODATA_value  -9999'   

    write (DFile,999) 'Depth.txt'
    OPEN(unit=300,file=DFile,status='unknown')
    Write (300,*) 'ncols ',NCols
    Write (300,*) 'nrows ',NRows
    Write (300,*) 'xllcorner ',xllcor 
    Write (300,*) 'yllcorner ',yllcor 
    Write (300,*) 'cellsize  ',CellSize 
    Write (300,*) 'NODATA_value  -9999'   
    
!    write (FEFile,999) 'PureErosion.txt'
!    OPEN(unit=400,file=FEFile,status='unknown')
!    Write (400,*) 'ncols ',NCols
!    Write (400,*) 'nrows ',NRows
!    Write (400,*) 'xllcorner ',xllcor 
!    Write (400,*) 'yllcorner ',yllcor 
!    Write (400,*) 'cellsize  ',CellSize 
!    Write (400,*) 'NODATA_value  -9999'  
    
    write (FEFile,999) 'FinalErosion.txt'
    OPEN(unit=500,file=FEFile,status='unknown')
    Write (500,*) 'ncols ',NCols
    Write (500,*) 'nrows ',NRows
    Write (500,*) 'xllcorner ',xllcor 
    Write (500,*) 'yllcorner ',yllcor 
    Write (500,*) 'cellsize  ',CellSize 
    Write (500,*) 'NODATA_value  -9999'   
    
    write (FEFile,999) 'FinalNetEro.txt'
    OPEN(unit=1600,file=FEFile,status='unknown')
    Write (1600,*) 'ncols ',NCols
    Write (1600,*) 'nrows ',NRows
    Write (1600,*) 'xllcorner ',xllcor 
    Write (1600,*) 'yllcorner ',yllcor 
    Write (1600,*) 'cellsize  ',CellSize 
    Write (1600,*) 'NODATA_value  -9999'  
    
!    write (FEFile,999) 'FracGT2mm.txt'
!    OPEN(unit=1700,file=FEFile,status='unknown')
!    Write (1700,*) 'ncols ',NCols
!    Write (1700,*) 'nrows ',NRows
!    Write (1700,*) 'xllcorner ',xllcor 
!    Write (1700,*) 'yllcorner ',yllcor 
!    Write (1700,*) 'cellsize  ',CellSize 
!    Write (1700,*) 'NODATA_value  -9999'  
    
    write (FEFile,999) 'Junctions.txt'
    OPEN(unit=3210,file=FEFile,status='unknown')
    Write (3210,*) 'ncols ',NCols
    Write (3210,*) 'nrows ',NRows
    Write (3210,*) 'xllcorner ',xllcor 
    Write (3210,*) 'yllcorner ',yllcor 
    Write (3210,*) 'cellsize  ',CellSize 
    Write (3210,*) 'NODATA_value  -9999'  


996  Format (a,I0,a)    
     Do i=0, NumLayers ! Creating layers D50 output files
        unit1 = 100+i
        write (D50aFile,996) 'D50aL',i,'.txt'
        OPEN(unit=unit1,file=D50aFile,status='unknown')
        Write (unit1,*) 'ncols ',NCols
        Write (unit1,*) 'nrows ',NRows
        Write (unit1,*) 'xllcorner ',xllcor 
        Write (unit1,*) 'yllcorner ',yllcor 
        Write (unit1,*) 'cellsize  ',CellSize 
        Write (unit1,*) 'NODATA_value  -9999' 
    End Do
  
 !Write Phyton script file for converting outputs to ArcGIS format 
    If (ClimatFluc.NE.0) then
        OutPCv=REAL(OutputTS)/REAL(NumIter)*100
        Call PythonScriptWriter (NumLayers,NumOutputs,OutPCv(:),ClimatFluc,NumGrade)
    Else
!       Call PythonScriptWriter (NumLayers,ClimateRec,1.0,ClimatFluc,NumGrade) 
    End If
!    
    Do i=1, NumGrade
        if (GradingSize(i).GE.2) Then
        GT2grade=i
        Exit
        End If
    End Do 
    MeanGrading = (GradingSize+CSHIFT(GradingSize,1,1))/2
    MeanGrading(NumGrade)=GradingSize(NumGrade)
    Qd = Qf*(FlowAcc)/CellSize**2 
    SizeTh=((((0.1*Qd)**0.6) *(Slope**0.7))/0.074)*1000 
    !Calculate minimum entrainable sediment size
   
    tmpj=MAXVAL(FlowDir,MASK=FlowDir.LT.200) !find the maximum flow direction value
    If (tmpj<9) Call Junctions_TauDEM(FlowDir,NCols,NRows,Junc) 
    If (tmpj>8)Call Junctions_ArcGIS(FlowDir,NCols,NRows,Junc) 
    d50map=-9999
    Junc0=0
    tmpJuncsions=0
    tmpJunc0=0
    counter=0
    indx=0
    OUTcount=1
    TotalErosion=-9999
    FinalErosion=-9999
    FinalNetErosion=-9999
    Depth=-9999
    tmpWeatheringChange = 0 
    tmpErosionChange = 0
    tmpDischargeChange = 0
    StepWeatheringChange  = 0
    StepErosionChange = 0
    StepDischargeChange = 0
    StepAeolianChange = 0
	StepCreepChange = 0
    tmpWeatheringChange = WeatheringChange(1)
    tmpErosionChange = ErosionChange(1)
    tmpDischargeChange = DischargeChange(1)
    tmpAeolianChange = AeolianChange(1)
	tmpCreepChange = CreepChange(1)
    UpSlopFlowGrad=0

    Do cj=1,NCols
        Do rj=1,NRows
            counter=counter+1
            If (Junc(cj,rj)==0) then
                 indx=indx+1
                 Junc0(indx,1)=cj
                 Junc0(indx,2)=rj
            End if
            Do i=1,NumGrade
               if  (GradingSize(i)<=SizeTh(cj,rj)) Then
                   Dmin(cj,rj)=i;
               end if
            End Do
        End Do
    End Do
    
    Call WeatheringRateMatrixes (Mode,WeatherAlpha,LayerDepth,NumLayers,NumGrade,GradingSize,SplitProp &
    ,LayerWeatheringM,LayerWeatherAlpha,BRLayerWeatherAlpha)
    tmpBRLayerWeatherAlpha=BRLayerWeatherAlpha
    tmpLayerWeatherAlpha=LayerWeatherAlpha
    Call DepositionMatrixCalc (NumGrade,MeanGrading(:),DepositionMatrix)
	Call CreepRateCalc(CreepRate,NumLayers,LayerDepth,ArmourDepth,Slope,CreepDepth,SurfCreepProp,SurfCreepRate,NCols,NRows,dt,CellSize,xllcor,yllcor)
    Do ig=1,NumGrade
        SurfaceGrad(:,:,ig)=InitGrad(ig);
    end do
    ProfileGrading=0
    ProfileGrading(:,:,:,NumGrade)=100.0

!@#$$$$$$$$$$ Iteration loop starts here $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$           
 count=1
 ChangeCount=1
 Do j=0, NumIter  !number of loops (Iterations)  
   tmpJunc0=Junc0
   tmpJuncsions=Junc
   UpSlopFlowGrad=0
   QsUS=0
   tmpUSGrad=0
   if ((j*100)/NumIter==count)then
     tmpp=(Real(j*100)/Real(NumIter))
     Write (*,"(1F6.2,a)") tmpp,'%'
     count=count+1
   End If
    indx=1
    c=tmpJunc0(indx,1)
    r=tmpJunc0(indx,2)
! Climate change implemantaion !
!   If (ClimatFluc.NE.0) Then
! Randomly fluctuating runoff at a range of +/-50%
       !RandChange=1 !no random
       RandChange=(ran2(seed1)+0.5)
       tmpQf = Qf*RandChange
    !   Do i=1,ClimateRec 
	If (j.EQ.(NumIter*IterPrecent(ChangeCount)/100)+1) then
       StepWeatheringChange = (WeatheringChange(ChangeCount)-WeatheringChange(ChangeCount+1)) / ((NumIter*(IterPrecent(ChangeCount)/100))-(NumIter*(IterPrecent(ChangeCount+1)/100)))
       StepErosionChange = (ErosionChange(ChangeCount)-ErosionChange(ChangeCount+1)) / ((NumIter*(IterPrecent(ChangeCount)/100))-(NumIter*(IterPrecent(ChangeCount+1)/100)))
       StepDischargeChange = (DischargeChange(ChangeCount)-DischargeChange(ChangeCount+1)) / ((NumIter*(IterPrecent(ChangeCount)/100))-(NumIter*(IterPrecent(ChangeCount+1)/100)))
       StepAeolianChange = (AeolianChange(ChangeCount)-AeolianChange(ChangeCount+1)) / ((NumIter*(IterPrecent(ChangeCount)/100))-(NumIter*(IterPrecent(ChangeCount+1)/100)))
	   StepCreepChange = (CreepChange(ChangeCount)-CreepChange(ChangeCount+1)) / ((NumIter*(IterPrecent(ChangeCount)/100))-(NumIter*(IterPrecent(ChangeCount+1)/100)))
       ChangeCount=ChangeCount+1
	End If
  !     End Do 
       tmpWeatheringChange = tmpWeatheringChange + StepWeatheringChange 
       tmpErosionChange = tmpErosionChange + StepErosionChange  
       tmpDischargeChange = tmpDischargeChange + StepDischargeChange
	   tmpAeolianChange = tmpAeolianChange + StepAeolianChange
	   tmpCreepChange = tmpCreepChange + StepCreepChange
       tmpkErode = kErode * tmpErosionChange
       WeatherChange = RandChange*tmpWeatheringChange
	   tmpAeolianRate=AeolianRate*tmpAeolianChange
	 ! tmpCreepRate = 1*tmpCreepChange
	   tmpCreepRate = 1*tmpCreepChange*RandChange
 
!!!!!^^^^^^^^^^^^^^^^^^^^^^^^^ Do While loop^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    
    Do while (tmpJunc0(indx,1)>0)
        If (tmpJunc0(indx,1)>0) then
            dir=FlowDir(c,r)
            to_c=c+ND(dir,1)
            to_r=r+ND(dir,2)
            If (dir.lt.0) then
                 d50map(c,r)=-9999
                 cycle
            End If
            TempSurfGrad=SurfaceGrad(c,r,:)
            TempLayerGrad=ProfileGrading(c,r,:,:)
            Call Discharge ((tmpQf*tmpDischargeChange), FlowAcc(c,r), CellSize, Slope(c,r) ,GradingSize(:),NumGrade, Qd(c,r), SizeTh(c,r), Dmin(c,r) &
				,TempLayerGrad,Numlayers,LayerDepth,ArmourDepth)
            cd50 = d50(TempSurfGrad,GradingSize ,NumGrade)
           
		    If (Mode.EQ.10.OR.Mode.EQ.11) Then !The soil reverse exponential modes
                If (TempLayerGrad(1,NumGrade).EQ.100) then
                    WeatherFront=2
                Else     
                Do l=2, Numlayers
                    If (TempLayerGrad(l,NumGrade).EQ.100.AND.TempLayerGrad(l-1,NumGrade).LT.100) then
                        WeatherFront=l
                        Exit
                    End IF
                End Do
                End If
                i=Numlayers+1
                tmpLayerWeatherAlpha=0
                Do l=WeatherFront,1,-1
                    tmpLayerWeatherAlpha(l)= LayerWeatherAlpha(i) 
                    i=i-1
                End do    
            End If
!!##################################################################################################################  
!! Surface Weathering 
            if (WeatherAlpha > 0) then
                TempGrad=0
                Do i=1, NumGrade-1
                    Do ii=i, NumGrade-1
                          TempGrad= TempGrad + (TempSurfGrad(ii)* LayerWeatheringM(1,ii,i)*WeatherChange*tmpLayerWeatherAlpha(1))      
                    End Do
                    TempSurfGrad(i)=(TempSurfGrad(i)*(1-(tmpLayerWeatherAlpha(1)* WeatherChange))) + TempGrad
                    TempGrad=0
                End Do 
                Do i=1, NumGrade
                    Do ii=NumGrade, NumGrade
                        TempGrad= TempGrad + (TempSurfGrad(ii)* LayerWeatheringM(1,ii,i)*WeatherChange*tmpBRLayerWeatherAlpha(1))
                    End Do
                    if (i.LT.NumGrade) then 
                        TempSurfGrad(i)=TempSurfGrad(i) + TempGrad
                    else
                        TempSurfGrad(NumGrade)=(TempSurfGrad(NumGrade)*(1-(tmpBRLayerWeatherAlpha(1)* WeatherChange))) + TempGrad
                    End If
                    TempGrad=0
                End Do  
!! End of surface Weathering
!!#####################################################################################################################
!!#####################################################################################################################
!! Sub-Surface Weathering 
                Do l=1, Numlayers
                    If (l.EQ.1) Then
                        TempGrad=0
                        Do i=1, NumGrade-1
                            Do ii=i, NumGrade-1
                               TempGrad= TempGrad + (TempLayerGrad(l,ii)* LayerWeatheringM(l+1,ii,i)*WeatherChange*tmpLayerWeatherAlpha(l+1))    
                            End Do
                            TempLayerGrad(l,i)=(TempLayerGrad(l,i)*(1-(tmpLayerWeatherAlpha(l+1)* WeatherChange))) + TempGrad
                            TempGrad=0
                        End Do 
                        Do i=1, NumGrade
                            Do ii=NumGrade, NumGrade
                                TempGrad= TempGrad + (TempLayerGrad(l,ii)* LayerWeatheringM(l+1,ii,i)*WeatherChange*tmpBRLayerWeatherAlpha(l+1))
                            End Do
                            if (i.LT.NumGrade) then 
                                TempLayerGrad(l,i)=TempLayerGrad(l,i) + TempGrad
                            Else
                                TempLayerGrad(l,NumGrade)=(TempLayerGrad(l,NumGrade)*(1-(tmpBRLayerWeatherAlpha(l+1)* WeatherChange))) + TempGrad
                            End If
                            TempGrad=0
                        End Do  
                    Else If (TempLayerGrad(l-1,NumGrade).LT.20) Then !!!!! Threshold for weathering transition down the profile (added v5.3) 

                        TempGrad=0
                        Do i=1, NumGrade-1
                            Do ii=i, NumGrade-1
                                TempGrad= TempGrad + (TempLayerGrad(l,ii)* LayerWeatheringM(l+1,ii,i)*WeatherChange*tmpLayerWeatherAlpha(l+1))  
                            End Do
                            TempLayerGrad(l,i)=(TempLayerGrad(l,i)*(1-(tmpLayerWeatherAlpha(l+1)* WeatherChange))) + TempGrad
                            TempGrad=0
                        End Do
                        Do i=1, NumGrade
                            Do ii=NumGrade, NumGrade
                                TempGrad= TempGrad + (TempLayerGrad(l,ii)* LayerWeatheringM(l+1,ii,i)*WeatherChange*tmpBRLayerWeatherAlpha(l+1))
                            End Do
                            if (i.LT.NumGrade) then 
                                TempLayerGrad(l,i)=TempLayerGrad(l,i) + TempGrad
                            Else
                                TempLayerGrad(l,NumGrade)=(TempLayerGrad(l,NumGrade)*(1-(tmpBRLayerWeatherAlpha(l+1)* WeatherChange))) + TempGrad
                            End If
                            TempGrad=0
                        End Do  
                    End If                    
                End Do
             End If
!! End of Weathering
!!#########################################################################################################################  
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!Start Aeolian Deposition
        If (tmpAeolianRate.GT.0) then
            TempSurfGrad=TempSurfGrad+(AeolDepGrad*tmpAeolianRate)
            DepoProp=(sum(TempSurfGrad)-100)/100
            tmpBRprop=TempLayerGrad(1,NumGrade)
            PushProp=DepoProp*(ArmourDepth/LayerDepth)
            If (tmpBRprop.LT.100) then !Layer 1 has soil
                TempSurfGrad=TempSurfGrad*(100/sum(TempSurfGrad))
                TempLayerGrad(1,:)=TempLayerGrad(1,:)+TempSurfGrad*PushProp
            Else !layer 1 is all bedrock
                If (TempSurfGrad(NumGrade).GT.(DepoProp*100)) then ! more bedrock at the surface then deposition material
                    TempSurfGrad(NumGrade)=TempSurfGrad(NumGrade)-(DepoProp*100)
                Else ! less bedrock then deposeted
                    TempSurfGrad(NumGrade) = 0
                    DepoProp =( sum(TempSurfGrad)-100)/100
                    PushProp = DepoProp*(ArmourDepth/LayerDepth)
                    TempSurfGrad = TempSurfGrad*(100/sum(TempSurfGrad))
                    TempLayerGrad(1,:) = TempLayerGrad(1,:)+TempSurfGrad*PushProp
                    TempLayerGrad(1,:) =  TempLayerGrad(1,:)*(100/sum(TempLayerGrad(1,:)))
                End If
            End If
            Do l=1, NumLayers-1
                DepoProp=(sum(TempLayerGrad(l,:))-100)/100
                tmpBRprop=TempLayerGrad(l+1,NumGrade)
                If (tmpBRprop.LT.100) then !Layer l+1 has soil
                    TempLayerGrad(l,:) = TempLayerGrad(l,:)*(100/sum(TempLayerGrad(l,:))) 
                    TempLayerGrad(l+1,:) = TempLayerGrad(l+1,:) + TempLayerGrad(l,:)*DepoProp 
                Else !Layer l+1 is all bedrock
                    If (TempLayerGrad(l,NumGrade).GT.(DepoProp*100)) then ! more bedrock at the layer l then deposition material
                        TempLayerGrad(l,NumGrade) = TempLayerGrad(l,NumGrade)-(DepoProp*100)
                    Else ! less bedrock in layer l then deposeted (need to remove some bedrock and some soil)
                        TempLayerGrad(l,NumGrade) = 0
                        DepoProp =( sum(TempLayerGrad(l,:))-100)/100
                        TempLayerGrad(l,:) =  TempLayerGrad(l,:)*(100/sum(TempLayerGrad(l,:)))
                        TempLayerGrad(l+1,:) = TempLayerGrad(l+1,:) + TempLayerGrad(l,:)*DepoProp
                        TempLayerGrad(l+1,:) =  TempLayerGrad(l+1,:)*(100/sum(TempLayerGrad(l+1,:)))
                    End If                
                    Exit
                End If
            End Do 
            If (TempLayerGrad(NumLayers,NumGrade).LT.100) then
                TempLayerGrad(NumLayers,:) =  TempLayerGrad(NumLayers,:)*(100/sum(TempLayerGrad(NumLayers,:)))
            End If
          End If !aeolian Deposition
!End Aeolian Deposition            
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@         
!EROSION Calculation 
			a=1/(sum(((TempSurfGrad/100)**2)/((MeanGrading/1000)**1)))
            Qs = ((((tmpkErode*(Qd(c,r))**1)* Slope(c,r)**1.2)/((((1.65**2)*(cd50/1000)**0.5)))))*dt ![m3/iter]
            E = Qs/(ArmourDepth/100)!Erosion rate [proportion of surface layer removed]
            Erosion = E-QsUS(c,r) !delta Z
            ErodedGrad=0
            DepositedGrad=0
            PotenDepo=0
			If (Dmin(c,r).EQ.NumGrade) Dmin(c,r)=NumGrade-1
            If (Erosion.GT.0) Then !EROSION case
                ErodedGrad(1:Dmin(c,r))=((TempSurfGrad(1:Dmin(c,r))/100))*(Erosion)*((a*(TempSurfGrad(1:Dmin(c,r))/100))/((MeanGrading(1:Dmin(c,r))/1000)**1))
                ErosionGrad(1:Dmin(c,r))= 100*(((TempSurfGrad(1:Dmin(c,r))/100)) - ErodedGrad(1:Dmin(c,r)))
               
                !calculating last erodible class
                if (Dmin(c,r)<NumGrade) then
                    ErosionGrad(Dmin(c,r))=(100*(((TempSurfGrad(Dmin(c,r))/100))-((TempSurfGrad(Dmin(c,r))/100))* (Erosion)*(((SizeTh(c,r)-GradingSize(Dmin(c,r)))&
                        /(GradingSize(Dmin(c,r)+1)-GradingSize(Dmin(c,r)))))*((a*(TempSurfGrad(Dmin(c,r))/100))/((MeanGrading(Dmin(c,r))/1000)**1))))
                end If
                ! Larger then erodible
                ErosionGrad(Dmin(c,r)+1:NumGrade)=TempSurfGrad(Dmin(c,r)+1:NumGrade)
                Do i=1,NumGrade
                   if (ErosionGrad(i)<0) ErosionGrad(i)=0 
                End Do
                UpSlopFlowGrad(to_c,to_r,:)=UpSlopFlowGrad(to_c,to_r,:) + UpSlopFlowGrad(c,r,:)+(TempSurfGrad-ErosionGrad)
                NetErosion=100-sum(ErosionGrad) !Positive for erosion negative for deposition
                QsUS(to_c,to_r)=QsUS(to_c,to_r)+ (NetErosion/100) !Qs
            Else 
                If (Erosion.LT.0) then !DEPOSITION case
                    tmpUSGrad=UpSlopFlowGrad(c,r,:)
                    If (sum(tmpUSGrad).GT.0) then
                        b=1/(sum((tmpUSGrad/100)* DepositionMatrix))
                        PotenDepo = -(Erosion * b*(tmpUSGrad/100)* DepositionMatrix)*100
                        Do i=1,NumGrade
                            If (PotenDepo(i).LT.tmpUSGrad(i)) then
                                DepositedGrad(i)=PotenDepo(i)
                            else
                                DepositedGrad(i)=tmpUSGrad(i)
                            End If
                        End Do
                    else
                    DepositedGrad=0
                    End If
                    ErosionGrad = TempSurfGrad + DepositedGrad
                    Do i=1,NumGrade
                        if (ErosionGrad(i)<0) ErosionGrad(i)=0 
                    End Do
                    UpSlopFlowGrad(to_c,to_r,:)=UpSlopFlowGrad(to_c,to_r,:) + (tmpUSGrad-DepositedGrad) 
					NetErosion=100-sum(ErosionGrad) !Positive for erosion negative for deposition
                    QsUS(to_c,to_r)=QsUS(to_c,to_r)+ (sum(tmpUSGrad-DepositedGrad))/100
                End If
            End If 
            
            If  (Erosion.EQ.0) then
                NetErosion=0
                QsUS(to_c,to_r)=QsUS(to_c,to_r)+ 0
            End If        
!!End of EROSION Calculation
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@      
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$           
!Resupply Calculation            
        If (NetErosion.NE.0.AND.Erosion.NE.0) then
            If (Erosion.GT.0) Then ! Erosion case^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    ! Resupplying sediments from underlying layer according to the amount of sediment that have been removed by erosion
                ResupplyGrad=0
                tmpBRprop=TempLayerGrad(1,NumGrade)/100 !proportion of bedrock in layer 1
                PushProp=(NetErosion/100)*(ArmourDepth/LayerDepth) !proportion of material transition from layer1 to surface
                If ((1-tmpBRprop).LT.PushProp) then !If there is less soil in layer1 then have been eroded
                    ResupplyGrad(1:NumGrade-1)=TempLayerGrad(1,1:NumGrade-1)
                    ResupplyGrad(NumGrade)=NetErosion-sum(TempLayerGrad(1,1:NumGrade-1))
                    TempLayerGrad(1,NumGrade)=100
                    TempLayerGrad(1,1:NumGrade-1)=0 !all soil particels been transfered to surface
                Else !Layer1 has more soil then eroded
                    ResupplyGrad(1:NumGrade-1)= TempLayerGrad(1,1:NumGrade-1)*(100/sum(TempLayerGrad(1,1:NumGrade-1)))*(NetErosion/100)
                    TempLayerGrad(1,1:NumGrade-1)=TempLayerGrad(1,1:NumGrade-1)-(TempLayerGrad(1,1:NumGrade-1)*(100/sum(TempLayerGrad(1,1:NumGrade-1)))*PushProp)
                    Do l=1, Numlayers-1
                        tmpBRprop=TempLayerGrad(l+1,NumGrade)/100
                        If (tmpBRprop.LT.1) then 
                            If ((1-tmpBRprop).LT.PushProp) then !Layer l+1 has less soil then erosion proportion
                                TempLayerGrad(l,1:NumGrade-1) = TempLayerGrad(l,1:NumGrade-1) + TempLayerGrad(l+1,1:NumGrade-1) !adding all the soil to from l+1 to l
                                TempLayerGrad(l,NumGrade) = TempLayerGrad(l,NumGrade) + ((PushProp*100)-(sum(TempLayerGrad(l+1,1:NumGrade-1)))) !adding bedrock to layer l
                                TempLayerGrad(l+1,1:NumGrade-1) = 0
                                TempLayerGrad(l+1,NumGrade) = 100
                            Else ! Layer l+1 has more soil the erosion proportion
                                TempLayerGrad(l,1:NumGrade-1) = TempLayerGrad(l,1:NumGrade-1) + TempLayerGrad(l+1,1:NumGrade-1)* PushProp
                                TempLayerGrad(l+1,1:NumGrade-1) = TempLayerGrad(l+1,1:NumGrade-1) - TempLayerGrad(l+1,1:NumGrade-1)* PushProp
                            End If    
                        Else
                            TempLayerGrad(l,NumGrade)=TempLayerGrad(l,NumGrade) + 100*PushProp
							TempLayerGrad(l,:) =  TempLayerGrad(l,:)*(100/sum(TempLayerGrad(l,:)))
                            Exit
                        End If
                    End Do
                    If (TempLayerGrad(Numlayers,NumGrade).LT.100) then
                        TempLayerGrad(Numlayers,NumGrade)=  TempLayerGrad(Numlayers,NumGrade)+ 100*PushProp
                        TempLayerGrad(Numlayers,:)=TempLayerGrad(Numlayers,:)*(100/sum(TempLayerGrad(Numlayers,:)))
                    End If
                End If
     ! Adding the resupplyed material to the eroded armour layer
                TempSurfGrad=ErosionGrad+ResupplyGrad

            Else !Deposition case ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                TempSurfGrad=ErosionGrad !*(100/sum(ErosionGrad))
                DepoProp=(sum(TempSurfGrad)-100)/100
                tmpBRprop=TempLayerGrad(1,NumGrade)
                PushProp=DepoProp*(ArmourDepth/LayerDepth)
                If (tmpBRprop.LT.100) then !Layer 1 has soil
                    TempSurfGrad=TempSurfGrad*(100/sum(TempSurfGrad))
                    TempLayerGrad(1,:)=TempLayerGrad(1,:)+TempSurfGrad*PushProp
                Else !layer 1 is all bedrock
                    If (TempSurfGrad(NumGrade).GT.(DepoProp*100)) then ! more bedrock at the surface then deposition material
                        TempSurfGrad(NumGrade)=TempSurfGrad(NumGrade)-(DepoProp*100)
                    Else ! less bedrock then deposeted
                        TempSurfGrad(NumGrade) = 0
                        DepoProp =( sum(TempSurfGrad)-100)/100
                        PushProp = DepoProp*(ArmourDepth/LayerDepth)
                        TempSurfGrad = TempSurfGrad*(100/sum(TempSurfGrad))
                        TempLayerGrad(1,:) = TempLayerGrad(1,:)+TempSurfGrad*PushProp
                        TempLayerGrad(1,:) =  TempLayerGrad(1,:)*(100/sum(TempLayerGrad(1,:)))
                    End If
                End If
                Do l=1, NumLayers-1
                    DepoProp=(sum(TempLayerGrad(l,:))-100)/100
                    tmpBRprop=TempLayerGrad(l+1,NumGrade)
                    If (tmpBRprop.LT.100) then !Layer l+1 has soil
                        TempLayerGrad(l,:) = TempLayerGrad(l,:)*(100/sum(TempLayerGrad(l,:))) 
                        TempLayerGrad(l+1,:) = TempLayerGrad(l+1,:) + TempLayerGrad(l,:)*DepoProp 
                    Else !Layer l+1 is all bedrock
                        If (TempLayerGrad(l,NumGrade).GT.(DepoProp*100)) then ! more bedrock at the layer l then deposition material
                            TempLayerGrad(l,NumGrade) = TempLayerGrad(l,NumGrade)-(DepoProp*100)
                        Else ! less bedrock in layer l then deposeted (need to remove some bedrock and some soil)
                            TempLayerGrad(l,NumGrade) = 0
                            DepoProp =( sum(TempLayerGrad(l,:))-100)/100
                            TempLayerGrad(l,:) =  TempLayerGrad(l,:)*(100/sum(TempLayerGrad(l,:)))
                            TempLayerGrad(l+1,:) = TempLayerGrad(l+1,:) + TempLayerGrad(l,:)*DepoProp
                            TempLayerGrad(l+1,:) =  TempLayerGrad(l+1,:)*(100/sum(TempLayerGrad(l+1,:)))
                        End If                
                        Exit
                    End If
                End Do  
			End If
        End If
!End Resupply
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Start Soil Creep
        ! Displacement by Creep
        If (SurfCreepProp(c,r).GT.0) then
            TempSurfGrad = TempSurfGrad - (TempSurfGrad*SurfCreepProp(c,r)*tmpCreepRate*CreepDepth(1)) !Displacing a portion of the surface layer
            SurfaceGrad(to_c,to_r,:) = SurfaceGrad(to_c,to_r,:) + TempSurfGrad*SurfCreepProp(c,r)*tmpCreepRate*CreepDepth(1) !Adding this portion to downstream node
            tmpBRprop=TempLayerGrad(1,NumGrade)
            If (tmpBRprop.LT.100) then
               TempLayerGrad(1,1:NumGrade-1)= TempLayerGrad(1,1:NumGrade-1)- TempLayerGrad(1,1:NumGrade-1)*SurfCreepProp(c,r)*tmpCreepRate*CreepDepth(2)
               If (ProfileGrading(to_c,to_r,1,NumGrade).LT.100) Then
                   ProfileGrading(to_c,to_r,1,1:NumGrade-1) = ProfileGrading(to_c,to_r,1,1:NumGrade-1) + TempLayerGrad(1,1:NumGrade-1)*SurfCreepProp(c,r)*tmpCreepRate*CreepDepth(2)
               Else
                   SurfaceGrad(to_c,to_r,:) = SurfaceGrad(to_c,to_r,:) + TempLayerGrad(1,1:NumGrade-1)*SurfCreepProp(c,r)*tmpCreepRate*CreepDepth(2)
               End If
            End IF
            Do l=2, NumLayers !Displacing material from layers with soil and adding it to the layer downstream
               tmpBRprop=TempLayerGrad(l,NumGrade)
               If (tmpBRprop.LT.100) then
                   TempLayerGrad(l,1:NumGrade-1)= TempLayerGrad(l,1:NumGrade-1)- TempLayerGrad(l,1:NumGrade-1)*SurfCreepProp(c,r)*tmpCreepRate*CreepDepth(l+1)
                   Do i=l, 1, -1 !Add material to downstram layer which is not badrock
                       If (ProfileGrading(to_c,to_r,i,NumGrade).LT.100) Then
                           ProfileGrading(to_c,to_r,i,1:NumGrade-1) = ProfileGrading(to_c,to_r,i,1:NumGrade-1) + TempLayerGrad(l,1:NumGrade-1)*SurfCreepProp(c,r)*tmpCreepRate*CreepDepth(l+1)
                           Exit
                       End If
                  End Do 
               Else 
                   Exit
               End If
             End Do 
             
 ! Resupplying local layers       
            tmpDeficit = 100-sum(TempSurfGrad) !how much deficit at the surface layer due to creep
            ResupplyGrad=0
            tmpBRprop=TempLayerGrad(1,NumGrade)/100 !proportion of bedrock in layer 1
            PushProp=(tmpDeficit/100)*(ArmourDepth/LayerDepth) !proportion of material transition from layer1 to surface
            If ((1-tmpBRprop).LT.PushProp) then !If there is less soil in layer1 then have been eroded
                ResupplyGrad(1:NumGrade-1)=TempLayerGrad(1,1:NumGrade-1)
                ResupplyGrad(NumGrade)=tmpDeficit-sum(TempLayerGrad(1,1:NumGrade-1))
                TempLayerGrad(1,NumGrade)=TempLayerGrad(1,NumGrade)-sum(TempLayerGrad(1,1:NumGrade-1))
				TempLayerGrad(1,NumGrade) = 100
                TempLayerGrad(1,1:NumGrade-1)=0 !all soil particels been transfered to surface
            Else !Layer1 has more soil then eroded
                ResupplyGrad(1:NumGrade-1)= TempLayerGrad(1,1:NumGrade-1)*(100/sum(TempLayerGrad(1,1:NumGrade-1)))*(tmpDeficit/100)
                TempLayerGrad(1,1:NumGrade-1)=TempLayerGrad(1,1:NumGrade-1)-(TempLayerGrad(1,1:NumGrade-1)*(100/sum(TempLayerGrad(1,1:NumGrade-1)))*PushProp)
                Do l=1, Numlayers-1
                    tmpBRprop=TempLayerGrad(l+1,NumGrade)/100
                    tmpDeficit=100-sum(TempLayerGrad(l,:))
                    PushProp=(tmpDeficit/100)
                    If (tmpBRprop.LT.1) then 
                        If ((1-tmpBRprop).LT.PushProp) then !Layer l+1 has less soil then erosion proportion
                            TempLayerGrad(l,1:NumGrade-1) = TempLayerGrad(l,1:NumGrade-1) + TempLayerGrad(l+1,1:NumGrade-1) !adding all the soil to from l+1 to l
                            TempLayerGrad(l,NumGrade) = TempLayerGrad(l,NumGrade) + ((PushProp*100)-(sum(TempLayerGrad(l+1,1:NumGrade-1)))) !adding bedrock to layer l
                            TempLayerGrad(l+1,1:NumGrade-1) = 0
                            TempLayerGrad(l+1,NumGrade) = 100
                        Else ! Layer l+1 has more soil the erosion proportion
                            TempLayerGrad(l,1:NumGrade-1) = TempLayerGrad(l,1:NumGrade-1) + (TempLayerGrad(l+1,1:NumGrade-1)*(100/sum(TempLayerGrad(l+1,1:NumGrade-1)))* PushProp)
                            TempLayerGrad(l+1,1:NumGrade-1) = TempLayerGrad(l+1,1:NumGrade-1) - (TempLayerGrad(l+1,1:NumGrade-1)*(100/sum(TempLayerGrad(l+1,1:NumGrade-1)))* PushProp)
                        End If    
                    Else
                        TempLayerGrad(l,NumGrade)=TempLayerGrad(l,NumGrade) + 100*PushProp
                        Exit
                    End If
                End Do
                If (TempLayerGrad(Numlayers,NumGrade).LT.100) then
                    TempLayerGrad(Numlayers,NumGrade)=  TempLayerGrad(Numlayers,NumGrade)+ 100*PushProp
                    TempLayerGrad(Numlayers,:)=TempLayerGrad(Numlayers,:)*(100/sum(TempLayerGrad(Numlayers,:)))
                End If
            End If
            TempSurfGrad=TempSurfGrad+ResupplyGrad
            SurfaceGrad(c,r,:)=TempSurfGrad
            ProfileGrading(c,r,:,:)=TempLayerGrad

    ! Resupplying Downstream layers 
            TempSurfGrad2 = SurfaceGrad(to_c,to_r,:)
            TempLayerGrad2 = ProfileGrading(to_c,to_r,:,:) 
            
            DepoProp=(sum(TempSurfGrad2)-100)/100
            tmpBRprop=TempLayerGrad2(1,NumGrade)
            PushProp=DepoProp*(ArmourDepth/LayerDepth)
            If (tmpBRprop.LT.100) then !Layer 1 has soil
                TempSurfGrad2=TempSurfGrad2*(100/sum(TempSurfGrad2))
                TempLayerGrad2(1,:)=TempLayerGrad2(1,:)+TempSurfGrad2*PushProp
            Else !layer 1 is all bedrock
                If (TempSurfGrad2(NumGrade).GT.(DepoProp*100)) then ! more bedrock at the surface then deposition material
                    TempSurfGrad2(NumGrade)=TempSurfGrad2(NumGrade)-(DepoProp*100)
                Else ! less bedrock then deposeted
                    TempSurfGrad2(NumGrade) = 0
                    DepoProp =( sum(TempSurfGrad2)-100)/100
                    PushProp = DepoProp*(ArmourDepth/LayerDepth)
                    TempSurfGrad2 = TempSurfGrad2*(100/sum(TempSurfGrad2))
                    TempLayerGrad2(1,:) = TempLayerGrad2(1,:)+TempSurfGrad2*PushProp
                    TempLayerGrad2(1,:) =  TempLayerGrad2(1,:)*(100/sum(TempLayerGrad2(1,:)))
                End If
            End If
            Do l=1, NumLayers-1
                DepoProp=(sum(TempLayerGrad2(l,:))-100)/100
                tmpBRprop=TempLayerGrad2(l+1,NumGrade)
                If (tmpBRprop.LT.100) then !Layer l+1 has soil
                    TempLayerGrad2(l,:) = TempLayerGrad2(l,:)*(100/sum(TempLayerGrad2(l,:))) 
                    TempLayerGrad2(l+1,:) = TempLayerGrad2(l+1,:) + TempLayerGrad2(l,:)*DepoProp 
                Else !Layer l+1 is all bedrock
                    If (TempLayerGrad2(l,NumGrade).GT.(DepoProp*100)) then ! more bedrock at the layer l then deposition material
                        TempLayerGrad2(l,NumGrade) = TempLayerGrad2(l,NumGrade)-(DepoProp*100)
                    Else ! less bedrock in layer l then deposeted (need to remove some bedrock and some soil)
                        TempLayerGrad2(l,NumGrade) = 0
                        DepoProp =( sum(TempLayerGrad2(l,:))-100)/100
                        TempLayerGrad2(l,:) =  TempLayerGrad2(l,:)*(100/sum(TempLayerGrad2(l,:)))
                        TempLayerGrad2(l+1,:) = TempLayerGrad2(l+1,:) + TempLayerGrad2(l,:)*DepoProp
                        TempLayerGrad2(l+1,:) =  TempLayerGrad2(l+1,:)*(100/sum(TempLayerGrad2(l+1,:)))
                    End If                
                    Exit
                End If
            End Do  
            SurfaceGrad(to_c,to_r,:)=TempSurfGrad2
            ProfileGrading(to_c,to_r,:,:)=TempLayerGrad2     
        End If !soil creep    
! End Soil Creep
!PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
        SurfaceGrad(c,r,:)=TempSurfGrad
        
     !   if (sum(TempLayerGrad) .GT. 2000.0) then
      !      write(*,*) sum(TempLayerGrad)
      !  end if   
        ProfileGrading(c,r,:,:)=TempLayerGrad
        d50map(c,r)=cd50
        FinalNetErosion(c,r)=(NetErosion)
        TotalErosion(c,r)=TotalErosion(c,r)+ (NetErosion)
         
        If (j==NumIter) then
            FinalErosion(c,r)=(Erosion/dt)
        End If

            If (tmpJuncsions(to_c,to_r)==1.0) then
                c=to_c
                r=to_r
            Else
                tmpJunc0(indx,:)=0
                indx=indx+1 
                c=tmpJunc0(indx,1)
                r=tmpJunc0(indx,2)
            End If
            tmpJuncsions(to_c,to_r)=tmpJuncsions(to_c,to_r)-1
        !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        End if !>0
    End Do !while
    
    !Write intermediate d50 layers
    If (OUTcount.LE.NumOutputs) then
        If (j.EQ.OutputTS(OUTcount))then
            Write ((590+OUTcount),1001) d50map     
            Write ((790+OUTcount),1009) FinalNetErosion
            ! Calculate intermedeate SOIL DEPTH
            Do cc=1, NCols
                Do rr=1, NRows
                    If (d50map(cc,rr)==-9999) then
                        LayerD50(cc,rr)=-9999
                        cycle
                    End if
                    Depth(cc,rr)=(NumLayers*LayerDepth)-(sum((ProfileGrading(cc,rr,:,NumGrade)/100)*LayerDepth))
                End Do
            End Do        
            Write ((690+OUTcount),1001) depth
            OUTcount=OUTcount+1
        End If
    End If
 End Do ! End of Iteration loop^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! Calculating soil depth for the output file
    Do c=1, NCols
        Do r=1, NRows
            If (d50map(c,r)==-9999) then
                LayerD50(c,r)=-9999
                SurfaceGrad(c,r,:)=-9999
                cycle
            End if
            Depth(c,r)=(NumLayers*LayerDepth)-(sum((ProfileGrading(c,r,:,NumGrade)/100)*LayerDepth))
        End Do
    End Do
        
1009    Format (<NCols>F15.8)
1001    Format (<NCols>F10.3)
1771    Format (<NCols>E12.4)        
    Write (100,1001)d50map   !Final d50 map
    Write (3210,1001)Junc
    Write (500,1009) FinalErosion  
    Write (1600,1009) FinalNetErosion 
    Write (200,1771) TotalErosion 
    Write (300,1001) Depth    
    Do ii=1, NumGrade
        Write (1500+ii,1001) SurfaceGrad(:,:,ii)
    end Do 

!! Calculate and write D50 for all layers
    Do l=1, NumLayers
        Do c=1, NCols
            Do r=1, NRows
                If (d50map(c,r)==-9999) then
                    LayerD50(c,r)=-9999
                    cycle
                End if

                LayerD50(c,r) = d50(ProfileGrading(c,r,l,:),GradingSize ,NumGrade)
            End Do
        End Do
        unit1=100+l
        Write (unit1,1001) LayerD50   
    End Do
    
    CALL GETDAT(tmpyear, tmpmonth, tmpday)
    CALL GETTIM(tmphour, tmpminute, tmpsecond, tmphund)
    Write (1,22) 'End: ',tmpday,'/',tmpmonth,'/',tmpyear,' ',tmphour,':',tmpminute,':',tmpsecond  
    Close (unit=1, status='keep')
    Close (unit=99, status='keep')    
    Close (unit=100, status='keep')
    Close (unit=200, status='keep')
    Close (unit=300, status='keep')
    Close (unit=400, status='keep')
    Close (unit=500, status='keep')
    Close (unit=1600, status='keep')
    
    Do l=1, NumLayers
        unit1=100+l
        Close (unit1, status='keep') 
    End Do
    
    Close (unit=10, status='keep') 
    i=0
    Write(*,*) 'Do you want to convert all maps to ArcGIS layers? (1=Yes,0=No)' 
    Read (*,*) i
    If (i==1) then
        Write (*,*) 'Converting ASCII outputs to ArcGIS rasters'
        i = SYSTEMQQ ('ASCIItoRasterMrARM.py') !Execute the Python script file which converts the output files to Rasters   
    End If
End
 
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!Routine to calculate the Flow Junction from TauDEM flow direction   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   Subroutine Junctions_TauDEM(FlowDir,NCols,NRows,Junc) 
   Integer cj,rj,dir
   Real(8) Junc(NCols,NRows),FlowDir(NCols,NRows)
   Junc=0
   Do cj=1,NCols
        Do rj=1,NRows
          dir=FlowDir(cj,rj)
          if (FlowDir(cj,rj)==-9999) then
                Junc(cj,rj)=-9999
                cycle
          End if
          if (cj==1.AND.dir==4.OR.cj==1.AND.dir==5.OR.cj==1.AND.dir==6)dir=-9999
          if (rj==1.AND.dir==4.OR.rj==1.AND.dir==3.OR.rj==1.AND.dir==2)dir=-9999
          if (rj==NRows.AND.dir==6.OR.rj==NRows.AND.dir==7.OR.rj==NRows.AND.dir==8)dir=-9999
          if (cj==NCols.AND.dir==2.OR.cj==NCols.AND.dir==1.OR.cj==NCols.AND.dir==8)dir=-9999
          
          SELECT CASE (dir)
            case(-9999)
            case (1)
                Junc(cj+1,rj)=Junc(cj+1,rj)+1
            case (2)
                Junc(cj+1,rj-1)=Junc(cj+1,rj-1)+1
            case (3)
                Junc(cj,rj-1)=Junc(cj,rj-1)+1
            case (4)
                Junc(cj-1,rj-1)=Junc(cj-1,rj-1)+1
            case (5)
                Junc(cj-1,rj)= Junc(cj-1,rj)+1                
            case (6)
                Junc(cj-1,rj+1)=Junc(cj-1,rj+1)+1
            case (7)
                Junc(cj,rj+1)=Junc(cj,rj+1)+1            
            case (8)
                Junc(cj+1,rj+1)=Junc(cj+1,rj+1)+1
           End Select
        End Do
   End Do
   End Subroutine
    
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!Routine to calculate the Flow Junction from ArcGIS flow direction   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   Subroutine Junctions_ArcGIS(FlowDir,NCols,NRows,Junc) 
   Integer cj,rj,dir
   Real(8) Junc(NCols,NRows),FlowDir(NCols,NRows)
   Junc=0
   Do cj=1,NCols
        Do rj=1,NRows
          dir=FlowDir(cj,rj)
          if (cj==1.AND.dir==4.OR.cj==1.AND.dir==5.OR.cj==1.AND.dir==6)dir=-9999
          if (rj==1.AND.dir==4.OR.rj==1.AND.dir==3.OR.rj==1.AND.dir==2)dir=-9999
          if (rj==NRows.AND.dir==6.OR.rj==NRows.AND.dir==7.OR.rj==NRows.AND.dir==8)dir=-9999
          if (cj==NCols.AND.dir==2.OR.cj==NCols.AND.dir==1.OR.cj==NCols.AND.dir==8)dir=-9999
          
          SELECT CASE (dir)
            case(-9999) 

            case (1)
                Junc(cj+1,rj)=Junc(cj+1,rj)+1
            case (128)
                Junc(cj+1,rj-1)=Junc(cj+1,rj-1)+1
            case (64)
                Junc(cj,rj-1)=Junc(cj,rj-1)+1
            case (32)
                Junc(cj-1,rj-1)=Junc(cj-1,rj-1)+1
            case (16)
                Junc(cj-1,rj)= Junc(cj-1,rj)+1                
            case (8)
                Junc(cj-1,rj+1)=Junc(cj-1,rj+1)+1
            case (4)
                Junc(cj,rj+1)=Junc(cj,rj+1)+1            
            case (2)
                Junc(cj+1,rj+1)=Junc(cj+1,rj+1)+1
           End Select
        End Do
   End Do
   End Subroutine
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!Routine TO find the d50 of a grading profile assumes the the profile is such that the reading
!is that amount not passing that grading size. Smallest grading must be zero. 
!For best results the largest grading should pass 100% of the material. 
!From the ARMOUR code!!!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    REAL*8 function d50(grade,dgrade,NumGrade)
    IMPLICIT NONE
        INTEGER NumGrade
        REAL*8 grade(NumGrade),dgrade(NumGrade)
        INTEGER i,low,high
        REAL*8 sum,sum1

        sum=0
        sum1=0
        DO 1000 i=1,NumGrade
        sum=sum+grade(i)
        1000 CONTINUE
        IF (sum.eq.0.0) THEN
        d50=1.0*(dgrade(2)+dgrade(3))/2.0
                
        RETURN
        END IF
        sum=sum*0.5
        DO 1010 i=1,NumGrade
        sum1=sum1+grade(i)
        IF (sum1.ge.sum) THEN
          high=i+1
          sum1=sum1-sum
          GO TO 1020
        END IF
    1010 CONTINUE
        !CALL runerror(1,'d50')
        WRITE (*,*) sum,grade
        STOP
    1020 low=high-1
        IF (high.gt.NumGrade) THEN
        d50=dgrade(NumGrade)
        ELSE
        d50=dgrade(high)-sum1*(dgrade(high)-dgrade(low))/grade(low)
        END IF
        RETURN
    END function d50
      
!  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!  Routine TO calculate the Weathering transition matrix   
!  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    Subroutine WeatheringTransitionMatrix (NumGrade, GradingSize,WeatherAlpha, WeatheringTrans, SplitProp)
	IMPLICIT NONE
        Integer NumGrade, xi, xii, ii,sum
        Real(8) WeatherAlpha, GradingSize(NumGrade),MinBr,MaxBr,P,RevPara1,RevPara2
        Real(8) Volume(NumGrade),MeanGrading(NumGrade),RevBreakDia1(NumGrade),RevBreakDia2(NumGrade),SplitProp
        Real(8) WeatheringTrans(NumGrade,NumGrade)
        !xi-Lines
        !xii-Columes
        WeatherAlpha=WeatherAlpha/1000
        sum=0

        RevPara1=(1/SplitProp)**0.333
        RevPara2=(1/(1-SplitProp))**0.333
        Do ii=1, NumGrade 
            RevBreakDia1(ii)=GradingSize(ii)*RevPara1
            RevBreakDia2(ii)=GradingSize(ii)*RevPara2
        End Do    
        
        WeatheringTrans=0
        
        Do xi=NumGrade-1, 1, -1
            Do xii=xi, 1, -1
                MaxBr= RevBreakDia1(xii+1)
                MinBr= RevBreakDia1(xii)  
                If (MaxBr>GradingSize(xi+1)) MaxBr=GradingSize(xi+1)
                If (MinBr>GradingSize(xi+1)) MinBr=GradingSize(xi+1)
                If (MinBr<GradingSize(xi))   MinBr=GradingSize(xi)
                If (MaxBr<GradingSize(xi))   MaxBr=GradingSize(xi)
                P = (MaxBr-MinBr)/(GradingSize(xi+1)-GradingSize(xi))
                !WeatheringTrans(xi,xii) = P * WeatherAlpha*SplitProp
                WeatheringTrans(xi,xii) = P * SplitProp
            End Do
         End Do
         
         Do xi=NumGrade-1, 1, -1
            Do xii=xi, 1, -1
                MaxBr= RevBreakDia2(xii+1)
                MinBr= RevBreakDia2(xii)  
                If (MaxBr>GradingSize(xi+1)) MaxBr=GradingSize(xi+1)
                If (MinBr>GradingSize(xi+1)) MinBr=GradingSize(xi+1)
                If (MinBr<GradingSize(xi))   MinBr=GradingSize(xi)
                If (MaxBr<GradingSize(xi))   MaxBr=GradingSize(xi)
                P = (MaxBr-MinBr)/(GradingSize(xi+1)-GradingSize(xi))
                !WeatheringTrans(xi,xii) = WeatheringTrans(xi,xii) + P * WeatherAlpha*(1-SplitProp)
                WeatheringTrans(xi,xii) = WeatheringTrans(xi,xii) + P *(1-SplitProp)
            End Do
         End Do
         
         !! Calculation of the biggest grading
         Do xii=1, NumGrade-1
                If (GradingSize(xii+1) > GradingSize(NumGrade)*0.8) sum=sum+1
         End Do
         
         Do xii=1, NumGrade-1
                If (GradingSize(xii+1) > GradingSize(NumGrade)*0.8) then
                    WeatheringTrans(NumGrade,xii) = 0.5/sum !* WeatherAlpha
                End If
         End Do       

         WeatheringTrans(NumGrade,NumGrade) = 0.5 !* WeatherAlpha !!! Assumes that half of the coursest class break into itself
         Return
     End Subroutine WeatheringTransitionMatrix
    
!###############################################################################################      
! A function to calculate the decline in weathering rate as a function of depth (new from v4.2)
! It return the WeatheringAlpha which is different for every layer 
! The exponential function is taken from Minasny & McBratney (2006)
! de/dt=Po{Exp(-k1h)-Exp(-k2h)}+Pa ; 
! Their original values are: Po-potential WR=0.25; k1=4; k2=6; Pa- steady state WR=0.005 
! The function was changed in version 4.5.1 to account to close to zero weathering in lower layers
!################################################################################################ 
    REAL*8 Function DepthWeatheringAlpha(WeatherAlpha,i,LayerDepth)
    IMPLICIT NONE
   
    Integer i
    Real*8 WeatherAlpha, LayerDepth, Ratio, Depth

    If (i==1) Then
    Depth=0.5
    Else
    Depth=(i-1)*LayerDepth-(LayerDepth/2)
    End If
    Depth=Depth/100 !convert to meters
    Ratio=(0.25*((EXP(-4*Depth+0.02))-(EXP(-6*Depth))+0))/0.04 !We devide by 0.042 to normalize it

    DepthWeatheringAlpha = WeatherAlpha*Ratio

    Return
    End Function DepthWeatheringAlpha
 
 !###############################################################################################      
! A function to calculate the decline in weathering rate as a function of depth (new from v4.5.0)
! It return the WeatheringAlpha which is different for every layer 
! WR=Wo*Exp^(-k*D) ; 
! Their values are:Wo=1; k=2.5; D-depth (calibrated to the hump equation) 
!################################################################################################ 
    REAL*8 Function DepthWeatheringAlpha2(WeatherAlpha,i,LayerDepth)
    IMPLICIT NONE
        Integer i
        Real*8 WeatherAlpha, LayerDepth, Ratio, Depth

        If (i==1) Then
        Depth=0.5
        Else
        Depth=(i-1)*LayerDepth-(LayerDepth/2)
        End If
        Depth=Depth/100 !convert to meters
        Ratio = 1*EXP(-1.738*Depth)

        DepthWeatheringAlpha2 = WeatherAlpha*Ratio

        Return
    End Function DepthWeatheringAlpha2 

 !###############################################################################################      
! A function to calculate the a reverse exponential soil production function.
! The weatehring rate increase with depth rather then decrease; 
! Adress an assumption of a maximum chemical weathering rate at the bedrock weathering front. 
! It return the WeatheringAlpha which is different for every layer 
! WR=Wo*Exp^(-k*D) ; 
! Their values are:Wo=1; k=2.5; D-depth (calibrated to the hump equation) 
!################################################################################################ 
    REAL*8 Function DepthWeatheringAlpha3(WeatherAlpha,i,LayerDepth)
    IMPLICIT NONE
        Integer i
        Real*8 WeatherAlpha, LayerDepth, Ratio, Depth

        If (i==1) Then
        Depth=0.5
        Else
        Depth=(i-1)*LayerDepth-(LayerDepth/2)
        End If
        Depth=Depth/100 !convert to meters
        Ratio = 0.033739*EXP(1.738*Depth)

        DepthWeatheringAlpha3 = WeatherAlpha*Ratio

        Return
    End Function DepthWeatheringAlpha3 
        
!########################################################################################
!This subroutine calculate the weatheruing rate; account for climate changes.
!Added in 7/10/08 version 4.7
!Modified in version 4.9.5 to include separete bedrock weathering function 
!########################################################################################        
 Subroutine WeatheringRateMatrixes (Mode,tmpWeatherAlpha,LayerDepth,NumLayers,NumGrade,GradingSize,SplitProp,LayerWeatheringM,LayerWeatherAlpha,BRLayerWeatherAlpha)
 IMPLICIT NONE
    Real(8) tmpWeatherAlpha,WeatheringTrans(NumGrade,NumGrade),SplitProp,WI,LayerDepth
    Real(8) GradingSize(NumGrade), LayerWeatherAlpha(NumLayers+1),BRLayerWeatherAlpha(NumLayers+1),LayerWeatheringM(NumLayers+1,NumGrade,NumGrade)
    Integer NumLayers,NumGrade,Mode,i,g,gg
    Real(8) DepthWeatheringAlpha2,DepthWeatheringAlpha,DepthWeatheringAlpha3
    WI=1
        If (Mode.eq.0) Then ! Pure (WR constant in depth and space) 
            If (tmpWeatherAlpha > 0) Then
                Do i=1, NumLayers+1
                   LayerWeatherAlpha(i)=tmpWeatherAlpha ! Pure (WR constant in depth and space) 
                   Call WeatheringTransitionMatrix (NumGrade, GradingSize, LayerWeatherAlpha(i), WeatheringTrans,SplitProp)
                   LayerWeatheringM(i,:,:)=WeatheringTrans
                End Do
                BRLayerWeatherAlpha= LayerWeatherAlpha
            End If 
        End If
        If (Mode.eq.1) Then ! Exponential WR only
             if (tmpWeatherAlpha > 0) Then
                Do i=1, NumLayers+1
                 LayerWeatherAlpha(i)= DepthWeatheringAlpha2(tmpWeatherAlpha,i,LayerDepth) ! Exponential WR only
                 Call WeatheringTransitionMatrix (NumGrade, GradingSize, LayerWeatherAlpha(i), WeatheringTrans,SplitProp)
                 LayerWeatheringM(i,:,:)=WeatheringTrans
                End Do
                BRLayerWeatherAlpha= LayerWeatherAlpha
            End If 
        End If
        
        If (Mode.eq.2) Then ! Hump WR only
             if (tmpWeatherAlpha > 0) Then
                Do i=1, NumLayers+1
                 LayerWeatherAlpha(i)= DepthWeatheringAlpha(tmpWeatherAlpha,i,LayerDepth) ! Hump WR only
                 Call WeatheringTransitionMatrix (NumGrade, GradingSize, LayerWeatherAlpha(i), WeatheringTrans,SplitProp)
                 LayerWeatheringM(i,:,:)=WeatheringTrans
                End Do
                BRLayerWeatherAlpha= LayerWeatherAlpha
            End If 
        End If
        
        If (Mode.eq.3) Then !WI only
            if (tmpWeatherAlpha > 0) Then
                Do i=1, NumLayers+1
                    LayerWeatherAlpha(i) = tmpWeatherAlpha * WI !WI only
                    Call WeatheringTransitionMatrix (NumGrade, GradingSize, LayerWeatherAlpha(i), WeatheringTrans,SplitProp)
                    LayerWeatheringM(i,:,:)=WeatheringTrans
                End Do
                BRLayerWeatherAlpha= LayerWeatherAlpha
            End If 
        End IF
        
       If (Mode.eq.4) Then !WI + Exponential WR
            if (tmpWeatherAlpha > 0) Then
                Do i=1, NumLayers+1
                    LayerWeatherAlpha(i) = DepthWeatheringAlpha2(tmpWeatherAlpha,i,LayerDepth) * WI !WI and exponential WR
                    Call WeatheringTransitionMatrix (NumGrade, GradingSize, LayerWeatherAlpha(i), WeatheringTrans,SplitProp)
                    LayerWeatheringM(i,:,:)=WeatheringTrans
                End Do
                BRLayerWeatherAlpha= LayerWeatherAlpha
            End If 
        End IF
        
        If (Mode.eq.5) Then !WI + Hump WR
            if (tmpWeatherAlpha > 0) Then
                Do i=1, NumLayers+1
                    LayerWeatherAlpha(i) = DepthWeatheringAlpha(tmpWeatherAlpha,i,LayerDepth) * WI !WI and hump WR
                    Call WeatheringTransitionMatrix (NumGrade, GradingSize, LayerWeatherAlpha(i), WeatheringTrans,SplitProp)
                    LayerWeatheringM(i,:,:)=WeatheringTrans
                End Do
                BRLayerWeatherAlpha= LayerWeatherAlpha
            End If 
        End IF
        
        If (Mode.eq.6) Then ! Pure except bedrock exponential 
            If (tmpWeatherAlpha > 0) Then
                Do i=1, NumLayers+1
                   LayerWeatherAlpha(i)=tmpWeatherAlpha*0.28 ! Pure (Weathering constant in depth) 0.28 to have similar avarage rates with the expo and hump function 
                   Call WeatheringTransitionMatrix (NumGrade, GradingSize, LayerWeatherAlpha(i), WeatheringTrans,SplitProp)
                   LayerWeatheringM(i,:,:)=WeatheringTrans
                End Do
                Do i=1, NumLayers+1 ! Assinge a weathering rate only to the bedrock class
                   BRLayerWeatherAlpha(i)= DepthWeatheringAlpha2(tmpWeatherAlpha,i,LayerDepth) ! Humped WR only
                   Call WeatheringTransitionMatrix (NumGrade, GradingSize, BRLayerWeatherAlpha(i), WeatheringTrans,SplitProp)
                   Do g=1,NumGrade
                      Do gg=1,NumGrade
                          If (g.eq.NumGrade.or.gg.eq.NumGrade) LayerWeatheringM(i,g,gg)=WeatheringTrans(g,gg)
                      End Do
                   End Do
                End Do
            End If 
        End If
        
        If (Mode.eq.7) Then ! Pure except bedrock humped 
            If (tmpWeatherAlpha > 0) Then
                Do i=1, NumLayers+1
                   LayerWeatherAlpha(i)=tmpWeatherAlpha*0.28 ! Pure (Weathering constant in depth) 0.28 to have similar avarage rates with the expo and hump function 
                   Call WeatheringTransitionMatrix (NumGrade, GradingSize, LayerWeatherAlpha(i), WeatheringTrans,SplitProp)
                   LayerWeatheringM(i,:,:)=WeatheringTrans
                End Do
                Do i=1, NumLayers+1 ! Assinge a weathering rate only to the bedrock class
                   BRLayerWeatherAlpha(i)= DepthWeatheringAlpha(tmpWeatherAlpha,i,LayerDepth) ! Exponential WR only
                   Call WeatheringTransitionMatrix (NumGrade, GradingSize, BRLayerWeatherAlpha(i), WeatheringTrans,SplitProp)
                   Do g=1,NumGrade
                      Do gg=1,NumGrade
                          If (g.eq.NumGrade.or.gg.eq.NumGrade) LayerWeatheringM(i,g,gg)=WeatheringTrans(g,gg)
                      End Do
                   End Do
                End Do
            End If 
        End If
       
        If (Mode.eq.8) Then ! Reverse Exponential
             if (tmpWeatherAlpha > 0) Then
                Do i=1, NumLayers+1
                 LayerWeatherAlpha(i)= DepthWeatheringAlpha3(tmpWeatherAlpha,i,LayerDepth) ! Exponential
                 Call WeatheringTransitionMatrix (NumGrade, GradingSize, LayerWeatherAlpha(i), WeatheringTrans,SplitProp)
                 LayerWeatheringM(i,:,:)=WeatheringTrans
                End Do
                BRLayerWeatherAlpha= LayerWeatherAlpha
            End If 
        End If
        
       If (Mode.eq.9) Then ! Pure except bedrock which is Reverse Exponential 
            If (tmpWeatherAlpha > 0) Then
                Do i=1, NumLayers+1
                   LayerWeatherAlpha(i)=tmpWeatherAlpha*0.28 ! Pure (Weathering constant in depth) 0.28 to have similar avarage rates with the expo and hump function 
                   Call WeatheringTransitionMatrix (NumGrade, GradingSize, LayerWeatherAlpha(i), WeatheringTrans,SplitProp)
                   LayerWeatheringM(i,:,:)=WeatheringTrans
                End Do
                Do i=1, NumLayers+1 ! Assinge a weathering rate only to the bedrock class
                   BRLayerWeatherAlpha(i)= DepthWeatheringAlpha3(tmpWeatherAlpha,i,LayerDepth) ! Exponential WR only
                   Call WeatheringTransitionMatrix (NumGrade, GradingSize, BRLayerWeatherAlpha(i), WeatheringTrans,SplitProp)
                   Do g=1,NumGrade
                      Do gg=1,NumGrade
                          If (g.eq.NumGrade.or.gg.eq.NumGrade) LayerWeatheringM(i,g,gg)=WeatheringTrans(g,gg)
                      End Do
                   End Do
                End Do
            End If 
        End If
        
        If (Mode.eq.10) Then ! Bedrock Exponential and soil Reverse Exponential
            If (tmpWeatherAlpha > 0) Then
                Do i=1, NumLayers+1
                   LayerWeatherAlpha(i)=DepthWeatheringAlpha3(tmpWeatherAlpha,i,LayerDepth) ! Reverse exponential for soil
                   Call WeatheringTransitionMatrix (NumGrade, GradingSize, LayerWeatherAlpha(i), WeatheringTrans,SplitProp)
                   LayerWeatheringM(i,:,:)=WeatheringTrans
                End Do
                Do i=1, NumLayers+1 ! Assinge a weathering rate only to the bedrock class
                   BRLayerWeatherAlpha(i)= DepthWeatheringAlpha2(tmpWeatherAlpha,i,LayerDepth) ! Exponential WR for bedrock
                   Call WeatheringTransitionMatrix (NumGrade, GradingSize, BRLayerWeatherAlpha(i), WeatheringTrans,SplitProp)
                   Do g=1,NumGrade
                      Do gg=1,NumGrade
                          If (g.eq.NumGrade.or.gg.eq.NumGrade) LayerWeatheringM(i,g,gg)=WeatheringTrans(g,gg)
                      End Do
                   End Do
                End Do
            End If 
        End If
        
        If (Mode.eq.11) Then ! Bedrock Humped and soil Reverse Exponential
            If (tmpWeatherAlpha > 0) Then
                Do i=1, NumLayers+1
                   LayerWeatherAlpha(i)=DepthWeatheringAlpha3(tmpWeatherAlpha,i,LayerDepth) ! Reverse exponential for soil
                   Call WeatheringTransitionMatrix (NumGrade, GradingSize, LayerWeatherAlpha(i), WeatheringTrans,SplitProp)
                   LayerWeatheringM(i,:,:)=WeatheringTrans
                End Do
                Do i=1, NumLayers+1 ! Assinge a weathering rate only to the bedrock class
                   BRLayerWeatherAlpha(i)= DepthWeatheringAlpha(tmpWeatherAlpha,i,LayerDepth) ! Humped WR for bedrock
                   Call WeatheringTransitionMatrix (NumGrade, GradingSize,BRLayerWeatherAlpha(i), WeatheringTrans,SplitProp)
                   Do g=1,NumGrade
                      Do gg=1,NumGrade
                          If (g.eq.NumGrade.or.gg.eq.NumGrade) LayerWeatheringM(i,g,gg)=WeatheringTrans(g,gg)
                      End Do
                   End Do
                End Do
            End If 
        End If
        
 End Subroutine WeatheringRateMatrixes
!#################################################################################################
!This subroutine calculate the deposition matrix as a function of the particle settling velocity
!Added in 6/4/10 version 5.3
!##################################################################################################        
 Subroutine DepositionMatrixCalc (NumGrade,Dmean,DepositionMatrix)
 IMPLICIT NONE      
 Integer NumGrade,id
 Real(8) DepositionMatrix(NumGrade),VGD,VG,Vs(NumGrade),g,sg,nu ,Dmean(NumGrade)
 Data g,sg,nu / 9.81,2.6,0.000001/
 !Using Einstein's approach  (Ref. Henderson (1966) and HEC-6 (1993) )
     Do id=1, NumGrade
          IF (Dmean(id).gt.0.) THEN
             VGD=36.0*nu**2/(g*Dmean(id)**3.0*(sg-1))
             VG=(2.0/3.0+VGD)**0.5-VGD**0.5
             Vs(id)=VG*(g*Dmean(id)*(sg-1))**0.5
          ELSE
             Vs(id)=0.
          END IF
      End Do
      DepositionMatrix=Vs !/sum(Vs)
 End Subroutine DepositionMatrixCalc     
!########################################################################################
!This subroutine write a Phython script file which convert the text ASCII output files 
!to Raster using ArcGIS
!Added in 14/7/08 version 4.4.2
!########################################################################################     
 Subroutine PythonScriptWriter (NumLayers,ClimateRec,IterPrecent,ClimatFluc,NumGrade)
 USE IFPORT
 IMPLICIT NONE
    
    Integer NumLayers,i,y,temp,ClimateRec,ClimatFluc,NumGrade
    Real(8) IterPrecent(ClimateRec)
    Character*100 D50aFile,TEFile,D50aOut,DFile,Name
    character*100 dirname,RevDirName
!   variable dirname must be long enough to hold entire string
    integer(4) istat
    ISTAT = GETCWD (dirname)
 !   IF (ISTAT == 0) write *, 'Current directory is ',dirname
    
    Do i=1, 50
        If (IACHAR(dirname(i:i))==92) then
        dirname(i:i)='/'
        End If
    End Do

999  Format (a,I0,a,F0.2,a)
996  Format (a,a,a,I0,a)    
888  Format (a,a,I0)  
    ! Creating script file
    OPEN(unit=70,file="ASCIItoRasterMrARM.py",status='unknown')
    
    !Writing header
    Write (70,'(A)')'#This file was created by MrARM program'
    Write (70,'(A)')'#to convert ASCII outputs to ArcGIS Raster'
    Write (70,'(A)')'import sys, string, os, arcgisscripting'
    Write (70,'(A)')'gp = arcgisscripting.create()'
    Write (70,'(A)')'gp.AddToolbox("C:/Programs/ESRI/ArcGIS/ArcToolbox/Toolboxes/Conversion Tools.tbx")'
    Write (70,'(A)')'try:'
    
 777 Format (a,a,a)   
    Do i=0,NumLayers
        write (D50aFile,996) TRIM(dirname),'/','D50aL',i,'.txt'
        Write (70,777)'    InAsciiFile = "',TRIM(D50aFile),'"'
        write (D50aOut,888)  TRIM(dirname),'/D50a_L',i
        
        Do y=1,100
            If (IACHAR(D50aOut(y:y))==32) then
                D50aOut(y:100)= D50aOut(y+1:100)
            End If
        End Do
        
        Write (70,777)'    OutRaster = "',TRIM(D50aOut),'"'
        Write (70,'(A)')'    gp.ASCIIToRaster_conversion(InAsciiFile, OutRaster, "FLOAT")'
    End Do
   
    !Adding the TotalErosion File
    Write (TEFile,996) TRIM(dirname),'/','TotalErosion.txt'
    Write (70,777)'    InAsciiFile = "',TRIM(TEFile),'"'
    write (D50aOut,888)  TRIM(dirname),'/Total_Erosion'
    Write (70,777)'    OutRaster = "',TRIM(D50aOut),'"'
    Write (70,'(A)')'    gp.ASCIIToRaster_conversion(InAsciiFile, OutRaster, "FLOAT")'
    
   !Adding the Depth File
    Write (DFile,996) TRIM(dirname),'/','Depth.txt'
    Write (70,777)'    InAsciiFile = "',TRIM(DFile),'"'
    write (D50aOut,888)  TRIM(dirname),'/Depth'
    Write (70,777)'    OutRaster = "',TRIM(D50aOut),'"'
    Write (70,'(A)')'    gp.ASCIIToRaster_conversion(InAsciiFile, OutRaster, "FLOAT")'
    
!        !Adding the Pure Erosion file
!        Write (DFile,996) TRIM(dirname),'/','PureErosion.txt'
!        Write (70,777)'    InAsciiFile = "',TRIM(DFile),'"'
!        write (D50aOut,888)  TRIM(dirname),'/Pure_Erosion'
!        Write (70,777)'    OutRaster = "',TRIM(D50aOut),'"'
!        Write (70,'(A)')'    gp.ASCIIToRaster_conversion(InAsciiFile, OutRaster, "FLOAT")'
    
    !Adding the Final Erosion file
    Write (DFile,996) TRIM(dirname),'/','FinalErosion.txt'
    Write (70,777)'    InAsciiFile = "',TRIM(DFile),'"'
    write (D50aOut,888)  TRIM(dirname),'/Final_Erosion'
    Write (70,777)'    OutRaster = "',TRIM(D50aOut),'"'
    Write (70,'(A)')'    gp.ASCIIToRaster_conversion(InAsciiFile, OutRaster, "FLOAT")'
        
    !Adding the Final NetErosion file
    Write (DFile,996) TRIM(dirname),'/','FinalNetEro.txt'
    Write (70,777)'    InAsciiFile = "',TRIM(DFile),'"'
    write (D50aOut,888)  TRIM(dirname),'/Final_Net_E'
    Write (70,777)'    OutRaster = "',TRIM(D50aOut),'"'
    Write (70,'(A)')'    gp.ASCIIToRaster_conversion(InAsciiFile, OutRaster, "FLOAT")'
    
    !Adding climate change d50 files
9996  Format (a,a,a,F5.2,a)    
    If (ClimatFluc.NE.0) Then        
        Do i=1, ClimateRec
            Write (DFile,9996) TRIM(dirname),'/','d50aL0',IterPrecent(i),'pc.txt'  !!!
            Write (70,777)'    InAsciiFile = "',TRIM(DFile),'"'
            write (D50aOut,996)  TRIM(dirname),'/','d50a',INT(IterPrecent(i)),'pc' !!!
            Write (70,777)'    OutRaster = "',TRIM(D50aOut),'"'
            Write (70,'(A)')'    gp.ASCIIToRaster_conversion(InAsciiFile, OutRaster, "FLOAT")'
        End Do
        
        Do i=1, ClimateRec
            Write (DFile,9996) TRIM(dirname),'/','Depth',IterPrecent(i),'pc.txt'  !!!
            Write (70,777)'    InAsciiFile = "',TRIM(DFile),'"'
            write (D50aOut,996)  TRIM(dirname),'/','Depth',INT(IterPrecent(i)),'pc' !!!
            Write (70,777)'    OutRaster = "',TRIM(D50aOut),'"'
            Write (70,'(A)')'    gp.ASCIIToRaster_conversion(InAsciiFile, OutRaster, "FLOAT")'
        End Do
        
        Do i=1, ClimateRec
            Write (DFile,9996) TRIM(dirname),'/','FinalNetE',IterPrecent(i),'pc.txt'  !!!
            Write (70,777)'    InAsciiFile = "',TRIM(DFile),'"'
            write (D50aOut,996)  TRIM(dirname),'/','FinalNetE',INT(IterPrecent(i)),'pc' !!!
            Write (70,777)'    OutRaster = "',TRIM(D50aOut),'"'
            Write (70,'(A)')'    gp.ASCIIToRaster_conversion(InAsciiFile, OutRaster, "FLOAT")'
        End Do
    End If
    
!Adding final surface grading files
1502 Format (a,a,a,I2,a) 
1503 format (a,a,a,a)
    Do i=1,NumGrade
        if (i<10) then
            write (Name,'(a,I1)')'FinSurfGrad',i
        else 
            write (Name,'(a,I2)')'FinSurfGrad',i
        End if
        Write (DFile,1502) TRIM(dirname),'/','FinSurfGrad',i,'.txt'  
        Write (70,777)'    InAsciiFile = "',TRIM(DFile),'"'
        write (D50aOut,1503)  TRIM(dirname),'/',TRIM(Name)
        Write (70,777)'    OutRaster = "',TRIM(D50aOut),'"'
        Write (70,'(A)')'    gp.ASCIIToRaster_conversion(InAsciiFile, OutRaster, "FLOAT")'
    End Do
         
    Write (70,'(A)')'except:'
    Write (70,*)'   print gp.GetMessages()'
    Close (70, status='keep') 
End Subroutine PythonScriptWriter
    
    
!########################################################################################
!This subroutine calculate discharge, SizeTH, and minsize; account for climate changes.
!Added in 7/10/08 version 4.7
!########################################################################################        
  Subroutine Discharge (tmpQf, FlowAcc, CellSize, Slope, GradingSize,NumGrade, QdNew, SizeThNew, DminNew,TempLayerGrad,Numlayers,LayerDepth,ArmourDepth)
        Real(8) FlowAcc,tmpQf,QdNew, SizeThNew, DminNew, CellSize, Slope, GradingSize(NumGrade),TempLayerGrad(Numlayers,NumGrade), tmpDepth,LayerDepth,ArmourDepth 
        Integer ii,Numlayers

        tmpDepth=(NumLayers*LayerDepth)-(sum((TempLayerGrad(:,NumGrade)/100)*LayerDepth))
     !   QdNew = tmpQf*(FlowAcc/CellSize**2)**0.5 !!!!discharge per unit width [m^3/iter/m]
      !  QdNew = tmpQf*(FlowAcc/CellSize**2)**0.1
      QdNew = tmpQf*(FlowAcc/CellSize**2)
        tmpDepth=tmpDepth+ArmourDepth
        !QdNew = (tmpQf*(FlowAcc)/CellSize**2)/((tmpDepth/100)*10)!calibrated to 10cm
     !   QdNew = tmpQf*((FlowAcc/CellSize**2)**0.5)*(((tmpDepth/100)*10)**-0.5)!calibrated for 10cm
     !   QdNew = tmpQf*((FlowAcc/CellSize**2)**0.5)*(((tmpDepth/100)*5)**-0.5)!calibrated for 20cm
        !QdNew = tmpQf*((FlowAcc/CellSize**2)**0.5)*(((tmpDepth/100)*100)**-0.5)!calibrated for 1cm
        SizeThNew = ((((0.1*QdNew)**0.6) *(Slope**0.7))/0.074)*1000 !*1000 converet to mm
    
        Do ii=1, NumGrade
            if (GradingSize(ii).LE.SizeThNew) then
                DminNew=ii
            End if
        end Do                               
    End Subroutine Discharge

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!   Function to generate random number between 0.0 and 1.0
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
Real*4 Function ran2(idum)
!  random number from Numerical recipes RAN2
 ! integer idum
  PARAMETER (m=714025, ia=1366, ic=150889, rm=1./m)
  DIMENSION ir(97)
  SAVE 
  DATA iff /0/
  IF (idum.lt.0.or.iff.eq.0) THEN
    iff=1
    idum=mod(ic-idum,m)
  !  idum=mod(ic-idum,m)
    DO 11 j=1,97
     idum=mod(ia*idum+ic,m)
     ir(j)=idum
11     CONTINUE
    idum=mod(ia*idum+ic,m)
    iy=idum
  END IF
  j=1+(97*iy)/m
  IF (j.gt.97.or.j.lt.1) THEN
    WRITE(*,*) 'error in RAN',j
    STOP
  END IF
  iy=ir(j)
  ran2=iy*rm
  idum=mod(ia*idum+ic,m)
  ir(j)=idum
  RETURN
END
!########################################################################################
!This subroutine calculate the rate of soil creep at the surface anjd the relative change
! with depth. Added in Version 5.6 17/5/2010
!########################################################################################        
Subroutine CreepRateCalc(CreepRate,NumLayers,LayerDepth,ArmourDepth,Slope,CreepDepth,SurfCreepProp,SurfCreepRate,NCols,NRows,dt,CellSize,xllcor,yllcor)
    Real(8) CreepRate,LayerDepth,ArmourDepth,Slope(NCols,NRows),CreepDepth(NumLayers+1),SurfCreepRate(NCols,NRows),Depth,SurfCreepProp(NCols,NRows)
    Real(8) CellSize,xllcor,yllcor 
    Integer i,NumLayers,dt
    !The surface creep rate (ks[cm/year])is the maximum rate; Creep rate for a slope of 20%
    CreepRate=CreepRate/(10000) !Convert from cm/year to cm/iteration. 1year=10iteration*dt
    SurfCreepRate=((Slope/0.2)*CreepRate)*dt !was devided by 10,000; Creep rate for a slope of 20% (0.2)
	!SurfCreepRate=(((Slope/0.2)**0.1)*CreepRate)*dt
    SurfCreepProp=SurfCreepRate/(CellSize*100) !the proportion of layer displaced; Convert pixel size to cm
    Do i=1,NumLayers+1 !Calculate relative change in creep rate with depth
        If (i==1) Then
            Depth=0.5
        Else
            Depth=(i-1)*LayerDepth-(LayerDepth/2)
        End If
        CreepDepth(i)=1*EXP(-0.02*i*Depth)
    End Do
1001    Format (<NCols>F10.3)
    OPEN(unit=94339,file='SurfaceCreepRate.txt',status='unknown')
    Write (94339,*) 'ncols ',NCols
    Write (94339,*) 'nrows ',NRows
    Write (94339,*) 'xllcorner ',xllcor 
    Write (94339,*) 'yllcorner ',yllcor 
    Write (94339,*) 'cellsize  ',CellSize 
    Write (94339,*) 'NODATA_value  -9999'  
    Write (94339,1001) SurfCreepRate
    Close (unit=94339, status='keep')
End Subroutine CreepRateCalc