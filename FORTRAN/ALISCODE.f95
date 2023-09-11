!                                          **********************************************************************
!                                          **                                                                  **
!                                          **                  A Code for Theta - Phi Distribution             **
!                                          **                                                                  **
!                                          **********************************************************************

Program Theta_Phi_Distribution
    
    Implicit None

!    Integer:: iseed
    Integer:: N_sshot,Natom_max,Nwater_max,Norg_max
	Integer:: N_counter,N_count
	Integer:: i,ii,m,mm,j,jj,kk
	Integer:: np,nt,nz,tot_nz
	Integer:: Mol_No,Atom_No
	Integer:: RDF_Switch
	
	Real(8):: del_z,Z_Lower,Z_Upper
	Real(8):: mass_O,mass_H
	Real(8):: arg,pphi,ptheta,V1,V2,V3
	Real(8):: Tot_water
	Real(8):: Delt_X,Delt_Y,Delt_Z
	
    Character(20):: nfile1,nfile2,nfile3,nfile4,nfile5
    Character(20):: nfile6,nfile7,nfile8,nfile9,nfile10
	Character(20):: nfile11,nfile12,nfile13,nfile14,nfile15
	Character(30):: Gromax_Name
    
    Integer,Allocatable:: Unit_No(:),N_atom(:),N_water(:)
	
!	Integer,Allocatable:: Mol_No(:,:),Atom_No(:,:)
	Character(3),Allocatable:: Mol_Type(:,:),Atom_Type(:,:)
	Real(8),Allocatable:: X_atom(:,:),Y_atom(:,:),Z_atom(:,:)
	
	Real(8),Allocatable:: XO(:,:),YO(:,:),ZO(:,:)
    Real(8),Allocatable:: XH1(:,:),YH1(:,:),ZH1(:,:)
    Real(8),Allocatable:: XH2(:,:),YH2(:,:),ZH2(:,:)
	Real(8),Allocatable:: X_cm(:,:),Y_cm(:,:),Z_cm(:,:)
    Real(8),Allocatable:: mu_X(:,:),mu_Y(:,:),mu_Z(:,:)
	Real(8),Allocatable:: Mu_2(:,:)
	Real(8),Allocatable:: Theta(:,:),Phi(:,:)
	
	Integer,Allocatable:: N_H2O(:,:)
    Integer,Allocatable:: ntheta(:,:),nphi(:,:)
	
	Real(8),Allocatable:: Tot_H2O(:)

!    Call Srand(16)
!    nxsite(1)=0;    nxsite(2)=0;   nxsite(3)=1;   nxsite(4)=1

!	Open (99,File='nohup.out',Status='unknown')
!	Read (99,*) mcell

! Open the input file.
	Open (Unit=88, File='co.input', Status='Old')
	Open (Unit=50, File='Theta_Distribution' , Status='Unknown')
	Open (Unit=51, File='Phi_Distribution' , Status='Unknown')
	
	Read (88,9991) nfile1
	Read (88,9991) nfile2
	Read (88,9991) nfile3	
	Read (88,9991) nfile4	
	Read (88,9991) nfile5
	Read (88,9991) nfile6	
	Read (88,9991) nfile7
	Read (88,9991) nfile8	
	Read (88,9991) nfile9
	Read (88,9991) nfile10
	Read (88,9991) nfile11
	Read (88,9991) nfile12
	Read (88,9991) nfile13
	Read (88,9991) nfile14
	Read (88,9991) nfile15
9991	Format (a)

    Print*,nfile1
	Print*,nfile2
	Print*,nfile3
	Print*,nfile4
	Print*,nfile5
	Print*,nfile6
	Print*,nfile7
	Print*,nfile8
	Print*,nfile9
	Print*,nfile10
	Print*,nfile11
	Print*,nfile12
	Print*,nfile13
	Print*,nfile14
	Print*,nfile15

    Read (88,9992) N_sshot
9992 Format (I2)

    Read (88,9993),del_z,Z_Lower,Z_Upper
9993 Format (3F5.2)

    Read (88,9994),mass_O,mass_H
9994 Format (F6.3,F6.4)
 
    Read (88,9992) RDF_Switch

! Allocate the arrays that only depend on N_sshot.
    Allocate ( Unit_No(N_sshot),N_atom(N_sshot),N_water(N_sshot),N_org(N_Sshot) )
	
	Do mm=1,N_sshot
	     Unit_No(mm)=mm
	End Do
	
! Open all the Gromx snapshot files.
    Open (Unit=Unit_No(1), File=nfile1 , Status='UNKNOWN')
	Open (Unit=Unit_No(2), File=nfile2 , Status='UNKNOWN')
	Open (Unit=Unit_No(3), File=nfile3 , Status='UNKNOWN')
	Open (Unit=Unit_No(4), File=nfile4 , Status='UNKNOWN')
	Open (Unit=Unit_No(5), File=nfile5 , Status='UNKNOWN')
	Open (Unit=Unit_No(6), File=nfile6 , Status='UNKNOWN')
	Open (Unit=Unit_No(7), File=nfile7 , Status='UNKNOWN')
	Open (Unit=Unit_No(8), File=nfile8 , Status='UNKNOWN')
	Open (Unit=Unit_No(9), File=nfile9 , Status='UNKNOWN')
	Open (Unit=Unit_No(10), File=nfile10, Status='UNKNOWN')
    Open (Unit=Unit_No(11), File=nfile11, Status='UNKNOWN')
    Open (Unit=Unit_No(12), File=nfile12 , Status='UNKNOWN')
    Open (Unit=Unit_No(13), File=nfile13 , Status='UNKNOWN')
    Open (Unit=Unit_No(14), File=nfile14 , Status='UNKNOWN')
	Open (Unit=Unit_No(15), File=nfile15 , Status='UNKNOWN')
        
!   Read (Unit_No(1),*),Gromax_Name   ;   Read (Unit_No(1),*),N_atom(Unit_No(1))

!	Read (Unit_No(2),*),Gromax_Name   ;   Read (Unit_No(2),*),N_atom(Unit_No(2))

!	Read (Unit_No(3),*),Gromax_Name   ;   Read (Unit_No(3),*),N_atom(Unit_No(3))

!	Read (Unit_No(4),*),Gromax_Name   ;   Read (Unit_No(4),*),N_atom(Unit_No(4))
	
	Do ii=1,N_sshot
	     Read (Unit_No(ii),*),Gromax_Name
		 Read (Unit_No(ii),*),N_atom(Unit_No(ii))
	End Do

! Before allocation of the rest, calculate the maximum number of atoms in the snapshots.
    Natom_max=0
	Do mm=1,N_sshot
	     IF ( mm==1 ) Then
		     Natom_max=N_atom(mm)
		 Else 
		     IF ( N_atom(mm) > Natom_max ) Then
			     Natom_max=N_atom(mm)
			 End IF
		 End IF
	End Do

!	Allocate ( Mol_No(N_sshot,Natom_max),Atom_No(N_sshot,Natom_max) )
	Allocate ( Mol_Type(N_sshot,Natom_max),Atom_Type(N_sshot,Natom_max) )
	Allocate ( X_atom(N_sshot,Natom_max),Y_atom(N_sshot,Natom_max),Z_atom(N_sshot,Natom_max) )

! Calculation of the total number of water molecules in each snapshot.
	N_water(:)=0
	N_org(:)=0
	Do mm=1,N_sshot
	     N_counter=0
	     Do ii=1,N_atom(mm)
		 
              Read (mm,2001),Mol_No,Mol_Type(mm,ii),Atom_Type(mm,ii),Atom_No, &
	          &                        X_atom(mm,ii),Y_atom(mm,ii),Z_atom(mm,ii),V1,V2,V3                 
2001      Format (I5,A3,4x,A3,I5,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.4,1x,f7.4,1x,f7.4)

              IF (  Mol_Type(mm,ii)=="SOL" ) Then
		           N_counter=N_counter+1
			  End IF
			  IF ( N_counter==3 ) Then
			      N_counter=0
				  N_water(mm)=N_water(mm)+1
			  End IF
			  
			  IF ( Mol_Type(mm,ii)=="UNK" ) Then
		           N_count=N_count+1
			  End IF
			  IF ( N_count==20 ) Then
			      N_count=0
				  N_org(mm)=N_org(mm)+1
			  End IF
			  
         End Do
		 Read (mm,2008),Lbox_x(mm),Lbox_y(mm),Lbox_z(mm)                
2008      Format (3f10.5)
		 Print*,'mm=',mm,N_water(mm),N_org(mm)
	End Do

! Before allocation of the rest, calculate the maximum number of watrs and orgs in the snapshots.

	Nwater_max=0
	Norg_max=0
	Do mm=1,N_sshot
	     IF ( mm==1 ) Then
	         Nwater_max=N_water(mm)
			 Norg_max=N_org(mm)
		 Else 
		     IF ( N_water(mm) > Nwater_max ) Then
			     Nwater_max=N_water(mm)
			 End IF
			 IF ( N_org(mm) > Norg_max ) Then
			     Norg_max=N_org(mm)
			 End IF
		 End IF
	End Do
!	Print*,'Nwater_max=',Nwater_max
    
    Allocate ( XO(N_sshot,Natom_max),YO(N_sshot,Natom_max),ZO(N_sshot,Natom_max) )
    Allocate ( XH1(N_sshot,Natom_max),YH1(N_sshot,Natom_max),ZH1(N_sshot,Natom_max) )
    Allocate ( XH2(N_sshot,Natom_max),YH2(N_sshot,Natom_max),ZH2(N_sshot,Natom_max) )
	Allocate ( X_cm(N_sshot,Nwater_Max),Y_cm(N_sshot,Nwater_Max),Z_cm(N_sshot,Nwater_Max) )
	
! Decide weather calculate RDF or not. Thisis done just to save some memory and be able to run the program (hopefully!).
! RDF_Swittch==1 means that RDF calculations are set to be implemented.

	IF ( RDF_Switch == 0 ) Then
	     Allocate ( Mu_X(N_sshot,Nwater_Max),Mu_Y(N_sshot,Nwater_Max),Mu_Z(N_sshot,Nwater_Max) )
	     Allocate ( Mu_2(N_sshot,Nwater_Max) )
	     Allocate ( Theta(N_sshot,Nwater_Max),Phi(N_sshot,Nwater_Max) )
	Else IF ( RDF_Switch == 1 ) Then
	     Allocate ( R_WOr(Nwater_max,Norg_max),R_WW(Nwater_max,Nwater_max) )
		 Allocate ( R_OrW(Nwater_max,Norg_max),R_OrOr(Nwater_max,Nwater_max) )
	End IF
	
	Print*,"Checked ............... #1"

!		Read (88,9993),box_x,box_y,box_z
!9993 Format (F5.2,F5.2,F5.2)

    N_water(:)=0
    Do mm=1,N_sshot
	     N_counter=0
         Do ii=1,N_atom(mm)
	          IF (  Mol_Type(mm,ii)=="SOL" ) Then
		          N_counter=N_counter+1
		          IF (  Atom_Type(mm,ii)==" OW" ) Then 
			          XO(mm,ii)=X_atom(mm,ii)
				      YO(mm,ii)=Y_atom(mm,ii)
				      ZO(mm,ii)=Z_atom(mm,ii)	  
		          Else IF (  Atom_Type(mm,ii)=="HW1" ) Then 
			          XH1(mm,ii)=X_atom(mm,ii)
				      YH1(mm,ii)=Y_atom(mm,ii)
				      ZH1(mm,ii)=Z_atom(mm,ii)	  
		          Else IF (  Atom_Type(mm,ii)=="HW2" ) Then 
			          XH2(mm,ii)=X_atom(mm,ii)
				      YH2(mm,ii)=Y_atom(mm,ii)
				      ZH2(mm,ii)=Z_atom(mm,ii)	  
		          End IF

                  IF ( N_counter==3 ) Then
				      N_counter=0
! N_water(mm) is the total number of water molecules in each snapshot mm.
			          N_water(mm)=N_water(mm)+1
! N_water(mm) here is just an index to identify which water molecule in the snapshot.				  
				      X_cm(mm,N_water(mm))=(mass_O*XO(mm,ii-2)+mass_H*XH1(mm,ii-1) &
				      &                                       + mass_H*XH2(mm,ii))                                                   &
				      &                                        / (mass_O+2.0D0*mass_H)
				  
				      Y_cm(mm,N_water(mm))=(mass_O*YO(mm,ii-2)+mass_H*YH1(mm,ii-1) &
				      &                                       + mass_H*YH2(mm,ii))                                                   &
				      &                                        / (mass_O+2.0D0*mass_H)
				  
				      Z_cm(mm,N_water(mm))=(mass_O*ZO(mm,ii-2)+mass_H*ZH1(mm,ii-1) &
				      &                                       + mass_H*ZH2(mm,ii))                                                   &
				      &                                        / (mass_O+2.0D0*mass_H)
				  
				      IF ( RDF_Switch==0 ) Then
                          Mu_X(mm,N_water(mm))=2.0D0*XO(mm,ii-2)-(XH1(mm,ii-1)+XH2(mm,ii))
			              Mu_Y(mm,N_water(mm))=2.0D0*YO(mm,ii-2)-(YH1(mm,ii-1)+YH2(mm,ii))
			              Mu_Z(mm,N_water(mm))=2.0D0*ZO(mm,ii-2)-(ZH1(mm,ii-1)+ZH2(mm,ii))

			              Mu_2(mm,N_water(mm))=(mu_X(mm,N_water(mm)))**2+(mu_Y(mm,N_water(mm)))**2 &
					      &                                     +(mu_Z(mm,N_water(mm)))**2
                          Theta(mm,N_water(mm))=ACOS(mu_Z(mm,N_water(mm))/SQRT(Mu_2(mm,N_water(mm))))/0.017453292D0

                          IF ( ( ABS(Mu_X(mm,N_water(mm))) <= 10D-7 ) .AND. ( ABS(Mu_Y(mm,N_water(mm))) <= 10D-7 ) ) Then
                              Phi(mm,N_water(mm))=180.0D0*(2.0D0*Rand(0)-1.0D0)
                          Else    
                              arg=Mu_X(mm,N_water(mm))/SQRT(mu_X(mm,N_water(mm))**2+mu_Y(mm,N_water(mm))**2)
! To convert from radian to degrees multiply by 180/pi.
                              Phi(mm,N_water(mm))=ACOS(arg)*SIGN(1.0D0,mu_Y(mm,N_water(mm)))/0.017453292D0
                          End IF
                      End IF					  
			      End IF		  
		      End IF
		 End Do
		 Print*,'mm=',mm,N_water(mm)
	End Do
							
	IF ( RDF_Switch==0 ) Then
!	Deallocate (XO,YO,ZO,XH1,YH1,ZH1,XH2,YH2,ZH2,X_atom,Y_atom,Z_atom,Mu_X,Mu_Y,Mu_Z,Mu_2)	
	    Deallocate (Mu_X,Mu_Y,Mu_Z,Mu_2)
	End IF
	
	
!   ********************************************************************************
!   *                                                                                                           *
!   *                                         RDF Calculation                                         *
!   *                                                                                                           *
!   ********************************************************************************

!  I didn't calculate the centre of mass of organic molecules. Because that will increase the number of functions and 
!  that will cause more memory issues. So for now I'm assuming that we have all the centres of mass calculated before 
!  this part of the calculations.

! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! ^     STEP 1: Calculation of the distances between all particles, Rij       ^
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! We can claculate all the distances between all pairs of particles, even between atoms. Here I am just calculating the
! distances between molecules (between centres of masses of the water and organic molecules). But for these, we need to also 
! calculate the centres of mass of the organic molecules as well, but for the moment I am going to skip these calculations 
! and use the variables X_cm_org, Y_cm_org, and Z_cm_org for the coordinates of the centres of mass of organic molecules.

        Do mm=1,N_sshot
! N_water(mm) and N_org(mm) here are the number of water and org molecules in the whole layer of each snapshot "mm".
! We can consider them as the indeces of, in turn, water and org molecules in the whole layer of each snapshot "mm",
! which means that kk and jj are the indeces of, in turn, water and org molecules in the snapshot.
	         Do kk=1,N_water(mm)
		           Do jj=1,N_org(mm)
				        Delt_X= X_cm(mm,kk) - X_cm_org(mm,jj)
					    Delt_Y= Y_cm(mm,kk) - Y_cm_org(mm,jj)
					    Delt_Z= Z_cm(mm,kk) - Z_cm_org(mm,jj)
! First make sure that the distance doesn't exceed the simulation box length. This is done to avoid problems with some of the molecules.
					    Delt_X=Delt_X - LBox_x*ANINT(Delt_X/LBox_x)
					    Delt_Y=Delt_Y - LBox_y*ANINT(Delt_Y/LBox_y)
					    Delt_Z=Delt_Z - LBox_z*ANINT(Delt_Z/LBox_z)
					   
! It's better to store the distances this way (by the snapshot) because these will become handy in the future calculations.
                        R_WOr(mm,kk,jj) = ABS( SQRT ( Delt_X**2 + Delt_Y**2 + Delt_Z**2 ) )
                        R_OrW(mm,jj,kk) = R_WOr(mm,kk,jj)
				   End Do
			 End Do
			 
			 Do kk=1,N_water(mm)-1
		           Do jj=kk+1,N_water(mm)
				        Delt_X= X_cm(mm,kk) - X_cm(mm,jj)
					    Delt_Y= Y_cm(mm,kk) - Y_cm(mm,jj)
					    Delt_Z= Z_cm(mm,kk) - Z_cm(mm,jj)

					    Delt_X=Delt_X - LBox_x*ANINT(Delt_X/LBox_x)
					    Delt_Y=Delt_Y - LBox_y*ANINT(Delt_Y/LBox_y)
					    Delt_Z=Delt_Z - LBox_z*ANINT(Delt_Z/LBox_z)
				       
                        R_WW(mm,kk,jj) = ABS( SQRT ( Delt_X**2 + Delt_Y**2 + Delt_Z**2  ) )
				   End Do
			 End Do
			 
			 Do kk=1,N_org(mm)-1
		           Do jj=kk+1,N_org(mm)
				        Delt_X= X_cm_org(mm,kk) - X_cm_org(mm,jj)
					    Delt_Y= Y_cm_org(mm,kk) - Y_cm_org(mm,jj)
					    Delt_Z= Z_cm_org(mm,kk) - Z_cm_org(mm,jj)

					    Delt_X=Delt_X - LBox_x*ANINT(Delt_X/LBox_x)
					    Delt_Y=Delt_Y - LBox_y*ANINT(Delt_Y/LBox_y)
					    Delt_Z=Delt_Z - LBox_z*ANINT(Delt_Z/LBox_z)
				  
                        R_OrOr(mm,kk,jj) = ABS( SQRT ( Delt_X**2 + Delt_Y**2 + Delt_Z**2 ) )
				   End Do
			 End Do
			 
		End Do
			 
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! ^       STEP 2: Calculation of RDF functions in both Cases 1 and 2       ^
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
			 
! Now we calculate the RDF functions in two different cases/scenarios: 
! Case 1: When the whole layer is considered. All we need to do is to input del_z = Z_Upper - Z_Lower.
! This only gives us one layer, which is the whole water/organic layer. In this case, the RDF calculation is three-dimensional.

! Case 2: When the water/organic layer is sliced into sublayers. In this case, the RDF calculation is two-dimensional. To make sure that
! we are in this case, it is best for del_z to not exceed 4 angstrom (but to be more accurate, it should't exceed 2 angstrom).

! Otehr Cases: It is very comlicated to do RDF calculations for cases in between two- and three-dimensional regimes.

! Now let's caluclate them together, but we separate them whenever necessary.

! Total number of bins/sublayers/slices, "nz", in the snapshots (this number is the same in all snapshots).
! Also, del_z is the thickness of each slice.
! Note: In Case 1, del_z = Z_Upper - Z_Lower, and hence, tot_nz=1.
	tot_nz=IDINT((Z_Upper - Z_Lower)/del_z)

! First calculate the total number of water and organic molecules in each slice. 
! Once again, N_water(mm) is the total number of water molecules in the whole simulation box of each snapshot "mm". 
! Also, N_org(mm) is the total number of organic molecules in the whole simulation box of each snapshot "mm".

! But, Total_water(nz) and Total_org(nz) are the average of the total number of water and organic molecules 
! in each slice nz (for all snapshots), respectively.
 
! Note: In Case 1, tot_nz=1 and there is only 1 possible nz value (nz=1). So Total_water(nz) should be equal to N_water(mm), 
! also Total_org(nz)=N_org(mm).
	    
		Total_water(:)=0
        Do mm=1,N_sshot
		     Do kk=1,N_water(mm)
                  nz=IDINT((1.0D0/del_z)*(Z_cm(mm,kk) - Z_Lower)) + 1
		          Total_water(nz)=Total_water(nz)+1
			 End Do
		End Do
		Total_water(nz)=IDINT(Total_water(nz)/N_sshot)
		
		Total_org(:)=0
        Do mm=1,N_sshot
		     Do jj=1,N_org(mm)
                  nz=IDINT((1.0D0/del_z)*(Z_cm(mm,jj) - Z_Lower)) + 1
		          Total_org(nz)=Total_org(nz)+1
			 End Do
		End Do
		Total_org(nz)=IDINT(Total_org(nz)/N_sshot)
	
! Hist gives the number of particles in each radial shell of each sublayer/slice of the system.
! Note: In case 1, there is only one sublayer which is equal to the whole layer.

        Hist_WOr (:,:)=0    ;    Hist_OrW(:,:)=0
	    Hist_WW(:,:)=0     ;    Hist_OrOr (:,:)=0
		
        Do mm=1,N_sshot
	         Do kk=1,N_water(mm)
		           Do jj=1,N_org(mm)

! Calculate the number of pairs in each histogram/shell in the snapshot before moving to the next snapshot.
! "nz" is the index of sublayer/slice for the central (or first) molecule in the simulation box,
! and "nshell" is the index of the spherical shell around the central molecule.
! In Case 1, there is only one value for "nz" and it is equal to 1 (i.e., nz=1).
                        nz=IDINT((1.0D0/del_z)*(Z_cm(mm,kk) - Z_Lower)) + 1
! This next line may not be necessary, but it saves time in the calculations.
		                IF ( (nz>tot_nz) .OR. (nz<=0) ) Cycle
					   
! Before counting the number of pairs in each histogram, first determine if the second molecule also lies within 
! this same sublayer/slice.
! "del_z" is the thickness of each slice.
! In Case 1: del_z = Z_Upper - Z_Lower (or whole layer), and "nz_second" is equal to 1 as well.
                        nz_second=IDINT((1.0D0/del_z)*(Z_cm_org(mm,jj) - Z_Lower)) + 1
						
! The next if statement is a bit dangerous. Make sure that both "nz" and "nz_second" are integer.
! "DelR" is the thickness of the spherical shell around the central atom.
                        IF ( nz == nz_second ) Then
                            IF (   (R_WOr(mm,kk,jj) < LBox_x/2.0D0) .AND. (R_WOr(mm,kk,jj) < LBox_y/2.0D0)  &
							& .AND. (R_WOr(mm,kk,jj) < LBox_z/2.0D0)    ) Then
	                            nshell=IDINT( R_WOr(mm,kk,jj) / DelR ) + 1
! Now store the pair in each histogram of each sublayer.
		                        Hist_WOr(nz,nshell)=Hist_WOr(nz,nshell)+1
						        Hist_OrW(nz,nshell)=Hist_OrW(nz,nshell)+1
							End IF
					    End IF
					   
                   End Do
              End Do

              Do kk=1,N_water(mm)-1
		           Do jj=kk+1,N_water(mm)
                        nz=IDINT((1.0D0/del_z)*(Z_cm(mm,kk) - Z_Lower)) + 1
! This next line may not be necessary, but it saves time in the calculations.
		                IF ( (nz>tot_nz) .OR. (nz<=0) ) Cycle
                        nz_second=IDINT((1.0D0/del_z)*(Z_cm(mm,jj) - Z_Lower)) + 1

                        IF ( nz == nz_second ) Then
						     IF (   (R_WW(mm,kk,jj) < LBox_x/2.0D0) .AND. (R_WW(mm,kk,jj) < LBox_y/2.0D0)  &
                            & .AND. (R_WW(mm,kk,jj) < LBox_z/2.0D0)   ) Then
	                             nshell=IDINT( R_WW(mm,kk,jj) / DelR ) + 1
		                         Hist_WW(nz,nshell)=Hist_WW(nz,nshell)+2
							 End IF
	                    End IF
					   
                   End Do
              End Do
		  
		      Do kk=1,N_org(mm)-1
		           Do jj=kk+1,N_org(mm)
                        nz=IDINT((1.0D0/del_z)*(Z_cm_org(mm,kk) - Z_Lower)) + 1
! This next line may not be necessary, but it saves time in the calculations.
		                IF ( (nz>tot_nz) .OR. (nz<=0) ) Cycle
                        nz_second=IDINT((1.0D0/del_z)*(Z_cm_org(mm,jj) - Z_Lower)) + 1
						
                        IF ( nz == nz_second ) Then
						    IF (   (R_OrOr(mm,kk,jj) < LBox_x/2.0D0) .AND. (R_OrOr(mm,kk,jj) < LBox_y/2.0D0)  &
							& .AND. (R_OrOr(mm,kk,jj) < LBox_z/2.0D0)    ) Then
	                            nshell=IDINT( R_OrOr(mm,kk,jj) / DelR ) + 1
		                        Hist_OrOr(nz,nshell)=Hist_OrOr(nz,nshell)+2
							End IF
	                    End IF
					   
                   End Do
              End Do
			 
	    End Do
		
! Now take the average of the number of pairs in each histogram for all snapshots.
! Again, in Case 1, nz=1.
! Note: Here we must use IDNINT (nearest Integer in double precision).
        Hist_WOr(nz,nshell)=IDNINT(Hist_WOr(nz,nshell)/N_sshot)
		Hist_OrW(nz,nshell)=IDNINT(Hist_OrW(nz,nshell)/N_sshot)
		Hist_WW(nz,nshell)=IDNINT(Hist_WW(nz,nshell)/N_sshot)
		Hist_OrOr(nz,nshell)=IDNINT(Hist_OrOr(nz,nshell)/N_sshot)

! Now calculate the Radial Distribution Function in each sublayer based on the code described in Allen's book/GithHub code.
! To calculate the total number of shells go with the smallest dimension (for example, LBox_x in this system).
! IMPORTANT: We need to make distinguish between Cases 1 and 2. Because here the "V_shell" and "V_sublayer" are basically 
! defined based on if we are looking at the whole layer (Case 1, three-dimensional regime) or sublayers (Case 2, two-dimensional regime).
	   Tot_nshell= (LBox_x/2.0D0) / DelR
	   PI=3.14159265359D0
	   
	   Do nz = 1,tot_nz 
	        Do nshell = 1,Tot_nshell
			
			      R_RDF(nz,nshell) = DelR * (nshell-0.5)
				  
! We can easily differentiate between Cases 1 and 2 by using the following IF statement.
! IF the difference between them is almost zero ( or del_z = (Z_Upper - Z_Lower) ) Then ...
				  IF ( (del_z - (Z_Upper - Z_Lower)) < 1.0D-4 ) Then
! Case 1: we are in the three-dimensional regime.		      
				      V_shell = (4.0D0/3.0D0) * PI * ( ( (nshell+1) * DelR )**3 - ( nshell * DelR )**3 )
					  V_sublayer = Lbox_x * Lbox_y * del_z
				  Else 
! Case 2: we are in the two-dimensional regime.
                      V_shell = PI * ( ( (nshell+1) * DelR )**2 - ( nshell * DelR )**2 )
					  V_sublayer = Lbox_x * Lbox_y
                  Else IF

				  g_WOr(nz,R_RDF(nz,nshell)) = (Hist_WOr(nz,nshell) / V_shell) / ( Total_org(nz) / V_sublayer ) / Total_water(nz)
				  
				  g_OrW(nz,R_RDF(nz,nshell)) = (Hist_OrW(nz,nshell) / V_shell) / ( Total_water(nz) / V_sublayer ) / Total_org(nz)
				  
				  g_WW(nz,R_RDF(nz,nshell)) = (Hist_WW(nz,nshell) / V_shell) / ( Total_water(nz) / V_sublayer ) / Total_water(nz)
				  
				  g_OrOr(nz,R_RDF(nz,nshell)) = (Hist_OrOr(nz,nshell) / V_shell) / ( Total_org(nz) / V_sublayer ) / Total_org(nz)
				  
! Write every one of these new functions in a separate file with this format:
!                 Write (31,*) R_RDF(nz,nshell) , g_WOr(nz,R_RDF(nz,nshell))
! *               Format (2f16.7)
!                 Write (32,*) R_RDF(nz,nshell) , g_OrW(nz,R_RDF(nz,nshell))
!  etc.
       		
			End Do
	   End Do
	   
	End IF
	
!   ********************************************************************************
!   *                                                                                                           *
!   *                                   End of RDF Calculation                                    *
!   *                                                                                                           *
!   ********************************************************************************

	
	Allocate ( N_H2O(N_sshot,tot_nz) )
	
	IF ( RDF_Switch==0 ) Then
        Allocate ( ntheta(tot_nz,0:180),nphi(tot_nz,-180:180) )	
	
	    Print*,"Checked ............... #2"

	    N_H2O(1:N_sshot,1:tot_nz)=0
	    Do mm=1,N_sshot
	         Do nz=1,tot_nz
!		          Print*,'mm=',mm,'nz=',nz,'N_H2O=',N_H2O(mm,nz)
		     End Do
	    End Do
	
		Print*,''
		Print*,'******************************************'
		Print*,'******************************************'
		Print*,''
	
		Do ii=1,tot_nz
			 Do jj=0,180
				  ntheta(ii,jj)=0
!		          Print*,'mm=',mm,'nz=',nz,'N_H2O=',N_H2O(mm,nz)
			 End Do
		End Do
	
		Do ii=1,tot_nz
	         Do jj=-180,180
		          nphi(ii,jj)=0
!		          Print*,'mm=',mm,'nz=',nz,'N_H2O=',N_H2O(mm,nz)
		     End Do
		End Do
	
	Do mm=1,N_sshot
!	     N_H2O(mm,1:tot_nz)=0
	     Do kk=1,N_water(mm)
!		       IF ( kk==1 ) N_H2O(mm,1:tot_nz)=0
	           nz=IDINT((1.0D0/del_z)*(Z_cm(mm,kk) - Z_Lower)+1)
		       IF ( (nz>tot_nz) .OR. (nz<=0) ) Cycle
! Total number of water molecules in each bin of each snapshot.
               N_H2O(mm,nz)=N_H2O(mm,nz)+1
! Define a new theta and phi for the water molecules in each nz bin. N_H2O(mm,nz) here 
! is just an index to identify which water molecule is in the bin.
! Count the number of water molecules having a specific theta and phi in each bin of each snapshot.
               nt=IDNINT(Theta(mm,kk))
			   np=IDNINT(Phi(mm,kk))
			   ntheta(nz,nt)=ntheta(nz,nt)+1
			   nphi(nz,np)=nphi(nz,np)+1
         End Do
	End Do
	
	Allocate ( Tot_H2O(tot_nz) )
	
	Print*,"Checked ............... #3"

! The average of total water molecules in each bin.
	Do kk=1,tot_nz
	     Tot_H2O(kk)=0.0D0
	     Do mm=1,N_sshot
! The sum of all H2O molecules in each bin of all snapshots.
		        Tot_H2O(kk)=Tot_H2O(kk)+FLOAT(N_H2O(mm,kk))
		  End Do
		  Tot_H2O(kk)=Tot_H2O(kk)/N_sshot
	End Do
	
	Print*,"Checked ............... #4"
	
! Calculate the average theta and phi probability distribution in each bin among all snapshots.
    Tot_water=0.0D0
	Do kk=1,tot_nz
	     Write(50,441),kk,Tot_H2O(kk)
 441   Format("Bin Number=",I5,F15.8)
          Tot_water=Tot_water+Tot_H2O(kk)
	     Do ii=0,180
!               IF ( ntheta(kk,ii)==0 ) Cycle
               ptheta=FLOAT(ntheta(kk,ii))/(Tot_H2O(kk)*FLOAT(N_sshot))
               Write(50,636),FLOAT(ii),ptheta
636          Format(f7.3,f15.8)
        End Do
		
		Write(51,441),kk,Tot_H2O(kk)
		Do jj=-180,180
!           IF ( nphi(kk,jj)==0 ) Cycle
           pphi=FLOAT(nphi(kk,jj))/(Tot_H2O(kk)*FLOAT(N_sshot))
           Write(51,637),FLOAT(jj),pphi
637          Format(f8.3,f15.8)
        End Do
	End Do
	
	End IF
	
	Write(50,488) Tot_water
	Write(51,488) Tot_water
488  Format ('Total Number of Water Molecules=',F15.8)

    Print*,'Total Number of Water Molecules=',Tot_water
	
	Print*,"End of Program"
		
End Program Theta_Phi_Distribution
