!     last update 17.1.2011
      module rdens

      contains 

!************************************************************************
      subroutine roo(jamax,jbmax,
     *ihnmn,ihnmx,ihpmn,ihpmx,ipnmn,
     *ipnmx,ippmn,ippmx,isimax)

      use anglib   ! angular momentum staff

      implicit double precision(a-h,o-z)

      include 'types_phon_dens.inc'
!      include 'input_phon_dens.inc'
      include 'formats_phon_dens.inc'


      type (amp2_typ), dimension(:,:), allocatable :: cc!,cc
      type (ro_typ), dimension(:), allocatable :: roc
      integer, dimension (:), allocatable :: ndx,ndc

      type (phon_typ), dimension (:,:), allocatable :: phonbs
      integer, dimension (:), allocatable :: nphon,iphous
      double precision, dimension (:,:,:), allocatable :: xx
      double precision, dimension (:,:,:,:,:,:), allocatable :: csixj
      double precision, dimension (:), allocatable :: rog

      character(len=30) fnamex,fnamec
      character(len=30) fnamer1,fnamer2,fnamer,fnamern,fnameo


!      open(22,file='roo.dat',status='unknown',form='unformatted')



      open(23,file='2phonon/2f_rnh.dat',status='unknown'
     *,form='unformatted') 
      open(32,file='2phonon/2f_rnp.dat',status='unknown'
     *,form='unformatted')
      open(34,file='2phonon/2f_rph.dat',status='unknown'
     *,form='unformatted')
      open(43,file='2phonon/2f_rpp.dat',status='unknown'
     *,form='unformatted')



      jjmx=14

      allocate (csixj(0:jjmx,0:jjmx,0:jjmx,0:jjmx,0:jjmx,0:jjmx))

      csixj=0.d0

      do i1=0,jjmx
       write(*,*)'i1 =',i1
       do i2=0,jjmx
        do i3=0,jjmx
         do i4=0,jjmx
          do i5=0,jjmx
           do i6=0,jjmx             
             csixj(i1,i2,i3,i4,i5,i6)=
     *sixj(2*i1,2*i2,2*i3,2*i4,2*i5,2*i6)
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo

      ndimroc=1000000000
      allocate(roc(ndimroc))

      roc%rho=0.d0
      roc%ib=0
      roc%is=0
      roc%i1=0
      roc%i2=0

      call loadphon(phonbs,nphon)

      write(*,*)' Number of 1phonon states ',nphon(1)
      write(*,*)' Number of 2phonon states ',nphon(2)


      nff=2

      call selphon(nff,phonbs,nphon,iphous)

!     stop

      ndamp=10000
      ndrho=1000000     


      if (nff.eq.2) then 
       fnamex='2phonon/2f_x.dat'
       fnamec='2phonon/2f_c.dat'
      endif

      allocate(rog(nphon(1)))
      rog=0.d0

!     call readc(fnamec,cc,ndc)
 
      call readc2(fnamex,xx,nphon(nff),nphon(1),nphon(1))

      do ipa=-1,1,2
       do ja=0,jamax

        write(*,*)'Ipar, J ',ipa,ja

        call readx(fnamec,cc,ndc,ipa,ja,ilax,nb,nbb)

        do iaaa=ilax+1,ilax+nb


       enerf=phonbs(nff,iaaa)%enf


       iroc=1

        if (iphous(iaaa).eq.1) then 

        write(*,*)'iaaa',iaaa

       write(*,*)'Phonon energy ',phonbs(nff,iaaa)%enf


        ibb=ilac

         do ibbb=1,nphon(nff)

!         write(*,*)'ibbb ',ibbb

          jb=phonbs(nff,ibbb)%jj

          do isi=0,isimax

           do il=1,nphon(1) ! lambda'

             rog=0.d0
             jilap=phonbs(1,il)%jj

             do i=1,nbb


            il1=cc(iaaa,i)%is 
            ila=cc(iaaa,i)%ig
            jila=phonbs(1,ila)%jj
            jil1=phonbs(1,il1)%jj 

            ifaz=(-1)**(isi+jilap+ja+jil1)
            fact=dfloat(ifaz)*(dfloat(2*ja+1))**0.5d0
             
   
       rog(ila)=rog(ila)+fact*
     *cc(iaaa,i)%am*xx(ibbb,il1,il)*csixj(isi,jilap,jila,jil1,ja,jb)

            ifaz=(-1)**(isi+jila+jb+jil1)
            fact=dfloat(ifaz)*(dfloat(2*ja+1))**0.5d0
             
       rog(il1)=rog(il1)+fact*
     *cc(iaaa,i)%am*xx(ibbb,il,ila)*csixj(isi,jilap,jil1,jila,ja,jb)
      
            enddo ! over i

            do ila=1,nphon(1)

            if (dabs(rog(ila)).gt.1.d-8) then 
             roc(iroc)%rho=rog(ila)       
             roc(iroc)%ib=ibbb
             roc(iroc)%is=isi
             roc(iroc)%i1=ila
             roc(iroc)%i2=il

!c             write(998,17)iaaa,ibbb,isi,ila,il,roc(iroc)%rho

            iroc=iroc+1

            if (iroc.gt.ndimroc) then 
            write(*,*)' Increase dimension of array roc '
            stop
            endif
            endif

            enddo

           enddo  !over il

           enddo   ! over isi

         enddo   !over iaaa

         endif
          

        do ityp=-1,1,2 
          do iph=-1,1,2

       if (ityp.eq.-1.and.iph.eq.1) then
!        write(*,*)' Calculating proton particle densities'
        fnamer='1phonon/1f_rpp.dat'
        fnamern='1phonon/1f_rpp.dat'
        ifile=43
        imin=ippmn
        imax=ippmx
!        fnameo='2phonon/2f_rpp.dat'
       endif

       if (ityp.eq.-1.and.iph.eq.-1) then
!        write(*,*)' Calculating proton hole densities'
        fnamer='1phonon/1f_rph.dat'
        fnamern='1phonon/1f_rph.dat'
        ifile=34
        imin=ihpmn
        imax=ihpmx

!        fnameo='2phonon/2f_rph.dat'
       endif


       if (ityp.eq.1.and.iph.eq.1) then
!        write(*,*)' Calculating neutron particle densities'
        fnamer='1phonon/1f_rnp.dat'
        fnamern='1phonon/1f_rnp.dat'
!        fnameo='2phonon/2f_rnp.dat'
        imin=ipnmn
        imax=ipnmx

        ifile=32
       endif

       if (ityp.eq.1.and.iph.eq.-1) then
!        write(*,*)' Calculating neutron hole densities'
        fnamer='1phonon/1f_rnh.dat'
        fnamern='1phonon/1f_rnh.dat'
        ifile=23
        imin=ihnmn
        imax=ihnmx

!        fnameo='2phonon/2f_rnh.dat'
       endif
       
        call ropar2(imin,imax,
     *ityp,iph,jamax,jbmax,isimax,roc,iroc-1,ifile
     *,phonbs,nphon,fnamern)


        enddo
        enddo



!         write(*,*)'iroc', iroc-1
!         write(22)iaaa,iroc-1
!         write(22)(roc(iii)%ib,roc(iii)%is,roc(iii)%i1,roc(iii)%i2,
!     *roc(iii)%rho,iii=1,iroc-1)

!          if ((iroc-1).gt.0) then 
!           roc%rho=0.d0
!           roc%ib=0
!           roc%is=0
!           roc%i1=0
!           roc%i2=0
!          endif

         continue

         enddo

       enddo

      enddo

      return
      end subroutine roo

!************************************************************************


      subroutine test_amp1(jamax,jbmax,isimax)

      use anglib   ! angular momentum staff

      implicit double precision(a-h,o-z)

      include 'types_phon_dens.inc'
      include 'input_phon_dens.inc'
      include 'formats_phon_dens.inc'


      type (amp2_typ), dimension(:,:), allocatable :: cc!,cc
      type (ro_typ), dimension(:), allocatable :: roc
      integer, dimension (:), allocatable :: ndx,ndc

      type (phon_typ), dimension (:,:), allocatable :: phonbs
      integer, dimension (:), allocatable :: nphon
      double precision, dimension (:,:,:), allocatable :: xx,cx
      double precision, dimension (:,:,:,:,:,:), allocatable :: csixj
      double precision, dimension (:), allocatable :: rog



      character(len=30) fnamex,fnamec

      open(22,file='test_amp.dat',status='unknown',form='unformatted')

      jjmx=10

      allocate (csixj(0:jjmx,0:jjmx,0:jjmx,0:jjmx,0:jjmx,0:jjmx))

      csixj=0.d0

      do i1=0,jjmx
       write(*,*)'i1 =',i1
       do i2=0,jjmx
        do i3=0,jjmx
         do i4=0,jjmx
          do i5=0,jjmx
           do i6=0,jjmx             
             csixj(i1,i2,i3,i4,i5,i6)=
     *sixj(2*i1,2*i2,2*i3,2*i4,2*i5,2*i6)
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo


      ndimroc=10000000
      allocate(roc(ndimroc))

      roc%rho=0.d0
      roc%ib=0
      roc%is=0
      roc%i1=0
      roc%i2=0

      call loadphon(phonbs,nphon)

      ndamp=10000
      ndrho=1000000     
      nff=2

      if (nff.eq.2) then 
       fnamex='2phonon/2f_x.dat'
       fnamec='2phonon/2f_c.dat'
      endif

      write(*,*)' Number of 1phonon states ',nphon(1)
      write(*,*)' Number of 2phonon states ',nphon(2)



c      call readc(fnamec,cc,ndc)
 
      call readc2(fnamex,xx,nphon(nff),nphon(1),nphon(1))
      call readc2(fnamec,cx,nphon(nff),nphon(1),nphon(1))
   
      do iaaa=1181,nphon(nff)

       ja=phonbs(nff,iaaa)%jj

        write(*,*)'iaaa',iaaa

         do ibbb=1181,nphon(nff)

          jb=phonbs(nff,ibbb)%jj

           do il=1,nphon(1) ! lambda
           jil=phonbs(1,il)%jj

            do ilp=1,nphon(1) ! lambda'
            jilp=phonbs(1,ilp)%jj


            do isi=0,isimax
            

            rsu=0.d0
            do il1=1,nphon(1)

            jil1=phonbs(1,il1)%jj

             rsixj=csixj(isi,jilp,jil,jil1,ja,jb)

 
             facta=(dfloat(2*ja+1))**0.5d0
             faz1=dfloat((-1)**(isi+jil1+jil+jb))
             faz2=dfloat((-1)**(isi+jil1+jilp+ja))
         
             rsu=rsu+facta*rsixj*(
     *faz1*cx(iaaa,il,il1)*xx(ibbb,ilp,il1)+
     *faz2*cx(iaaa,il1,il)*xx(ibbb,il1,ilp))

             if (dabs(rsu).gt.1.d-6) then 
      write(994,20)iaaa,ibbb,isi,il,ilp,il1,
     *cx(iaaa,il,il1),xx(ibbb,ilp,il1)
      write(994,20)iaaa,ibbb,isi,il,ilp,il1,
     *cx(iaaa,il1,il),xx(ibbb,il1,ilp)        
             endif




                          

            enddo

             if (dabs(rsu).gt.1.d-6)  
     *write(998,17)iaaa,ibbb,isi,il,ilp,rsu


            enddo

             

            enddo

           enddo  !over il

           enddo   ! over isi



         enddo   !over iaaa


         continue





      return
      end subroutine test_amp1


************************************************************************

      subroutine test_amp2(jamax,jbmax,isimax)

      use anglib   ! angular momentum staff

      implicit double precision(a-h,o-z)

      include 'types_phon_dens.inc'
      include 'input_phon_dens.inc'
      include 'formats_phon_dens.inc'


      type (amp2_typ), dimension(:,:), allocatable :: cc!,cc
      type (ro_typ), dimension(:), allocatable :: roc
      integer, dimension (:), allocatable :: ndx,ndc

      type (phon_typ), dimension (:,:), allocatable :: phonbs
      integer, dimension (:), allocatable :: nphon
      double precision, dimension (:,:,:), allocatable :: xx
      double precision, dimension (:,:,:,:,:,:), allocatable :: csixj
      double precision, dimension (:), allocatable :: rog



      character(len=30) fnamex,fnamec

      open(22,file='test_amp.dat',status='unknown',form='unformatted')

      jjmx=10


      ndimroc=10000000
      allocate(roc(ndimroc))

      roc%rho=0.d0
      roc%ib=0
      roc%is=0
      roc%i1=0
      roc%i2=0

      call loadphon(phonbs,nphon)

      ndamp=10000
      ndrho=1000000     
      nff=2

      if (nff.eq.2) then 
       fnamex='2phonon/2f_c.dat'
       fnamec='2phonon/2f_x.dat'
      endif

!      write(*,*)' Number of 1phonon states ',nphon(1)
!      write(*,*)' Number of 2phonon states ',nphon(2)

      allocate(rog(nphon(1)))
      rog=0.d0


c      call readc(fnamec,cc,ndc)
 
      call readc2(fnamex,xx,nphon(nff),nphon(1),nphon(1))

      do ipa=1,-1,-2
       do ja=0,jamax

        write(*,*)'Ipar, J ',ipa,ja

        call readx(fnamec,cc,ndc,ipa,ja,ilax,nb,nbb)

        do iaaa=ilax+1,ilax+nb
c        write(998,*)' ********* ',iaaa-1,iroc

        iroc=1

        write(*,*)'iaaa',iaaa

        ibb=ilac

         do ibbb=1,nphon(nff)

          jb=phonbs(nff,ibbb)%jj

           do ilp=1,nphon(1) ! lambda

c            do ilp=1,nphon(1) ! lambda'

            do i=1,nbb


             rog=0.d0
c             jilap=phonbs(1,il)%jj



            ila=cc(iaaa,i)%is 
            iga=cc(iaaa,i)%ig
c            jila=phonbs(1,ila)%jj
c            jil1=phonbs(1,il1)%jj 

            
           facta=(dfloat(2*ja+1))**0.5d0*cc(iaaa,i)%am*xx(ibbb,ilp,iga)
           if (dabs(facta).gt.1.d-10)  
     *write(997,18)iaaa,ibbb,ila,ilp,facta

c            factb=(dfloat(2*jb+1))**0.5d0
             
  

            

c             enddo


            enddo

           enddo  !over il

           enddo   ! over isi


         enddo   !over iaaa


         continue


         enddo


      enddo

      return
      end subroutine test_amp2


!************************************************************************

      subroutine ropar2(imin,imax,ityp,iph,jamax,jbmax,isimax,rom,irom,
     *ifile,phonbs,nphon,fnamern)

      use anglib   ! angular momentum staff

      implicit double precision(a-h,o-z)

      include 'types_phon_dens.inc'
      include 'input_phon_dens.inc'
      include 'formats_phon_dens.inc'


      type (amp2_typ), dimension(:,:), allocatable :: cc!,cc
      type (ro_typ), dimension(:), allocatable :: roc,rom,rh
      type(rho_typ), dimension(:), allocatable :: rop

      integer, dimension (:), allocatable :: ndx,ndc

      type (phon_typ), dimension (:,:), allocatable :: phonbs
      integer, dimension (:), allocatable :: nphon
      double precision, dimension (:,:,:), allocatable :: xx
      double precision, dimension (:,:,:,:,:,:), allocatable :: csixj
      double precision, dimension (:,:,:,:,:), allocatable :: rho,rhon
      double precision, dimension (:,:), allocatable :: rr



      character(len=30) fnamer1,fnamer2,fnamer,fnamern,fnameo

      iirg=0

      if (irom.ne.0) then  

      jjmx=12
      isimx=isimax



!      call loadphon(phonbs,nphon)

      ndamp=10000
      ndrho=10000000     
      nff=2


!      write(*,*)' Number of 1phonon states ',nphon(1)
!      write(*,*)' Number of 2phonon states ',nphon(2)
!      write(*,*)' Imin, Imax ',imin,imax


!c      open(88,file=fnamer1,status='old',form='unformatted')
!      open(77,file=fnamer2,status='old',form='unformatted')


      




      if (.not.allocated(rr)) allocate(rr(nphon(2),0:15))

      rr=0.d0


!      call readrho(fnamer,rho)
      call readrho(fnamern,rhon)


!c      imin=3
!c      imax=6

!      open(33,file=fnameo,status='unknown'
!     *,form='unformatted')

      ndrh=100000000
      if (.not.allocated(rh)) allocate(rh(ndrh))

 
!      do iaaa=1,nphon(nff)
        
        iirg=0
        

!        call readro(fnamer2,77,iaaa,rom,irom)

!        write(*,*)'***',iaaa,iror,irom

        do i1=imin,imax
         do i2=imin,imax

        rr=0.d0


         do ii=1,irom

         isir=rom(ii)%is
         iba=rom(ii)%ib
         ilam=rom(ii)%i1
         ilamp=rom(ii)%i2
 
       rr(iba,isir)=rr(iba,isir)+rhon(ilam,ilamp,isir,i1,i2)*rom(ii)%rho

             
            
         enddo ! over ii


         do iba=1,nphon(nff)
          do isir=0,isimx

          if (dabs(rr(iba,isir)).gt.1.d-10) then

          iirg=iirg+1
          rh(iirg)%ib=iba
          rh(iirg)%is=isir
          rh(iirg)%i1=i1
          rh(iirg)%i2=i2
          rh(iirg)%rho=rr(iba,isir)
          endif
          enddo
         enddo

          enddo    
         enddo

!      write(*,*)'iirg = ',iirg   
 666  continue
      write(ifile)iaaa,iirg
      if (iirg.gt.0) then
      write(ifile)(rh(iii)%ib,iii=1,iirg)
      write(ifile)(rh(iii)%is,iii=1,iirg)
      write(ifile)(rh(iii)%i1,iii=1,iirg)
      write(ifile)(rh(iii)%i2,iii=1,iirg)
      write(ifile)(rh(iii)%rho,iii=1,iirg)
      endif



      endif 

      if (iirg.eq.0) then
      write(ifile)iirg
      write(ifile)iirg
      write(ifile)iirg
      write(ifile)iirg
      write(ifile)dfloat(iirg)
      endif



!      enddo

!      close(33)
!      close(88)
!      close(77)

!      deallocate(rh,rr,rom)

      return
      end subroutine ropar2
!************************************************************************

      subroutine loadphon(phonbs,nphon)

      implicit double precision (a-h,o-z)

      include 'types_phon_dens.inc'

      character*30 namef

      type (phon_typ), dimension (:,:), allocatable :: phonbs
      integer, dimension (:), allocatable :: nphon


      allocate (nphon(0:4))
      nphon=0
      allocate (phonbs(0:4,1:1000000))

      do nfon=1,2

      if (nfon.eq.1) namef='1phonon/1f_states.dat'
      if (nfon.eq.2) namef='2phonon/2f_states.dat'
      
      open(1,file=namef,status='old',form='unformatted')

      do while (.not.eof(1))
       read(1)i,ipartt,ijj,en
         phonbs(nfon,i)%par=ipartt
         phonbs(nfon,i)%jj=ijj
         phonbs(nfon,i)%enf=en
      enddo
      nphon(nfon)=i

      close(1)

      enddo

      end subroutine loadphon

*********************************************************************
      subroutine selphon(nf,phonbs,nphon,iphous)

      implicit double precision (a-h,o-z)

      include 'types_phon_dens.inc'

      character*30 namef

      type (phon_typ), dimension (:,:), allocatable :: phonbs
      integer, dimension (:), allocatable :: nphon,iphous


      allocate(iphous(nphon(nf)))
      iphous=0


      write(*,*)'energy threshold for 2 phonon states'
      read(*,*)xthre
      write(*,*)xthre

!      xthre=30.0d0


      ii=0
      do i=1,nphon(nf)
       xene=phonbs(nf,i)%enf
       if (xene.le.xthre) then
        ii=ii+1
        iphous(i)=1
       endif
      enddo

      write(*,*)' Number of selected phonons ',ii

      end subroutine selphon

********************************************************************** 
    

      subroutine readx(fname,xcc,ndcc,ipcal,jcal,ilamps,no,idphon)

      implicit double precision (a-h,o-z)

      include 'types_phon_dens.inc'

      type (amp2_typ), dimension(:,:), allocatable :: xcc
      integer, dimension (:), allocatable :: ndcc

      character(len=30)fname

c      allocate(xcc(ndimi,ndimj))
c      allocate(ndcc(ndimi))

      open(2,file=fname,status='old',form='unformatted')

      ilamp=0

      ndimj=100000

     

      do while (.not.eof(2))

      read(2)ipar,ijj,no,idphon

      if (allocated(xcc).eq..TRUE.) deallocate(xcc,ndcc)

      allocate(xcc(ilamp+1:ilamp+no+1,ndimj))
      allocate(ndcc(ilamp+1:ilamp+1+no))

      xcc%is=0
      xcc%ig=0
      xcc%am=0.d0
      ndcc=0


      do ilam=1,no

       
         if (idphon.gt.ndimj) then 
      write(*,*)'Readcam: allocate bigger array in ndimj',idphon,ndimj
               stop
           endif
  
       read(2)(xcc(ilam+ilamp,i)%ig,xcc(ilam+ilamp,i)%is
     *,xcc(ilam+ilamp,i)%am,i=1,idphon)
       ndcc(ilam+ilamp)=idphon
      enddo

      ilamps=ilamp
      ilamp=ilamp+no


      if (ipar.eq.ipcal.and.ijj.eq.jcal) goto 11


      deallocate(xcc,ndcc)
      no=0
      idphon=0

      enddo

  11  continue     

      close(2)

      return
      end subroutine readx
************************************************************

      subroutine readro(fname,ifile,ia,roc,iror)

      implicit double precision (a-h,o-z)

      include 'types_phon_dens.inc'

      type (ro_typ), dimension(:), allocatable :: roc

      character(len=30)fname


c      open(22,file=fname,status='old',form='unformatted')

     

c      do while (.not.eof(22))


      read(ifile)iaaa,iror

      if (iaaa.ne.ia) then 
              write(*,*)' Error in reading in readro'
              stop

      endif

      if (allocated(roc)) deallocate(roc)

      allocate(roc(iror))


      read(ifile)(roc(iii)%ib,roc(iii)%is,roc(iii)%i1,roc(iii)%i2,
     *roc(iii)%rho,iii=1,iror)


c      if (iaaa.eq.ia) goto 11


c      enddo

c   11 continue   

c      close(22)

      return
      end subroutine readro
****************************************

      subroutine readc2(fname,cc,n1,n2,n3)

      implicit double precision (a-h,o-z)

      include 'types_phon_dens.inc'

      type (amp2_typ), dimension(:,:), allocatable :: xcc
      double precision, dimension (:,:,:), allocatable :: cc
      integer, dimension (:), allocatable :: ndcc

      character(len=30)fname

      ndimi=2000
      ndimj=100

      allocate (cc(n1,n2,n3))
      cc=0.d0


      open(2,file=fname,status='old',form='unformatted')

      ilamp=0

      ndimj=100000

     

      do while (.not.eof(2))

      read(2)ipar,ijj,no,idphon

      if (allocated(xcc).eq..TRUE.) deallocate(xcc,ndcc)

      allocate(xcc(1,ndimj))
      allocate(ndcc(1))


      do ilam=1,no

       
         if (idphon.gt.ndimj) then 
      write(*,*)'Readcam: allocate bigger array in ndimj',idphon,ndimj
               stop
           endif
  
       read(2)(xcc(1,i)%ig,xcc(1,i)%is
     *,xcc(1,i)%am,i=1,idphon)
       ndcc(1)=idphon

       do i=1,idphon
        ii=xcc(1,i)%is
        jj=xcc(1,i)%ig
        cc(ilam+ilamp,ii,jj)=xcc(1,i)%am
       enddo
      enddo


   

      ilamps=ilamp
      ilamp=ilamp+no


c      if (ipar.eq.ipcal.and.ijj.eq.jcal) goto 11


c      deallocate(xcc,ndcc)
c      no=0
c      idphon=0

      enddo

  11  continue     



      close(2)

      deallocate(xcc,ndcc)

      return
      end subroutine readc2
************************************************************
      subroutine readr(fname,xcc,ndcc,iaaa)

      implicit double precision (a-h,o-z)

      include 'types_phon_dens.inc'

      type (amp2_typ), dimension(:,:), allocatable :: xcc
      integer, dimension (:), allocatable :: ndcc

      character(len=30)fname


      open(2,file=fname,status='old',form='unformatted')

      ilamp=0

      ndimj=100000

     

      do while (.not.eof(2))

      read(2)ipar,ijj,no,idphon

      if (allocated(xcc).eq..TRUE.) deallocate(xcc,ndcc)

      allocate(xcc(ilamp:ilamp+no,ndimj))
      allocate(ndcc(ilamp:ilamp+no))


      do ilam=1,no

       
         if (idphon.gt.ndimj) then 
      write(*,*)'Readcam: allocate bigger array in ndimj',idphon,ndimj
               stop
           endif
  
       read(2)(xcc(ilam+ilamp,i)%ig,xcc(ilam+ilamp,i)%is
     *,xcc(ilam+ilamp,i)%am,i=1,idphon)
       ndcc(ilam+ilamp)=idphon
      enddo

      ilamps=ilamp
      ilamp=ilamp+no


      if (ipar.eq.ipcal.and.ijj.eq.jcal) goto 11


      deallocate(xcc,ndcc)
      no=0
      idphon=0

      enddo

  11  continue     

      close(2)

      return
      end subroutine readr
************************************************************
      subroutine readrho(fname,rho)

      implicit double precision (a-h,o-z)

c      include 'formats_eqm.inc'
      include 'types_phon_dens.inc'

      double precision, dimension (:,:,:,:,:), allocatable :: rho

      type(rho_typ), dimension(:), allocatable :: ron

      character(len=30)fname

      allocate (rho(104,104,0:10,10,10))

      rho=0.d0

      ndro=10000
      allocate(ron(ndro))


      open(2,file=fname,status='old',form='unformatted')

      do while (.not.eof(2))

       read(2)igg,ndgg

       if (ndgg.gt.ndro) then
                write(*,*)'WARNING: Increase dimension in readrho'
                stop
       endif

       read(2)(ron(ii)%ilap,ii=1,ndgg)
       read(2)(ron(ii)%j,ii=1,ndgg)
       read(2)(ron(ii)%i1,ii=1,ndgg)
       read(2)(ron(ii)%i2,ii=1,ndgg)
       read(2)(ron(ii)%ro,ii=1,ndgg)

       do i=1,ndgg
        ilapp=ron(i)%ilap
        jp=ron(i)%j
        i1p=ron(i)%i1
        i2p=ron(i)%i2
        rho(igg,ilapp,jp,i1p,i2p)=ron(i)%ro
       enddo

       enddo


      close(2)

      deallocate(ron)

      return
      end subroutine readrho

***********************************************************************








      end module rdens
      
