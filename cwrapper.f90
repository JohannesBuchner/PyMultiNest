! iso_c_binding interface for MultiNest_v3.2 Nested::nestRun
! Michele Vallisneri, 2013/09/20
!
! allows calling multinest with the "natural" C prototype
!
! void run(bool nest_IS,bool nest_mmodal,bool nest_ceff, \
!          int nest_nlive,double nest_tol,double nest_ef,int nest_ndims,int nest_totPar,int nest_nCdims,int maxClst, \
!          int nest_updInt,double nest_Ztol,char nest_root[],int seed,int nest_pWrap[], \
!          bool nest_fb,bool nest_resume,bool nest_outfile,bool initMPI,double nest_logZero,int nest_maxIter, \
!          double (*loglike)(double *,int,int,void *context), \
!          void (*dumper)(int,int,int,double *,double *,double *,double,double,double,void *context),void *context);
!
! as well as
!
! double loglike(double *Cube,int n_dim,int nPar,void *context);
! void dumper(int nSamples,int nlive,int nPar, \
!             double *physLive,double *posterior,double *paramConstr, \
!             double maxLogLike,double logZ,double logZerr,void *context);
!
! note that we are assuming that (void *) is the same size as int, but that's what multinest uses

module cnested

  contains

  subroutine run(nest_IS,nest_mmodal,nest_ceff,nest_nlive,nest_tol,nest_ef,nest_ndims,nest_totPar,nest_nCdims,maxClst, &
    nest_updInt,nest_Ztol,nest_root,seed,nest_pWrap,nest_fb,nest_resume,nest_outfile,initMPI,nest_logZero,nest_maxIter, &
    loglike,dumper,context) bind(c)

    use iso_c_binding, only: c_int, c_bool, c_double, c_char, c_funptr, c_ptr, C_NULL_CHAR
    use Nested, only: nestRun
    implicit none

    integer(c_int),  intent(in), value :: nest_ndims,nest_nlive,nest_updInt,seed
    integer(c_int),  intent(in), value :: maxClst,nest_totPar,nest_nCdims,nest_maxIter
    integer(c_int),  intent(in) :: nest_pWrap(nest_ndims)
    logical(c_bool), intent(in), value :: nest_IS,nest_mmodal,nest_fb,nest_resume,nest_ceff,nest_outfile,initMPI
    character(kind=c_char,len=1), dimension(1), intent(in) :: nest_root
    real(c_double),  intent(in), value :: nest_tol,nest_ef,nest_Ztol,nest_logZero
    type(c_funptr),  intent(in), value :: loglike, dumper
    type(c_ptr),     intent(in) :: context

    character(len=100) :: fnest_root
    integer :: i, context_f

    fnest_root = ' '
    do i = 1, 100
      if (nest_root(i) == C_NULL_CHAR) then
        exit
      else
        fnest_root(i:i) = nest_root(i)
      end if
    end do

    ! context_f = transfer(context,context_f)

    call nestRun(logical(nest_IS),logical(nest_mmodal),logical(nest_ceff), &
      nest_nlive,nest_tol,nest_ef,nest_ndims,nest_totPar,nest_nCdims,maxClst, &
      nest_updInt,nest_Ztol,fnest_root,seed,nest_pWrap, &
      logical(nest_fb),logical(nest_resume),logical(nest_outfile),logical(initMPI),nest_logZero,nest_maxIter, &
      loglike_f,dumper_f,context_f)

    contains

    subroutine loglike_f(Cube,n_dim,nPar,lnew,context_pass)
      use iso_c_binding, only: c_double, c_f_procpointer
      implicit none

      integer          :: n_dim,nPar,context_pass
      double precision :: Cube(nPar)
      double precision :: lnew

      interface
        real(c_double) function loglike_proto(Cube,n_dim,nPar,context)
          use iso_c_binding, only: c_int, c_double, c_ptr
          implicit none

          integer(c_int), intent(in), value :: n_dim,nPar
          real(c_double), intent(inout)     :: Cube(nPar)
          integer(c_int), intent(in) :: context
          ! better, but "transfer" is problematic:
          ! type(c_ptr),  intent(in) :: context
        end function loglike_proto
      end interface

      procedure(loglike_proto), pointer :: loglike_c
      call c_f_procpointer(loglike,loglike_c)

      ! type(c_ptr) :: context_c
      ! context_c = transfer(context_pass,context_c)

      lnew = loglike_c(Cube,n_dim,nPar,context_pass)

    end subroutine loglike_f

    subroutine dumper_f(nSamples,nlive,nPar,physLive,posterior,paramConstr,maxLogLike,logZ,logZerr,context_pass)
      use iso_c_binding, only: c_double, c_f_procpointer
      implicit none

      integer          :: nSamples, nlive, nPar, context_pass
      double precision, pointer :: physLive(:,:), posterior(:,:), paramConstr(:)
      double precision :: maxLogLike, logZ, logZerr

      interface
        subroutine dumper_proto(nSamples,nlive,nPar,physLive,posterior,paramConstr,maxLogLike,logZ,logZerr,context)
          use iso_c_binding, only: c_int, c_double, c_ptr
          implicit none

          integer(c_int), intent(in), value :: nSamples, nlive, nPar
          real(c_double), intent(in) :: physLive(nlive,nPar+1), posterior(nSamples,nPar+2),paramConstr(4*nPar)
          real(c_double), intent(in), value :: maxLogLike, logZ, logZerr
          integer(c_int), intent(in), value :: context
          ! better, but "transfer" is problematic:
          ! type(c_ptr),  intent(in) :: context
        end subroutine dumper_proto
      end interface

      procedure(dumper_proto), pointer :: dumper_c
      call c_f_procpointer(dumper,dumper_c)

      ! type(c_ptr) :: context_c
      ! context_c = transfer(context_pass,context_c)

      call dumper_c(nSamples,nlive,nPar,physLive,posterior,paramConstr,maxLogLike,logZ,logZerr,context_pass)

    end subroutine dumper_f

  end subroutine

end module
