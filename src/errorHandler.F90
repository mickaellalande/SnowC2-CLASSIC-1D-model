!> \file
!> Prints the name of the subroutine and an error code when
!! an error condition is encountered.
!! @author E. Chan
!
subroutine errorHandler (name, n) ! based on xit

  !     * Terminates a program by printing the program name and
  !     * a line across the page followed by a number n.
  !
  !     * n>=0 is for a normal end. the line is dashed.
  !     * n<0 is for an abnormal end. the line is dotted.
  !   
  !     * For CLASSIC offline, the routine "abandonCell" is called.
  !     * In coupled mode, the termination is handled via the routine "xit".

  use generalUtils, only : abandonCell

  implicit none

  character * ( * ), intent(in) :: name    !< Name of the subroutine in which the error was found
  integer, intent(in) :: n                 !< N: error code
  integer             :: i
  integer             :: length
  character*80        :: dash !<
  character*80        :: star !<
  !
  data dash(1:40)/'----------------------------------------'/
  data dash(41:80)/'----------------------------------------'/
  data star(1:40)/'****************************************'/
  data star(41:80)/'****************************************'/
  !---------------------------------------------------------------------
  !
  !>
  !! In CLASSIC, this subroutine is called when a test of ambient values
  !! of selected variables is performed and an abnormal condition is
  !! encountered. The name of the subroutine in which the condition
  !! arose is passed in, and is printed together with an error code,
  !! flagging the location of the error in the subroutine. 
  !! How termination is handled depends on whether the model is running
  !! in offline or coupled mode.
  !!
#if defined without_agcm_
  length = len_trim(name)
  if (n >= 0) write(6,6010) dash(1:8), name(1:length), dash(1:81-length),n
  if (n < 0)  write(6,6010) star(1:8), name(1:length), star(1:81-length),n
  call abandonCell
#else
  ! FLAG: Temporary if-block to bypass the call to xit in coupled mode for specific cases. 
  !       Remove if-block when resolved and just keep the call to xit.
  if ((name == 'balcar' .and. n == -14) .or. (name == 'luc' .and. n == -16)) then
    length = len_trim(name)
    write(6,6010) dash(1:8), name(1:length), dash(1:81-length),n
  else
    call xit(name,n)
  end if
#endif
  !
  !---------------------------------------------------------------------
  6010 format('0',a,'  END  ',a,1x,a,i8)
end subroutine errorHandler
