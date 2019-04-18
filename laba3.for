      program main
        implicit none
        integer mems(268435456), pt, N, jN, BFSZ ! 1 Gb
        common /memory/ mems, pt
        pt = 0
        call fillIArray("sizes.txt", mems(1), 9, 2, .FALSE.)
        N = mems(1)
        jN = mems(2)
        print *, N, jN
        pt = 2*jN + 3*N + 2
        BFSZ = N + 1
        call process(mems(2*jN + N + 1), mems(2*jN + 2*N + 1 + 1),
     *  mems(1), mems(jN + 1), mems(2*jN + 1),
     *  mems(pt), mems(pt + bfsz), mems(pt + bfsz * 2)
     *  , N, jN, BFSZ)
      end program main

      subroutine multRk(id, di, ja, el, vt, rs, N, jN)
        implicit none
        integer N, jN, i, j, k
        real el(jN), vt(N), rs(N), di(N)
        integer id(N + 1), ja(jN)
        do i = 1,N
            rs(i) = rs(i) + vt(i)*di(i)
        end do
        do i = 1,N
            do j = id(i), id(i + 1) - 1
                k = ja(j)
                rs(i) = rs(i) + vt(k) * el(j)
                rs(k) = rs(k) + vt(i) * el(j)
            end do
        end do
      end subroutine multRk

      subroutine clearArray(ms, N)
        implicit none
        real ms(N)
        integer  N, i
        do i = 1, N
            ms(i) = 0
        end do
      end subroutine clearArray

      subroutine printRMatrixToFile(id, di, ja, el, bf, N, jN, fn, Ln)
        implicit none
        integer N, jN, i, j, k, l, Ln
        real el(jN), di(N), bf(N)
        integer id(N + 1), ja(jN)
        character(Ln) fn
        open(57, file=fn, ERR=67, ACTION='WRITE')
        do i = 1, N
            call clearArray(bf, N)
            do j = id(i), id(i + 1) - 1
                k = ja(j)
                bf(k) = el(j)
            end do
            bf(i) = di(i)
            do l = i + 1, N
                do j = id(l), id(l + 1) - 1
                    k = ja(j)
                    if (k.eq.i) then
                        bf(l) = el(j)
                    end if
                end do
            end do
            write (57,*) bf
        end do
        close(57)
        goto 70
67      STOP "Error: Unknown error!"
70    end subroutine printRMatrixToFile

      subroutine printIMatrixToFile(id, di, ja, el, bf, N, jN, fn, Ln)
        implicit none
        integer N, jN, i, j, k, l, Ln
        integer id(N + 1), di(N), ja(jN), bf(N), el(jN)
        character(Ln) fn
        open(57, file=fn, ERR=67, ACTION='WRITE')
        do i = 1, N
            call clearArray(bf, N)
            do j = id(i), id(i + 1) - 1
                k = ja(j)
                bf(k) = el(j)
            end do
            bf(i) = di(i)
            do l = i + 1, N
                do j = id(l), id(l + 1) - 1
                    k = ja(j)
                    if (k.eq.i) then
                        bf(l) = el(j)
                    end if
                end do
            end do
            write (57,*) bf
        end do
        close(57)
        goto 70
67      STOP "Error: Unknown error!"
70    end subroutine printIMatrixToFile

      subroutine process(id, di, ja, el, vt, bf1, bf2, bf3, N, jN, BFSZ)
          implicit none
          integer N, jN, BFSZ
          real el(jN), vt(N), bf1(N), di(N)
          integer id(N + 1), ja(jN)
          integer  bf2(bfsz), bf3(N)
          call fillIArray("ja.txt", ja, 6, jN, .TRUE.)
          call fillRArray("elems.txt", el, 9, jN, .TRUE.)
          call fillRArray("vector.txt", vt, 10, N, .TRUE.)
          call fillIArray("indexes.txt", id, 11, N + 1, .TRUE.)
          call fillRArray("diag.txt", di, 8, N, .TRUE.)
          call multRk(id, di, ja, el, vt, bf1, N, jN)
          call printRarray(bf1, N, "res.txt", 7, .FALSE.)
          print *,  bf1
          print *,  vt
          call printRMatrixtoFile(id, di, ja, el,bf2, N, jN,"out.txt",7)
      end subroutine process

      subroutine readSizes(N, jN)
        implicit none
        integer N, jN
        open(56, file= "sizes.txt", ACTION='READ')
        read (56, *, END=66, ERR=66) N, jN
        goto 70
66      STOP "Error: the file must contain only 2 numbers! "
70      close(56)
      end subroutine readSizes

      subroutine fillIArray(fn, arr, Ln, N, b)
          implicit none
          logical b
          integer N, i, Ln
          integer arr(N)
          character(Ln) fn
          if (b) then
            open(57,file=fn,ERR=67,ACTION='READ'
     *       ,ACCESS = 'DIRECT',RECL=4)
            do i=1,N
                read (57, ERR=66, rec = i) arr(i)
            end do
          else
            open(57,file=fn,ERR=67,ACTION='READ')
            read (57, *) (arr(i), i=1,N)
          end if
          close(57)
          goto 70
66        write (*,'(A, I4, A)') "Error: file must contain only",
     *    N, " numbers!"
          STOP
67        write (*,'(A, A, A)') "Error: File ", fn, " not found!"
          STOP
70    end subroutine fillIArray

      subroutine fillRArray(fn, arr, Ln, N, b)
          implicit none
          logical b
          integer N, i, Ln
          real arr(N)
          character(Ln) fn
          if (b) then
            open(57,file=fn,ERR=67,ACTION='READ'
     *       ,ACCESS = 'DIRECT',RECL=4)
            do i=1,N
                read (57, ERR=66, rec = i) arr(i)
            end do
          else
            open(57,file=fn,ERR=67,ACTION='READ')
            read (57, *) (arr(i), i=1,N)
          end if
          close(57)
          goto 70
66        write (*,'(A, I4, A)') "Error: file must contain only",
     *    N, " numbers!"
          STOP
67        write (*,'(A, A, A)') "Error: File ", fn, " not found!"
          STOP
70    end subroutine fillRArray


      subroutine printIarray(di, N, fn, fnsz, b)
      implicit none
      logical b
      integer N, fnsz, i
      character (fnsz) fn
      integer di(N)
      if (b) then
        open (57, file=fn,ACTION='WRITE', access='DIRECT',RECL=4)
        do i = 1, N
            write (57, rec=i) di(i)
        end do
      else
        open (57, file=fn, ACTION='WRITE')
        write (57, *) di
      end if
      close(57)
      end subroutine printIarray

      subroutine printRarray(di, N, fn, fnsz, b)
      implicit none
      logical b
      integer N, fnsz, i
      character (fnsz) fn
      real di(N)
       if (b) then
        open (57, file=fn,ACTION='WRITE', access='DIRECT',RECL=4)
        do i = 1, N
            write (57, rec=i) di(i)
        end do
      else
        open (57, file=fn, ACTION='WRITE')
        write (57, *) di
      end if
      close(57)
      end subroutine printRarray
