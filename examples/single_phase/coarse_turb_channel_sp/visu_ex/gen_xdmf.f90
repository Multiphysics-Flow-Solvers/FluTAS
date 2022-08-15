program gen_xdmf
implicit none
include 'param.h90'
character(len=1), parameter :: lf = char(10)
character(len=400) :: buffer
character(len=6) :: ichar
character(len=3) :: fldname
character(len=4) :: islicechar
character(len=1) :: inormchar
integer :: ixdmf
integer :: e_io,indent
integer :: nflds
integer :: i,ii
!
ixdmf = 99 
nflds = (fldend-fldstart)/nskip + 1
indent = 0
open(unit =ixdmf, file = 'viewfld.xdmf',form='formatted')
write(unit=ixdmf,fmt='(A)') '<?xml version="1.0" ?>'
write(unit=ixdmf,fmt='(A)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
write(unit=ixdmf,fmt='(A)') '<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
write(unit=ixdmf,fmt='(A)') '<Domain>'
indent = indent + 4
  write(buffer,fmt='(A,3I5,A)') repeat(' ',indent)//'<Topology name="TOPO" TopologyType="3DCoRectMesh" Dimensions="',nz,ny,nx,'"/>'
  write(unit=ixdmf,fmt='(A)') trim(buffer)
  write(buffer,fmt='(A)') repeat(' ',indent)//'<Geometry name="GEO" Type="ORIGIN_DXDYDZ">'
  write(unit=ixdmf,fmt='(A)')trim(buffer)
  indent = indent + 4
    write(buffer,fmt='(A)') repeat(' ',indent)//'<DataItem Format="XML" Dimensions="3">'
    write(unit=ixdmf,fmt='(A)')trim(buffer)
    write(buffer,fmt='(A,3E15.6)') repeat(' ',indent),z0,y0,x0
    write(unit=ixdmf,fmt='(A)')trim(buffer)
    write(buffer,fmt='(A)') repeat(' ',indent)//'</DataItem>'
    write(unit=ixdmf,fmt='(A)')trim(buffer)
    write(buffer,fmt='(A)') repeat(' ',indent)//'<DataItem Format="XML" Dimensions="3">'
    write(unit=ixdmf,fmt='(A)')trim(buffer)
    write(buffer,fmt='(A,3E15.6)') repeat(' ',indent),dz,dy,dx
    write(unit=ixdmf,fmt='(A)')trim(buffer)
    write(buffer,fmt='(A)') repeat(' ',indent)//'</DataItem>'
    write(unit=ixdmf,fmt='(A)')trim(buffer)
    indent = indent - 4
  write(buffer,fmt='(A)') repeat(' ',indent)//'</Geometry>'
  write(unit=ixdmf,fmt='(A)')trim(buffer)
  write(buffer,fmt='(A)') repeat(' ',indent)//'<Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">'
  write(unit=ixdmf,fmt='(A)')trim(buffer)
  indent = indent + 4
    write(buffer,fmt='(A)') repeat(' ',indent)//'<Time TimeType="List">'
    write(unit=ixdmf,fmt='(A)')trim(buffer)
    indent = indent + 4
      write(buffer,fmt='(A,I6,A)') repeat(' ',indent)//'<DataItem Format="XML" NumberType="Float" Dimensions="',nflds*nscal,'">'
      write(unit=ixdmf,fmt='(A)')trim(buffer) 
      write(buffer,fmt='(A,I6)') repeat(' ',indent)
      write(ixdmf,fmt='(A)',advance='no') trim(buffer)
      do i = fldstart,fldend,nskip !1,nflds
        write(ixdmf,fmt='(E15.6)',advance='no') t0 + 1.d0*(i-1)*dt!1.*i
      enddo
      write(buffer,fmt='(A)') repeat(' ',indent)//'</DataItem>'
      write(unit=ixdmf,fmt='(A)')trim(buffer)
      indent = indent - 4
    write(buffer,fmt='(A)') repeat(' ',indent)//'</Time>'
    write(unit=ixdmf,fmt='(A)')trim(buffer)
    do i = fldstart,fldend,nskip
      write(ichar,fmt='(i6.6)') i
      write(buffer,fmt='(A)') repeat(' ',indent)//'<Grid Name="T'//ichar//'" GridType="Uniform">'
      write(unit=ixdmf,fmt='(A)')trim(buffer)
      indent = indent + 4
        write(buffer,fmt='(A)') repeat(' ',indent)//'<Topology Reference="/Xdmf/Domain/Topology[1]"/>'
        write(unit=ixdmf,fmt='(A)')trim(buffer)
        write(buffer,fmt='(A)') repeat(' ',indent)//'<Geometry Reference="/Xdmf/Domain/Geometry[1]"/>'
        write(unit=ixdmf,fmt='(A)')trim(buffer)
        do ii = 1,nscal
          write(buffer,fmt='(A)') repeat(' ',indent)//'<Attribute Name="'//scalname(ii)//'" Center="Node">'
          write(unit=ixdmf,fmt='(A)')trim(buffer)
          indent = indent + 4
            write(buffer,fmt='(A,3I5,A)') repeat(' ',indent)//'<DataItem Format="Binary"' // &
                                                                 ' DataType="Float" Precision="8" Endian="Native"' // &
                                                                 ' Dimensions="',nz,ny,nx,'">'
            write(unit=ixdmf,fmt='(A)')trim(buffer)
            indent = indent + 4
              write(buffer,fmt='(A,i9.9,A)') repeat(' ',indent)//scalname(ii)//'_fld_',i,'.bin'
              write(unit=ixdmf,fmt='(A)')trim(buffer)
              indent = indent - 4
            write(buffer,fmt='(A)') repeat(' ',indent)//'</DataItem>'
            write(unit=ixdmf,fmt='(A)')trim(buffer)
            indent = indent - 4
          write(buffer,fmt='(A)') repeat(' ',indent)//'</Attribute>'
          write(unit=ixdmf,fmt='(A)')trim(buffer)
        enddo
        do ii = 1,nscal
          write(buffer,fmt='(A)') repeat(' ',indent)//'<Attribute Name="'//scalname(ii)//'_0" Center="Node">'
          write(unit=ixdmf,fmt='(A)')trim(buffer)
          indent = indent + 4
            write(buffer,fmt='(A,3I5,A)') repeat(' ',indent)//'<DataItem Format="Binary"' // &
                                                                 ' DataType="Float" Precision="8" Endian="Native"' // &
                                                                 ' Dimensions="',nz,ny,nx,'">'
            write(unit=ixdmf,fmt='(A)')trim(buffer)
            indent = indent + 4
              write(buffer,fmt='(A,i9.9,A)') repeat(' ',indent)//scalname(ii)//'_fld_',fldinit,'.bin'
              write(unit=ixdmf,fmt='(A)')trim(buffer)
              indent = indent - 4
            write(buffer,fmt='(A)') repeat(' ',indent)//'</DataItem>'
            write(unit=ixdmf,fmt='(A)')trim(buffer)
            indent = indent - 4
          write(buffer,fmt='(A)') repeat(' ',indent)//'</Attribute>'
          write(unit=ixdmf,fmt='(A)')trim(buffer)
        enddo
        indent = indent - 4
      write(buffer,fmt='(A)') repeat(' ',indent)//'</Grid>'
      write(unit=ixdmf,fmt='(A)')trim(buffer)
    enddo
    indent = indent - 4
  write(buffer,fmt='(A)') repeat(' ',indent)//'</Grid>'
  write(unit=ixdmf,fmt='(A)')trim(buffer)
indent = indent - 4
write(buffer,fmt='(A)') repeat(' ',indent)//'</Domain>'
write(unit=ixdmf,fmt='(A)')trim(buffer)
write(buffer,fmt='(A)') repeat(' ',indent)//'</Xdmf>'
write(unit=ixdmf,fmt='(A)')trim(buffer)
!
close(ixdmf)
end program gen_xdmf
