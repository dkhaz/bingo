program gather

  integer::i,j,total_chains,pt_per_chain,r3_points,total_lines,countline
  
  character(len=128)::filename,numstr
  double precision, allocatable,dimension(:)::modek2,modek3,modek1
  double precision, allocatable,dimension(:)::f_nl_re
  double precision ::k1,k2,k3 
  
  total_chains=60
  total_lines=3600
  allocate(modek2(total_chains*total_lines),modek3(total_chains*total_lines),modek1(total_chains*total_lines),f_nl_re(total_chains*total_lines))


  
  countline=0
  do i=1,total_chains
   
    
     write(numstr,'(i10)') i
     
     filename=trim("F_nl")//"_"//trim(adjustl(numstr))//".txt"
     
       
    
   open(unit=900+i,file="plots/"//filename,action='read',status='old') 
     
     do k=1,total_lines
     countline=countline+1  
        read(900+i,*)modek1(countline),modek2(countline),modek3(countline),f_nl_re(countline)  
     end do
 
     
     close(900+i)

  end do
  
  open(unit=901,file='plots/F_nl_3d.txt',form='formatted',status='replace') 

  do i=1,total_chains*total_lines
       write(901,'(999E15.5)')log10(modek2(i)),log10(modek1(i)),log10(modek3(i)),real(f_nl_re(i))
  end do 

  close(901)


end program gather
