program main


 use parameters, only : x0,x,count_restart
 implicit none
 
 x0(:) = 2
 call optimisation(x0, 100, 1e-8)    
 print *, x
    
 print *, "Restart : ", count_restart
 end program  main

