program main
  use m_examples_pprts_box_cld

  call ex_pprts_box_cld()

  if (myid .eq. 0) then
    print *, ''
    print *, ''
    print *, 'Call this example e.g. with options as in: misc/box_cld_example.sh'
  end if
end program
