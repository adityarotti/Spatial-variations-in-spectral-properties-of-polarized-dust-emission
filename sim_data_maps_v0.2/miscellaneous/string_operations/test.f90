implicit none
character*10 :: s

s="abcdefghij"
print*,s(:3)
print*,s(4:7)
print*,s(7:10)

stop
end
