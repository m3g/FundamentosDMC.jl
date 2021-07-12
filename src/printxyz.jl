"""

```
printxyz(time,x,sys,traj_file_buffer,vel_file_buffer)
```

Subroutine that writes the current coordinates and velocities to 
output files, previously opened in respective buffers.

"""
function printxyz(time,x,v,sys,traj_file_buffer,vel_file_buffer)
  n = length(x)
  println(traj_file_buffer,n)
  println(traj_file_buffer," time = ", time)
  println(vel_file_buffer,n)
  println(vel_file_buffer," time = ", time)
  for i in 1:n
    xt = image(x[i],sys.sides)
    vt = v[i]
    if length(xt) == 2
      print(traj_file_buffer,@sprintf("He %12.5f %12.5f %12.5f\n", xt[1], xt[2], 0.))
      print(vel_file_buffer,@sprintf("He %12.5f %12.5f %12.5f\n", vt[1], vt[2], 0.))
    else
      print(traj_file_buffer,@sprintf("He %12.5f %12.5f %12.5f\n", xt[1], xt[2], xt[3]))
      print(vel_file_buffer,@sprintf("He %12.5f %12.5f %12.5f\n", vt[1], vt[2], vt[3]))
    end
  end
end


"""

```
printxyz(time,x,sys,traj_file_buffer)
```

Subroutine that writes the current coordinates (or velocities) to 
output files, previously opened in respective buffers.

"""
function printxyz(time,x,sys,traj_file_buffer)
  n = length(x)
  println(traj_file_buffer,n)
  println(traj_file_buffer," time = ", time)
  for i in 1:n
    xt = image(x[i],sys.sides)
    if length(xt) == 2
      print(traj_file_buffer,@sprintf("He %12.5f %12.5f %12.5f\n", xt[1], xt[2], 0.))
    else
      print(traj_file_buffer,@sprintf("He %12.5f %12.5f %12.5f\n", xt[1], xt[2], xt[3]))
    end
  end
end

"""

```
printxyz(x,sys,file::String) 
```

Prints one set of coordinates (or velocities), given in `x` to file `file`. Overwrites
the file.  

"""
function printxyz(x::Vector,sys::System,file::String)
  buffer = open(file,"w") 
  printxyz(0,x,sys,buffer)
  close(buffer)
  return nothing
end
